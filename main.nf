#!/usr/bin/env nextflow
/*
========================================================================================
                         lifebit-ai/gel-gwas
========================================================================================
 lifebit-ai/gel-gwas GWAS pipeline built for Genomics England using SAIGE
 #### Homepage / Documentation
 https://github.com/lifebit-ai/gel-gwas
----------------------------------------------------------------------------------------
*/

/*--------------------------------------------------
  Channel setup
---------------------------------------------------*/
Channel
  .fromPath(params.phenoFile)
  .ifEmpty { exit 1, "Pheno file not found: ${params.phenoFile}" }
  .into { phenoCh; phenoCh_gwas_filtering}
Channel
  .fromPath(params.sampleFile)
  .ifEmpty { exit 1, "Sample file not found: ${params.sampleFile}" }
  .set { sampleCh }
Channel
  .fromFilePairs("${params.plinkFile}",size:3, flat : true)
  .ifEmpty { exit 1, "PLINK files not found: ${params.plinkFile}" }
  .set { plinkCh }
Channel
  .fromPath(params.plink_keep_pheno)
  .set {plink_keep_pheno_ch}
Channel
  .fromPath(params.vcfsList)
  .ifEmpty { exit 1, "Cannot find CSV VCFs file : ${params.vcfsList}" }
  .splitCsv(skip:1)
  .map { chr, vcf, index -> [file(vcf).simpleName, chr, file(vcf), file(index)] }
  .set { vcfsCh }
Channel
  .fromPath(params.gwas_cat)
  .ifEmpty { exit 1, "Cannot find GWAS catalogue CSV  file : ${params.gwas_cat}" }
  .set { ch_gwas_cat }
/*--------------------------------------------------
  Pre-GWAS masking - download and mask vcfs
---------------------------------------------------*/
if ( params.skip_masking) { maskedVcfsCh = vcfsCh }
if (!params.skip_masking) {
process gwas_masking {
  tag "$name"
  publishDir "${params.outdir}/gwas_masking", mode: 'copy'
  
  input:
  set val(name), val(chr), file(vcf), file(index) from vcfsCh

  output:
  set val(name), val(chr), file("${name}.masked_filtered.vcf.gz"), file("${name}.masked_filtered.vcf.gz.csi") into maskedVcfsCh
  
  script:
 """ 
  bcftools +setGT ${name}.vcf.gz -Ou -- \
            -t q \
            -i \"FMT/DP<10 | FMT/GQ<20\" \
            -n . \
        | bcftools +setGT -Ou -- \
            -t \"b:AD<=0.001\" \
            -n . \
        | bcftools view \
            -Oz -o ${name}.masked.vcf.gz
tabix ${name}.masked.vcf.gz

bcftools view ${name}.masked.vcf.gz -Oz -o ${name}.masked_filtered.vcf.gz \
            -i 'F_MISSING<0.05'
bcftools index ${name}.masked_filtered.vcf.gz
"""
}
}
/*--------------------------------------------------
  Pre-GWAS filtering - download, filter and convert VCFs
---------------------------------------------------*/

if ( params.skip_gwas_filtering) { filteredVcfsCh = maskedVcfsCh }
if (!params.skip_gwas_filtering) {
process gwas_filtering {
  tag "$name"
  publishDir "${params.outdir}/gwas_filtering", mode: 'copy'

  input:
  set val(name), val(chr), file(vcf), file(index) from maskedVcfsCh
  each file(phenofile) from phenoCh_gwas_filtering
  each file(plink_keep_file) from plink_keep_pheno_ch

  output:
  set val(name), val(chr), file("${name}.filtered_final.vcf.gz"), file("${name}.filtered_final.vcf.gz.csi") into filteredVcfsCh
  set val(name), val(chr), file("${name}.filtered_final.bgen") into filteredbgenCh
  
  file("${name}_filtered.{bed,bim,fam}") into plinkTestCh

  script:
  // TODO: (High priority) Only extract needed individuals from VCF files with `bcftools -S samples.txt` - get from samples file?
  // TODO: (Not required) `bcftools -T sites_to_extract.txt`
  // Optional parameters
  extra_plink_filter_missingness_options = params.plink_keep_pheno != "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/nodata" ? "--keep ${plink_keep_file}" : ""
  """
  # Download, filter and convert (bcf or vcf.gz) -> vcf.gz
  bcftools view -q ${params.qFilter} $vcf -Oz -o ${name}_filtered.vcf.gz
  bcftools index ${name}_filtered.vcf.gz

  # Create PLINK binary from vcf.gz
  plink2 \
    --make-bed \
    --set-missing-var-ids '${params.plink_set_missing_var_ids}' \
    --vcf ${name}_filtered.vcf.gz \
    --out ${name}_filtered \
    --vcf-half-call ${params.plink_vcf_half_call} \
    --double-id \
    --set-hh-missing \
    --new-id-max-allele-len ${params.plink_new_id_max_allele_len} missing \
    --output-chr  ${params.plink_output_chr}
  
  #Filter missingness
  plink2 \
    --bfile ${name}_filtered \
    --pheno $phenofile \
    --pheno-name ${params.phenoCol} \
    --allow-no-sex \
    --test-missing midp \
    --out ${name} \
    --1 \
    --keep-allele-order \
    ${extra_plink_filter_missingness_options} \
    --output-chr ${params.plink_output_chr}

  awk -v thresm=${params.thres_m} '\$5 < thresm {print}'  ${name}.missing > ${name}.missing_FAIL

  #Filter HWE
  plink2 \
    --bfile ${name}_filtered \
    --pheno $phenofile \
    --pheno-name ${params.phenoCol} \
    --allow-no-sex \
    --hwe ${params.thres_HWE} midp \
    --out ${name}.filtered_final \
    --make-just-bim \
    --exclude ${name}.missing_FAIL \
    --1 \
    --keep-allele-order \
    ${extra_plink_filter_missingness_options} \
    --output-chr ${params.plink_output_chr} \
    --export bgen-1.2 bits=8 ref-first
    
  bcftools view ${name}_filtered.vcf.gz | awk -F '\\t' 'NR==FNR{c[\$1\$4\$6\$5]++;next}; c[\$1\$2\$4\$5] > 0' ${name}.filtered_final.bim - | bgzip > ${name}.filtered_temp.vcf.gz
  bcftools view -h ${name}_filtered.vcf.gz -Oz -o ${name}_filtered.header.vcf.gz
  cat ${name}_filtered.header.vcf.gz ${name}.filtered_temp.vcf.gz > ${name}.filtered_final.vcf.gz
  bcftools index ${name}.filtered_final.vcf.gz
  
  """
}
}

/*--------------------------------------------------
  GWAS Analysis 1 with SAIGE - Fit the null mixed-model
---------------------------------------------------*/

process gwas_1_fit_null_glmm {
  tag "$plink_GRM_snps"
  publishDir "${params.outdir}/gwas_1_fit_null_glmm", mode: 'copy'

  input:
  set val(plink_GRM_snps), file(bed), file(bim), file(fam) from plinkCh
  each file(phenoFile) from phenoCh

  output:
  file "*" into fit_null_glmm_results
  file ("step1_${params.phenoCol}_out.rda") into rdaCh
  file ("step1_${params.phenoCol}.varianceRatio.txt") into varianceRatioCh

  script:
  """
  step1_fitNULLGLMM.R \
    --plinkFile=${plink_GRM_snps} \
    --phenoFile=${phenoFile} \
    --phenoCol=${params.phenoCol} \
    --sampleIDColinphenoFile=IID \
    --traitType=${params.traitType} \
    --outputPrefix=step1_${params.phenoCol}_out \
    --outputPrefix_varRatio=step1_${params.phenoCol} \
    --nThreads=${task.cpus} ${params.saigeStep1ExtraFlags}
  """
}

/*--------------------------------------------------
  GWAS Analysis 2 with SAIGE - Perform mixed-model association testing
---------------------------------------------------*/

if (params.gwas_2_spa_tests_format_switch) {
process gwas_2_spa_tests_vcf {
  tag "$name"
  publishDir "${params.outdir}/gwas_2_spa_tests", mode: 'copy'

  input:
  set val(name), val(chr), file(vcf), file(index) from filteredVcfsCh
  each file(rda) from rdaCh
  each file(varianceRatio) from varianceRatioCh

  output:
  file "*" into results
  file("*.SAIGE.gwas.txt") into ch_report

  script:
  """
  step2_SPAtests.R \
    --vcfFile=${vcf} \
    --vcfFileIndex=${index} \
    --vcfField=GT \
    --chrom=${chr} \
    --minMAC=20 \
    --GMMATmodelFile=${rda} \
    --varianceRatioFile=${varianceRatio} \
    --SAIGEOutputFile=${params.phenoCol}.${name}.SAIGE.gwas.txt \
    --numLinesOutput=2 \
    --IsOutputAFinCaseCtrl=TRUE \
    --IsDropMissingDosages=TRUE \
    --IsOutputNinCaseCtrl=TRUE \
    --IsOutputHetHomCountsinCaseCtrl=TRUE
  """
}
}
if (!params.gwas_2_spa_tests_format_switch) {
process gwas_2_spa_tests_bgen {
  tag "$name"
  publishDir "${params.outdir}/gwas_2_spa_tests", mode: 'copy'

  input:
  set val(name), val(chr), file(vcf), file(index) from filteredVcfsCh
  each file(rda) from rdaCh
  each file(varianceRatio) from varianceRatioCh
  each file(sampleFile) from sampleCh

  output:
  file "*" into results
  file("*.SAIGE.gwas.txt") into ch_report

  script:
  """
  step2_SPAtests.R \
    --bgenFile=${vcf} \
    --bgenFileIndex=${index} \
    --chrom=${chr} \
    --minMAC=20 \
    --sampleFile=${sampleFile} \
    --GMMATmodelFile=${rda} \
    --varianceRatioFile=${varianceRatio} \
    --SAIGEOutputFile=${params.phenoCol}.${name}.SAIGE.gwas.txt \
    --numLinesOutput=2 \
    --IsOutputAFinCaseCtrl=TRUE \
    --IsDropMissingDosages=TRUE \
    --IsOutputNinCaseCtrl=TRUE \
    --IsOutputHetHomCountsinCaseCtrl=TRUE
  """
}
}
/*--------------------------------------------------
  GWAS Analysis 2 with SAIGE - Generate report
---------------------------------------------------*/
if (!params.skip_report) {
  process create_report {
    tag "report"
    publishDir "${params.outdir}/MultiQC/", mode: 'copy'

    input:
    file(saige_output) from ch_report.collect()
    file(gwas_cat) from ch_gwas_cat

    output:
    file "multiqc_report.html" into ch_report_outputs
    set file("*png"), file("*ipynb"), file("*csv") into ch_report_outputs_all

    script:

    """
    cp /opt/bin/* .
    
    # creates 2 .csv files, saige_results_<params.output_tag>.csv, saige_results_top_n.csv
    concat_chroms.R \
        --saige_output_name='saige_results' \
        --filename_pattern='${params.saige_filename_pattern}' \
        --output_tag='${params.output_tag}' \
        --top_n_sites=${params.top_n_sites} \
        --max_top_n_sites=${params.max_top_n_sites}

    # creates gwascat_subset.csv
    subset_gwascat.R \
        --saige_output='saige_results_${params.output_tag}.csv' \
        --gwas_cat='${gwas_cat}'

    # creates <params.output_tag>_manhattan.png with analysis.csv as input
    manhattan.R \
        --saige_output='saige_results_${params.output_tag}.csv' \
        --output_tag='${params.output_tag}'

    # creates <params.output_tag>_qqplot_ci.png with analysis.csv as input
    qqplot.R \
        --saige_output='saige_results_${params.output_tag}.csv' \
        --output_tag='${params.output_tag}'

    # Generates the report
    Rscript -e "rmarkdown::render('gwas_report.Rmd', params = list(manhattan='${params.output_tag}_manhattan.png',qqplot='${params.output_tag}_qqplot_ci.png', gwascat='gwascat_subset.csv', saige_results='saige_results_top_n.csv'))"
    mv gwas_report.html multiqc_report.html

    # Generates the ipynb
    jupytext --to ipynb gwas_report.Rmd
    """
  }
}
