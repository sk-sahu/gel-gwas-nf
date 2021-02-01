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
  echo true

  input:
  set val(name), val(chr), file(vcf), file(index) from vcfsCh

  output:
  set val(name), val(chr), file("${name}.masked_filtered.vcf.gz"), file("${name}.masked_filtered.vcf.gz.tbi") into maskedVcfsCh
  
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
            --threads 2 \
            -Oz -o ${name}.masked.vcf.gz
tabix ${name}.masked.vcf.gz

#rm \$(realpath ${vcf})
#rm \$(realpath ${index})

bcftools view ${name}.masked.vcf.gz -Oz -o ${name}.masked_filtered.vcf.gz --threads 2 -i 'F_MISSING<0.05'

tabix ${name}.masked_filtered.vcf.gz
du -h ${name}.masked_filtered.vcf.gz

#rm \$(realpath ${name}.masked.vcf.gz)
md5sum ${name}.masked_filtered.vcf.gz
md5sum ${name}.masked_filtered.vcf.gz.tbi
"""
}
}
