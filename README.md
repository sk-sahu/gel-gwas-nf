# GEL GWAS
To be used for local tests. 
## Example usage
```bash
nextflow run ~/gel_gwas_nf2/gel-gwas-nf/main.nf \
  --plinkFile "sampleA.{bed,bim,fam}" \
  --phenoFile "pheno.txt" \
  --phenoCol pheno \
  --vcfsList "vcfs_list.csv" \
  --skip_gwas_filtering false \
  --skip_masking false \
  --gwas_cat "gwascat.csv" \
  --sampleFile "samples.txt" \
  --skip_report false \
  --plink_output_chr 26 \
  -with-dag flowchart.png
```
