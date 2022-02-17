# spTotal Resources
Gene lists, accession IDs, and other useful info for total RNA analyses

## gene_lists
- Assorted lists of genes and gene metadata
- Command line code snippet for processing GENCODE M28 .gtf file into a parse-able table:
```
awk 'BEGIN{FS="\t"}{split($9,a,";"); if($3~"gene") print a[1]"\t"a[3]"\t"$1":"$4"-"$5"\t"a[2]"\t"$7}' <path to gtf file>/GRCm39_ReoT1L_merged_genes.gtf \
| sed 's/gene_id "//' | sed 's/gene_id "//' | sed 's/gene_type "//'| sed 's/gene_name "//' | sed 's/gene_type "//' | sed 's/"//g' | sed 's/ //g' \
| sed '1iGENEID\tGeneSymbol\tChromosome\tBiotype\tStrand' > resources/gene_lists/GRCm39_GENCODEm28_gene_info_gtf.tsv
```

### TargetScan8
- Contains miRNA target predictions from the [TargetScan8.0 website](http://www.targetscan.org/mmu_80/)
- Files are named according to the miRbase annotation for the miRNA they correspond to (default names from the website)
- R code snippet to load tables:
```
ts8.path <- "/path/to/spTotal/resources/TargetScan8/"
ts8_targets <- list(
  let7_5p = read.csv(paste0(ts8.path,"TargetScan8.0__let-7-5p_miR-98-5p.predicted_targets.txt"),sep = "\t",header = T),
  mir1_3p = read.csv(paste0(ts8.path,"TargetScan8.0__miR-1-3p_206-3p.predicted_targets.txt"),sep = "\t",header = T),
  mir133_3p = read.csv(paste0(ts8.path,"TargetScan8.0__miR-133-3p.predicted_targets.txt"),sep = "\t",header = T)
) %>%
  lapply( # Filter target lists based on number of binding sites
    FUN=function(X) X[X[["Aggregate.PCT"]]>0.9,]
  )
```
