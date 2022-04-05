# **Resources**
Gene lists, accession IDs, and other useful info for total RNA analyses

### cb_lists
Lists of spot barcodes for each sample assayed in the study. Some samples were assayed by placing 2 tissue slices in a single capture area - these barcodes make it much easier to select spots from each sample

### DGEA
Differential Gene Expression Analysis outputs

### gene_lists
- Assorted lists of genes and gene metadata
- Command line code snippet for processing GENCODE M28 .gtf file into a parse-able table:
```
awk 'BEGIN{FS="\t"}{split($9,a,";"); if($3~"gene") print a[1]"\t"a[3]"\t"$1":"$4"-"$5"\t"a[2]"\t"$7}' <path to gtf file>/GRCm39_ReoT1L_merged_genes.gtf \
| sed 's/gene_id "//' | sed 's/gene_id "//' | sed 's/gene_type "//'| sed 's/gene_name "//' | sed 's/gene_type "//' | sed 's/"//g' | sed 's/ //g' \
| sed '1iGENEID\tGeneSymbol\tChromosome\tBiotype\tStrand' > resources/gene_lists/GRCm39_GENCODEm28_gene_info_gtf.tsv
```

### metadata_sheets
Sample sheets for samples included in this study


## ***Other resources not used in our study...***
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

### Finding putative cis-NATs
Get positive strand transcripts
```
awk '$3 == "transcript"' GENCODE_M28/gencode.vM28.chr_patch_hapl_scaff.annotation.gtf | awk '$7 == "+"'  > GENCODE_M28/gencode.vM28.chr_patch_hapl_scaff.annotation_PLUS.gtf
```
Get negative strand transcripts
```
awk '$3 == "transcript"' GENCODE_M28/gencode.vM28.chr_patch_hapl_scaff.annotation.gtf | awk '$7 == "-"' > GENCODE_M28/gencode.vM28.chr_patch_hapl_scaff.annotation_MINUS.gtf
```
Find the intersection
```
bedtools intersect -S -wa -wb -a GENCODE_M28/gencode.vM28.chr_patch_hapl_scaff.annotation_PLUS.gtf -b GENCODE_M28/gencode.vM28.chr_patch_hapl_scaff.annotation_MINUS.gtf \
| awk 'BEGIN{FS="\t"}{split($9,a,";"); split($18,b,";"); print a[4]"\t"a[1]"\t"a[3]"\t"$1":"$4"-"$5"\t"$7"\t"b[4]"\t"b[1]"\t"b[3]"\t"$10":"$13"-"$14"\t"$16}' \
| sed 's/gene_id "//' | sed 's/gene_id "//' | sed 's/gene_type "//' | sed 's/gene_name "//' | sed 's/gene_name//' | sed 's/gene_type "//' | sed 's/"//g' | sed 's/ //g' \
| sort -u -k 4,4 \
| sed '1iGeneSymbol_Plus\tGENEID_Plus\tBiotype_Plus\tLocation_Plus\tStrand_Plus\tGeneSymbol_Minus\tGENEID_Minus\tBiotype_Minus\tLocation_Minus\tStrand_Minus' > /workdir/dwm269/totalRNA/spTotal/resources/gene_lists/overlap_M28.tsv
```
