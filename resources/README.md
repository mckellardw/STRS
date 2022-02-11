# spTotal Resources
Gene lists, accession IDs, and other useful info for total RNA analyses

## gene_lists
- Assorted lists of genes and gene metadata
- Command line code snippet for processing GENCODE M28 .gtf file:
```
awk 'BEGIN{FS="\t"}{split($9,a,";"); if($3~"gene") print a[1]"\t"a[3]"\t"$1":"$4"-"$5"\t"a[2]"\t"$7}' GRCm39_ReoT1L_merged_genes.gtf | sed 's/gene_id "//' | sed 's/gene_id "//' | sed 's/gene_type "//'| sed 's/gene_name "//' | sed 's/gene_type "//' | sed 's/"//g' | sed 's/ //g' | sed '1iGENEID\tGeneSymbol\tChromosome\tBiotype\tStrand'
```

### TargetScan8
- Contains miRNA target predictions from the [TargetScan8.0 website](http://www.targetscan.org/mmu_80/)
- Files are named according to the miRbase annotation for the miRNA they correspond to (default names from the website)
