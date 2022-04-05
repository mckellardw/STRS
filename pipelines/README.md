# **Pipelines**

Snakemake workflows used to align, quantify, and quality check sequencing libraries from STRS data. Brief descriptions below, see each pipeline for detailed descriptions on input/output files and package dependencies.

## 10x_kallisto
Read preprocessing, QC, and pseudoalignment/quantification with [kallisto and bustools](https://github.com/pachterlab/kallistobustools)

## 10x_STAR
Read preprocessing, QC, and alignment/quantification with [STARsolo](https://github.com/alexdobin/STAR)

## sp_mirge3
After alignment (w/ STAR), this pipeline splits the barcoded and UMI-deduplicated reads across spots/cells, then runs [miRge3.0]() on each spot individually to align/quantify small RNAs.
