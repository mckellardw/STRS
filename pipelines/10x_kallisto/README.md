# `10x_kallisto`
Snakemake pipeline for preprocessing, pseudoaligning, & quantifying single-cell and spatial RNAseq libraries w/ kallisto and BUStools

# Sources:
- kallisto/BUStools (`kb`) [link][https://www.kallistobus.tools/]
- [kallisto manual](http://pachterlab.github.io/kallisto/manual.html)
- [vsearch](https://github.com/torognes/vsearch)

# Barcode whitelist info:
[Name of barcodes list	Chemistry:](https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist-)
`3M-febrary-2018.txt.gz`: Single Cell 3' v3, Single Cell 3' v3.1, Single Cell 3' HT v3.1
`737k-august-2016.txt`:	Single Cell 3' v2, Single Cell 5' v1 and v2, Single Cell 5' HT v2
`737k-april-2014_rc.txt`:	Single Cell 3' v1
`737k-arc-v1.txt.gz`:	Single Cell Multiome (ATAC+GEX) v1
`9K-LT-march-2021.txt.gz`:	Single Cell 3' LT
`737k-fixed-rna-profiling.txt.gz`:	Fixed RNA Profiling (Present starting from Cell Ranger v7.0)

# Output for {sample}:
*Note*- comment out undesired outputs in `Snakefile` to reduce run time
```
{sample}
├── cutadapt_polyA_report.txt
├── cutadapt_polyG_report.txt
├── **kb_lamanno**
│   ├── **counts_filtered**
│   │   ├── adata.h5ad
│   │   ├── spliced.barcodes.txt
│   │   ├── spliced.genes.txt
│   │   ├── spliced.mtx
│   │   ├── unspliced.barcodes.txt
│   │   ├── unspliced.genes.txt
│   │   └── unspliced.mtx
│   ├── **counts_unfiltered**
│   │   ├── adata.h5ad
│   │   ├── spliced.barcodes.txt
│   │   ├── spliced.genes.txt
│   │   ├── spliced.mtx
│   │   ├── unspliced.barcodes.txt
│   │   ├── unspliced.genes.txt
│   │   └── unspliced.mtx
│   ├── filter_barcodes.txt
│   ├── inspect.json
│   ├── inspect.spliced.json
│   ├── inspect.unspliced.json
│   ├── kallisto_lamanno.log
│   ├── kb_info.json
│   ├── matrix.ec
│   ├── output.bus
│   ├── output.filtered.bus
│   ├── output.unfiltered.bus
│   ├── run_info.json
│   ├── spliced.filtered.bus
│   ├── spliced.unfiltered.bus
│   ├── transcripts.txt
│   ├── unspliced.filtered.bus
│   └── unspliced.unfiltered.bus
├── **kb_standard**
│   ├── **counts_filtered**
│   │   ├── adata.h5ad
│   │   ├── cells_x_genes.barcodes.txt
│   │   ├── cells_x_genes.genes.txt
│   │   └── cells_x_genes.mtx
│   ├── **counts_unfiltered**
│   │   ├── adata.h5ad
│   │   ├── adata.h5seurat
│   │   ├── cells_x_genes.barcodes.txt
│   │   ├── cells_x_genes.genes.txt
│   │   └── cells_x_genes.mtx
│   ├── filter_barcodes.txt
│   ├── inspect.json
│   ├── kallisto_standard.log
│   ├── kb_info.json
│   ├── matrix.ec
│   ├── output.bus
│   ├── output.filtered.bus
│   ├── output.unfiltered.bus
│   ├── run_info.json
│   └── transcripts.txt
├── **postTrim_fastqc_R2_out**
│   ├── {sample}_R2_final_fastqc.html
│   └── {sample}_R2_final_fastqc.zip
└── **preTrim_fastqc_R2_out**
    ├── {sample}_R2_fastqc.html
    └── {sample}_R2_fastqc.zip
```
