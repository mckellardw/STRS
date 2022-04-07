# TODO:
- Add Visium settings
- Write README
- Add custom .gtf filtering for intragenic features

#README TODO:
- Required packages & dependencies (add installation via .yml file)
- Write out pipeline details
- Info on sample_sheet format!


Sources: (#TODO)
- kallisto/BUStools (`kb`)
- https://github.com/torognes/vsearch

Output for {sample}:
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
