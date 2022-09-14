# 10x_STARsolo
Preprocessing, alignment, QC, and quantification workflow for 10x Genomics data (Chromium & Visium)
**David W. McKellar**

#README TODO:
- Required packages & dependencies (add installation via .yml file)
- Write out pipeline details
- Info on sample_sheet format

Dependencies & Sources:
- `cutadapt` [v3.4](https://cutadapt.readthedocs.io/en/stable/)
- `fastqc` [vv0.11.8](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- `STAR` [v2.7.10a](https://github.com/alexdobin/STAR)
- `qualimap` [v2.2.2a](http://qualimap.conesalab.org/)
- `vsearch` [v2.17.0_linux_x86_64](https://github.com/torognes/vsearch)
- `BLAST`

Outputs:
```
SAMPLE_ID
├── Aligned.sortedByCoord.dedup.out.bam
├── Aligned.sortedByCoord.dedup.out.bam.bai
├── Aligned.sortedByCoord.dedup.out.bed.gz
├── Aligned.sortedByCoord.dedup.out_log10_minus.bw
├── Aligned.sortedByCoord.dedup.out_log10_plus.bw
├── Aligned.sortedByCoord.dedup.out_merged.bw
├── Aligned.sortedByCoord.dedup.out_merged_log10.bw
├── Aligned.sortedByCoord.dedup.out_minus.bw
├── Aligned.sortedByCoord.dedup.out_plus.bw
├── Aligned.sortedByCoord.out.bam
├── Aligned.sortedByCoord.out.bam.bai
├── bam2splitBigWig.kill.warnings
├── cutadapt_polyA_report.txt
├── cutadapt_polyG_report.txt
├── Log.final.out
├── Log.out
├── Log.progress.out
├── mirbase
├── postTrim_FastQC_R2_lengthFiltered
│   ├── SAMPLE_ID_R2_final_filtered_fastqc.html
│   └── SAMPLE_ID_R2_final_filtered_fastqc.zip
├── postTrim_fastqc_R2_out
│   ├── SAMPLE_ID_R2_final_fastqc.html
│   └── SAMPLE_ID_R2_final_fastqc.zip
├── preTrim_fastqc_R2_out
│   ├── SAMPLE_ID_R2_fastqc.html
│   └── SAMPLE_ID_R2_fastqc.zip
├── qualimap_out
│   ├── css
│   │   ├── agogo.css
│   │   ├── ajax-loader.gif
│   │   ├── basic.css
│   │   ├── bgfooter.png
│   │   ├── bgtop.png
│   │   ├── comment-bright.png
│   │   ├── comment-close.png
│   │   ├── comment.png
│   │   ├── doctools.js
│   │   ├── down.png
│   │   ├── down-pressed.png
│   │   ├── file.png
│   │   ├── jquery.js
│   │   ├── minus.png
│   │   ├── plus.png
│   │   ├── pygments.css
│   │   ├── qualimap_logo_small.png
│   │   ├── report.css
│   │   ├── searchtools.js
│   │   ├── underscore.js
│   │   ├── up.png
│   │   ├── up-pressed.png
│   │   └── websupport.js
│   ├── images_qualimapReport
│   │   ├── Coverage Profile Along Genes (High).png
│   │   ├── Coverage Profile Along Genes (Low).png
│   │   ├── Coverage Profile Along Genes (Total).png
│   │   ├── Junction Analysis.png
│   │   ├── Reads Genomic Origin.png
│   │   └── Transcript coverage histogram.png
│   ├── qualimapReport.html
│   ├── raw_data_qualimapReport
│   │   ├── coverage_profile_along_genes_(high).txt
│   │   ├── coverage_profile_along_genes_(low).txt
│   │   └── coverage_profile_along_genes_(total).txt
│   └── rnaseq_qc_results.txt
├── SJ.out.tab
├── Solo.out
│   ├── Barcodes.stats
│   ├── Gene
│   │   ├── Features.stats
│   │   ├── filtered
│   │   │   ├── barcodes.tsv.gz
│   │   │   ├── features.tsv.gz
│   │   │   └── matrix.mtx.gz
│   │   ├── raw
│   │   │   ├── barcodes.tsv.gz
│   │   │   ├── features.tsv.gz
│   │   │   └── matrix.mtx.gz
│   │   ├── Summary.csv
│   │   └── UMIperCellSorted.txt
│   ├── GeneFull
│   │   ├── Features.stats
│   │   ├── filtered
│   │   │   ├── barcodes.tsv.gz
│   │   │   ├── features.tsv.gz
│   │   │   └── matrix.mtx.gz
│   │   ├── raw
│   │   │   ├── barcodes.tsv.gz
│   │   │   ├── features.tsv.gz
│   │   │   └── matrix.mtx.gz
│   │   ├── Summary.csv
│   │   └── UMIperCellSorted.txt
│   ├── SJ
│   │   ├── Features.stats
│   │   ├── raw
│   │   │   ├── barcodes.tsv.gz
│   │   │   ├── features.tsv.gz
│   │   │   └── matrix.mtx.gz
│   │   └── Summary.csv
│   └── Velocyto
│       ├── Features.stats
│       ├── filtered
│       │   ├── ambiguous.mtx.gz
│       │   ├── barcodes.tsv.gz
│       │   ├── features.tsv.gz
│       │   ├── spliced.mtx.gz
│       │   └── unspliced.mtx.gz
│       ├── raw
│       │   ├── ambiguous.mtx.gz
│       │   ├── barcodes.tsv.gz
│       │   ├── features.tsv.gz
│       │   ├── spliced.mtx.gz
│       │   └── unspliced.mtx.gz
│       └── Summary.csv
├── tmp
│   ├── SAMPLE_ID_R1_final_filtered.fq.gz
│   ├── SAMPLE_ID_R1_final.fq.gz
│   ├── SAMPLE_ID_R2_final_filtered.fq.gz
│   ├── SAMPLE_ID_R2_final.fq.gz
├── tmp.bam
├── tmp.bam.bai
├── umitools_dedup
│   ├── dedup.log
│   ├── SAMPLE_ID_edit_distance.tsv
│   ├── SAMPLE_ID_per_umi_per_position.tsv
│   └── SAMPLE_ID_per_umi.tsv
├── Unmapped_fastqc_out
│   ├── Unmapped.out.mate1_fastqc.html
│   ├── Unmapped.out.mate1_fastqc.zip
│   ├── Unmapped.out.mate2_fastqc.html
│   └── Unmapped.out.mate2_fastqc.zip
├── Unmapped.out.mate1.fastq.gz
├── Unmapped.out.mate2_blastResults.txt
└── Unmapped.out.mate2.fastq.gz

```
