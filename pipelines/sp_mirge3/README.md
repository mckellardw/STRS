# **sp_mirge3**
Snakemake pipeline to run [miRge3.0](https://github.com/mhalushka/miRge3.0) (gold standard small RNA alignment) on either spatial or single-cell RNAseq data

### Inputs
- Aligned/barcode-tagged .bam file (written for files aligned by STAR)
- Cell barcode list

### Outputs
- See miRge3.0 outputs [here](https://github.com/mhalushka/miRge3.0)

### Dependencies
- [samtools](http://www.htslib.org/)
- [miRge3.0](https://github.com/mhalushka/miRge3.0)
