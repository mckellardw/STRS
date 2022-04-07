#!/usr/bin/bash

# Align reads for C2C12 + PAP (MeOH-fixed) nuclei samples

DATADIR=/workdir/dwm269/totalRNA/data/c2c12_yeast_pap_bulk
OUTDIR=/workdir/dwm269/totalRNA/data/c2c12_yeast_pap_bulk/kallisto_mm10_out
TMPDIR=/workdir/dwm269/tmp/kallisto1

INDEXREF=/workdir/dwm269/genomes/mm10_all/mm10_kallisto/mm10_genome.fa.idx
GTFREF=/workdir/dwm269/genomes/mm10_all/mm10_kallisto/mm10_genes.gtf
CHROMINFO=/workdir/dwm269/genomes/mm10_all/mm10_chr.chromInfo

ADAPT_SEQ_R1="AAGCAGTGGTATCAACGCAGAGTGAATGGG"

NCORES=20

mkdir $OUTDIR
mkdir $TMPDIR

# conda activate kallisto1

## Remove adapter, UMI barcode, additional barcode, and low quality (q=20) base from 3prime end of reads. Keep read length >=15 after trimmming
# Remove adapter
SAMPLEID="12340_11006_133813_JGP4K_M_10_CGTACTAG"

cutadapt -a ${ADAPT_SEQ_R1} -e 0.10 --overlap 2 -q 20 \
--output=${TMPDIR}/${SAMPLEID}_q20trim_R1.fastq \
${DATADIR}/${SAMPLEID}_R1.fastq.gz

# --untrimmed-output=${TMPDIR}/${SAMPLEID}_untrim_R1.fastq \

# remove UMI2 and ADD_B2 from the 3 prime end of R1
 #n2=$[UMI2+ADD_B2]
 #cutadapt --cut -${n2} --minimum-length=10 ${TMPDIR}/${SAMPLEID}_trim_R1.fastq --output=${TMPDIR}/${SAMPLEID}_trim.${n2}Nremoved_R1.fastq -q 20 &
 #cutadapt --minimum-length=10 ${TMPDIR}/${SAMPLEID}_untrim_R1.fastq --output=${TMPDIR}/${SAMPLEID}_q20trim_R1.fastq -q 20 &
 #wait

 # cat ${TMPDIR}/${SAMPLEID}_q20trim_R1.fastq ${TMPDIR}/${name}_trim.${n2}Nremoved_R1.fastq | paste - - - - |LC_ALL=C sort --temporary-directory=${TMPDIR} --parallel=10 -k1,1 -S 10G | tr '\t' '\n' > ${TMPDIR}/noadapt/${name}_noadapt_R1.fastq &

# align w/ kallisto
kallisto quant \
-i ${INDEXREF} \
--output-dir ${OUTDIR} \
--single -l 76 \
-s 20 \
--threads ${NCORES} \
--genomebam \
--gtf ${GTFREF} \
${TMPDIR}/${SAMPLEID}_q20trim_R1.fastq

####
SAMPLEID="12340_11006_133814_JGP4K_M_19_CGAGGCTG"

cutadapt -a ${ADAPT_SEQ_R1} -e 0.10 --overlap 2 -q 20 \
--output=${TMPDIR}/${SAMPLEID}_q20trim_R1.fastq \
${DATADIR}/${SAMPLEID}_R1.fastq.gz

# align w/ kallisto
kallisto quant \
-i ${INDEXREF} \
--output-dir ${OUTDIR} \
--single -l 76 \
-s 20 \
--threads ${NCORES} \
--genomebam \
--gtf ${GTFREF} \
${TMPDIR}/${SAMPLEID}_q20trim_R1.fastq

####
SAMPLEID="12340_11006_134756_JGP4K_M_11_AGGCAGAA"

cutadapt -a ${ADAPT_SEQ_R1} -e 0.10 --overlap 2 -q 20 \
--output=${TMPDIR}/${SAMPLEID}_q20trim_R1.fastq \
${DATADIR}/${SAMPLEID}_R1.fastq.gz

# align w/ kallisto
kallisto quant \
-i ${INDEXREF} \
--output-dir ${OUTDIR} \
--single -l 76 \
-s 20 \
--threads ${NCORES} \
--genomebam \
--gtf ${GTFREF} \
${TMPDIR}/${SAMPLEID}_q20trim_R1.fastq

####
SAMPLEID="12340_11006_133815_JGP4K_M_20_AAGAGGCA"

cutadapt -a ${ADAPT_SEQ_R1} -e 0.10 --overlap 2 -q 20 \
--output=${TMPDIR}/${SAMPLEID}_q20trim_R1.fastq \
${DATADIR}/${SAMPLEID}_R1.fastq.gz

# align w/ kallisto
kallisto quant \
-i ${INDEXREF} \
--output-dir ${OUTDIR} \
--single -l 76 \
-s 20 \
--threads ${NCORES} \
--genomebam \
--gtf ${GTFREF} \
${TMPDIR}/${SAMPLEID}_q20trim_R1.fastq
