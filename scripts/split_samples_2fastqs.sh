#!/usr/bin/bash

DATADIR="/workdir/dwm269/totalRNA/data/fastq"
TMPDIR="/workdir/dwm269/totalRNA/tmp"
CBs_DIR="/workdir/dwm269/totalRNA/spTotal/resources/cb_lists"

WHITELIST="/home/dwm269/DWM_utils/align_pipes/10x_STARsolo/resources/barcodes_10x/visium-v1_coordinates.txt"
NCORES=16

##############################################################################
BAMPATH="/workdir/dwm269/totalRNA/data/STARsolo/GRCm39_GENCODEM28_Visium/Vis_yPAP_3C_STARsolo/Aligned.sortedByCoord.dedup.out.bam"

echo "Done extracting CBs for C1_D5SkM!"
samtools view -1 -b \
-@ ${NCORES} \
--tag-file CB:${CBs_DIR}/yPAP-Pro_SkM-D5.txt \
${BAMPATH} \
| samtools fastq \
> ${DATADIR}/split_by_section/3C_D5SkM/split_dedup_R2.fastq
pigz -p ${NCORES} ${DATADIR}/split_by_section/3C_D5SkM/split_dedup_R2.fastq

echo "Done extracting CBs for C1_D7SkM!"
samtools view -1 -b \
-@ ${NCORES} \
--tag-file CB:${CBs_DIR}/yPAP-Pro_SkM-D7.txt \
${BAMPATH} \
| samtools fastq \
-2 ${DATADIR}/split_by_section/3C_D7SkM/split_dedup_R2.fastq \
> ${DATADIR}/split_by_section/3C_D7SkM/split_dedup_R2.fastq
pigz -p ${NCORES} ${DATADIR}/split_by_section/3C_D7SkM/split_dedup_R2.fastq

##############################################################################
BAMPATH="/workdir/dwm269/totalRNA/data/STARsolo/GRCm39_GENCODEM28_Visium/Vis_yPAP_3D_STARsolo/Aligned.sortedByCoord.dedup.out.bam"

echo "Done extracting CBs for D1_D0SkM!"
samtools view -1 -b \
-@ ${NCORES} \
--tag-file CB:${CBs_DIR}/yPAP-Pro_SkM-D0.txt \
${BAMPATH} \
| samtools fastq \
-2 ${DATADIR}/split_by_section/3D_D0SkM/split_dedup_R2.fastq \
> ${DATADIR}/split_by_section/3D_D0SkM/split_dedup_R2.fastq
pigz -p ${NCORES} ${DATADIR}/split_by_section/3D_D0SkM/split_dedup_R2.fastq

echo "Done extracting CBs for C1_D7SkM!"
samtools view -1 -b \
-@ ${NCORES} \
--tag-file CB:${CBs_DIR}/yPAP-Pro_SkM-D5.txt \
${BAMPATH} \
| samtools fastq \
-2 ${DATADIR}/split_by_section/3D_D0SkM/split_dedup_R2.fastq \
> ${DATADIR}/split_by_section/3D_D0SkM/split_dedup_R2.fastq
pigz -p ${NCORES} ${DATADIR}/split_by_section/3D_D0SkM/split_dedup_R2.fastq
