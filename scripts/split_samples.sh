#!/usr/bin/bash

SUBSET_BAM="/home/dwm269/usr/bin/subset-bam_1.1.0"
BAM2FQ="/home/dwm269/usr/bin/bamtofastq_1.4.0"
NCORES=16
OUTDIR="/workdir/dwm269/totalRNA/data/fastq/Vis_yPAP3_2K/split_by_section"
TMPDIR="/workdir/dwm269/totalRNA/tmp"

WHITELIST="/home/dwm269/DWM_utils/align_pipes/10x_STARsolo/resources/barcodes_10x/visium-v1_coordinates.txt"

mkdir -p ${OUTDIR}

##############################################################################
# BAMPATH="/workdir/dwm269/totalRNA/data/STARsolo/GRCm39_GENCODEM27_Visium/Vis_yPAP_3C_STARsolo/Aligned.sortedByCoord.out.bam"
# JSONPATH="/workdir/dwm269/totalRNA/data/visium_images/Vis_yPAP3/V10T11-086-C1_D5SkM.json"
BAMPATH="/workdir/dwm269/totalRNA/data/spaceranger_count/Vis_yPAP_3C_TA-D5_spaceranger_count/outs/possorted_genome_bam.bam"

# python /home/dwm269/DWM_utils/sc_utils/TXG_helpers/json2CBs.py ${JSONPATH} ${OUTDIR}/V10T11-086-C1_D5SkM_CBs.txt ${WHITELIST}
# echo "Done extracting CBs for C1_D5SkM!"

# echo "Subsetting .bam file..."
#
# ${SUBSET_BAM} \
# --cores=${NCORES} \
# --bam=${BAMPATH} \
# --cell-barcodes=${OUTDIR}/V10T11-086-C1_D5SkM_CBs.txt \
# --out-bam=${TMPDIR}/C1_D5SkM.bam
#
# echo "Done."

# echo "Converting .bam to .fastq(s)..."
# ${BAM2FQ} \
# --reads-per-fastq=1000000000000 \
# ${BAMPATH} \
# ${OUTDIR}/C1_D5SkM

# ${TMPDIR}/C1_D5SkM.bam \
# echo "Done."

# ##############################################################################
# BAMPATH="/workdir/dwm269/totalRNA/data/spaceranger_count/Vis_yPAP_3C_TA-D7_spaceranger_count/outs/possorted_genome_bam.bam"
# ${BAM2FQ} \
# --reads-per-fastq=1000000000000 \
# ${BAMPATH} \
# ${OUTDIR}/C1_D7SkM
#
# ##############################################################################
# BAMPATH="/workdir/dwm269/totalRNA/data/spaceranger_count/Vis_yPAP_3D_TA-D0_spaceranger_count/outs/possorted_genome_bam.bam"
# ${BAM2FQ} \
# --reads-per-fastq=1000000000000 \
# ${BAMPATH} \
# ${OUTDIR}/D1_D0SkM
#
# ##############################################################################
# BAMPATH="/workdir/dwm269/totalRNA/data/spaceranger_count/Vis_yPAP_3D_TA-D2_spaceranger_count/outs/possorted_genome_bam.bam"
# ${BAM2FQ} \
# --reads-per-fastq=1000000000000 \
# ${BAMPATH} \
# ${OUTDIR}/D1_D2SkM
#
##############################################################################
##############################################################################
BAMPATH="/workdir/dwm269/totalRNA/data/spaceranger_count/Vis_yPAP_2A_TA-D0_spaceranger_count/outs/possorted_genome_bam.bam"
${BAM2FQ} \
--reads-per-fastq=1000000000000 \
${BAMPATH} \
${OUTDIR}/2A_SkMD0

##############################################################################
BAMPATH="/workdir/dwm269/totalRNA/data/spaceranger_count/Vis_yPAP_2A_TA-D5_spaceranger_count/outs/possorted_genome_bam.bam"
${BAM2FQ} \
--reads-per-fastq=1000000000000 \
${BAMPATH} \
${OUTDIR}/2A_SkMD5

##############################################################################
BAMPATH="/workdir/dwm269/totalRNA/data/spaceranger_count/Vis_yPAP_2B_TA-D0_spaceranger_count/outs/possorted_genome_bam.bam"
${BAM2FQ} \
--reads-per-fastq=1000000000000 \
${BAMPATH} \
${OUTDIR}/2B_SkMD0

##############################################################################
BAMPATH="/workdir/dwm269/totalRNA/data/spaceranger_count/Vis_yPAP_2B_TA-D2_spaceranger_count/outs/possorted_genome_bam.bam"
${BAM2FQ} \
--reads-per-fastq=1000000000000 \
${BAMPATH} \
${OUTDIR}/2B_SkMD2

##############################################################################
BAMPATH="/workdir/dwm269/totalRNA/data/spaceranger_count/Vis_yPAP_2C_Heart-T1L-D7PI_spaceranger_count/outs/possorted_genome_bam.bam"
${BAM2FQ} \
--reads-per-fastq=1000000000000 \
${BAMPATH} \
${OUTDIR}/2C_T1L

##############################################################################
BAMPATH="/workdir/dwm269/totalRNA/data/spaceranger_count/Vis_yPAP_2D_TA-D5_spaceranger_count/outs/possorted_genome_bam.bam"
${BAM2FQ} \
--reads-per-fastq=1000000000000 \
${BAMPATH} \
${OUTDIR}/2D_SkMD5

##############################################################################
BAMPATH="/workdir/dwm269/totalRNA/data/spaceranger_count/Vis_yPAP_2D_TA-D7_spaceranger_count/outs/possorted_genome_bam.bam"
${BAM2FQ} \
--reads-per-fastq=1000000000000 \
${BAMPATH} \
${OUTDIR}/2D_SkMD7
