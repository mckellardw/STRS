#!/usr/bin/bash

DATADIR="/workdir/dwm269/totalRNA/data/STARsolo/GRCm39_GENCODEM28_Visium"
TMPDIR="/workdir/dwm269/totalRNA/tmp"
CBs_DIR="/workdir/dwm269/totalRNA/spTotal/resources/cb_lists"

WHITELIST="/home/dwm269/DWM_utils/align_pipes/10x_STARsolo/resources/barcodes_10x/visium-v1_coordinates.txt"
NCORES=16

##############################################################################
BAMPATH="/workdir/dwm269/totalRNA/data/STARsolo/GRCm39_GENCODEM28_Visium/Vis_yPAP_3C_STARsolo/Aligned.sortedByCoord.dedup.out.bam"

# JSONPATH="/workdir/dwm269/totalRNA/data/visium_images/Vis_yPAP3/V10T11-086-C1_D5SkM.json"
# python /home/dwm269/DWM_utils/sc_utils/TXG_helpers/json2CBs.py ${JSONPATH} ${DATADIR}/V10T11-086-C1_D5SkM_CBs.txt ${WHITELIST}
echo "Done extracting CBs for C1_D5SkM!"

samtools view -1 -b \
-@ ${NCORES} \
--tag-file CB:${CBs_DIR}/yPAP-Pro_SkM-D5.txt \
${BAMPATH} \
> ${DATADIR}/Vis_yPAP_3C_STARsolo/D5SkM.Aligned.sortedByCoord.dedup.out.bam

samtools index -@ ${NCORES} ${DATADIR}/Vis_yPAP_3C_STARsolo/D5SkM.Aligned.sortedByCoord.dedup.out.bam


# JSONPATH="/workdir/dwm269/totalRNA/data/visium_images/Vis_yPAP3/V10T11-086-C1_D7SkM.json"
# python /home/dwm269/DWM_utils/sc_utils/TXG_helpers/json2CBs.py ${JSONPATH} ${DATADIR}/V10T11-086-C1_D7SkM_CBs.txt ${WHITELIST}
echo "Done extracting CBs for C1_D7SkM!"

samtools view -1 -b \
-@ ${NCORES} \
--tag-file CB:${CBs_DIR}/yPAP-Pro_SkM-D7.txt \
${BAMPATH} \
> ${DATADIR}/Vis_yPAP_3C_STARsolo/D7SkM.Aligned.sortedByCoord.dedup.out.bam

samtools index -@ ${NCORES} ${DATADIR}/Vis_yPAP_3C_STARsolo/D7SkM.Aligned.sortedByCoord.dedup.out.bam

##############################################################################
BAMPATH="/workdir/dwm269/totalRNA/data/STARsolo/GRCm39_GENCODEM28_Visium/Vis_yPAP_3D_STARsolo/Aligned.sortedByCoord.dedup.out.bam"

# JSONPATH="/workdir/dwm269/totalRNA/data/visium_images/Vis_yPAP3/V10T11-086-D1_D0SkM.json"
# python /home/dwm269/DWM_utils/sc_utils/TXG_helpers/json2CBs.py ${JSONPATH} ${DATADIR}/V10T11-086-D1_D0SkM_CBs.txt ${WHITELIST}
echo "Done extracting CBs for D1_D0SkM!"

samtools view -1 -b \
-@ ${NCORES} \
--tag-file CB:${CBs_DIR}/yPAP-Pro_SkM-D0.txt \
${BAMPATH} \
> ${DATADIR}/Vis_yPAP_3D_STARsolo/D0SkM.Aligned.sortedByCoord.dedup.out.bam

samtools index -@ ${NCORES} ${DATADIR}/Vis_yPAP_3D_STARsolo/D0SkM.Aligned.sortedByCoord.dedup.out.bam


# JSONPATH="/workdir/dwm269/totalRNA/data/visium_images/Vis_yPAP3/V10T11-086-D1_D2SkM.json"
# python /home/dwm269/DWM_utils/sc_utils/TXG_helpers/json2CBs.py ${JSONPATH} ${DATADIR}/V10T11-086-D1_D2SkM_CBs.txt ${WHITELIST}
echo "Done extracting CBs for C1_D7SkM!"

samtools view -1 -b \
-@ ${NCORES} \
--tag-file CB:${CBs_DIR}/yPAP-Pro_SkM-D5.txt \
${BAMPATH} \
> ${DATADIR}/Vis_yPAP_3D_STARsolo/D2SkM.Aligned.sortedByCoord.dedup.out.bam

samtools index -@ ${NCORES} ${DATADIR}/Vis_yPAP_3D_STARsolo/D2SkM.Aligned.sortedByCoord.dedup.out.bam
