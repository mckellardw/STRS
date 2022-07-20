#!/usr/bin/bash

# Usage:
#   kb.sh $OUTDIR $KB_IDX $WHITELIST $CHEMISTRY $LOG $THREADS $MEMLIMIT $R1FQ $R2FQ

# Get params
OUTDIR=$1
KB_IDX=$2
WHITELIST=$3
CHEMISTRY=$4
LOG=$5
THREADS=$6
MEMLIMIT=$7
R1FQ=$8
R2FQ=$9

# Check params
#TODO

# Set up output directory
mkdir -p ${OUTDIR}
cd ${OUTDIR}
touch ${LOG}

# Add params to log file
echo "~~ Parameters ~~~" >> ${LOG}
echo "Output directory: " ${OUTDIR}  >> ${LOG}
echo "Index used: " ${KB_IDX} >> ${LOG}
echo "Whitelist used: " ${WHITELIST} >> ${LOG}
echo "Chemistry: " ${CHEMISTRY} >> ${LOG}
echo "Barcode/UMI read file: " ${R1FQ} >> ${LOG}
echo "RNA read file:" ${R2FQ} >> ${LOG}
echo " " >> ${LOG}

# Pseudoalign and generate .bus file
echo "~~~Pseudoaligning with `kallisto bus`... " >> ${LOG}
kallisto bus \
-i ${KB_IDX} \
-x ${CHEMISTRY} \
--fr-stranded \
-o ${OUTDIR} \
-t ${THREADS} \
--verbose \
${R1FQ} ${R2FQ} 2>> ${LOG}
echo " " >> ${LOG}

# No barcode correction for VASAdrop - doesn't have a whitelist!!
# Correct cell/spot barcodes
# echo "~~~Correcting barcodes... " >> ${LOG}
# bustools correct \
# --whitelist ${WHITELIST} \
# -o output.corrected.bus \
# output.bus 2>> ${LOG}
# echo " " >> ${LOG}

# Sort .bus file
echo "~~~Sorting output bus... " >> ${LOG}
bustools sort \
-t ${THREADS} \
-m ${MEMLIMIT} \
-o output.sorted.bus \
output.bus 2>> ${LOG}
echo " " >> ${LOG}

# Inspect outputs
echo "~~~Inspecting sorted BUS file..." >> ${LOG}
bustools inspect \
--whitelist ${WHITELIST} \
--ecmap matrix.ec \
output.sorted.bus 2>> ${LOG}

echo "~~~Writing to inspect.sorted.bus.json" >> ${LOG}
bustools inspect \
--whitelist ${WHITELIST} \
--ecmap matrix.ec \
--output inspect.sorted.bus.json \
output.sorted.bus
echo "Done!" >> ${LOG}

# Convert bus file to text for easier counting
# echo "~~~Converting bus to text... " >> ${LOG}
# bustools text \
# -o output.corrected.txt \
# output.corrected.bus
