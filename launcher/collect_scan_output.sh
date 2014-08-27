#! /bin/bash

SIM_HEADER="pre_"

FINAL_LINE_NUMBER=146

DIAG_STEP_TO_BE_READ=10

OUTPUT_FILE="log2.txt"




SIMULATION_FOLDERS=($(find . -name "${SIM_HEADER}*" -type d))

rm -f ${OUTPUT_FILE}
touch ${OUTPUT_FILE}

#HEADER_LINE_NUMBER=45
#head -${HEADER_LINE_NUMBER} ${SIMULATION_FOLDERS[0]}/diag${DIAG_STEP_TO_BE_READ}.dat |tail -1 >> ${OUTPUT_FILE}

printf '# PreplasmaLength   Density   RampLength    BulkLength   ProtonMaxEnergy\n' >> ${OUTPUT_FILE} 

for sim in "${SIMULATION_FOLDERS[@]}"
do
 PREPLASMA_LENGTH=($(echo $sim | awk -F'_' '{print $2}'))
 RAMP_LENGTH=($(echo $sim | awk -F'_' '{print $6}'))
 DENSITY=($(echo $sim | awk -F'_' '{print $4}'))
 BULK_LENGTH=($(echo $sim | awk -F'_' '{print $8}'))
 # printf '%s\n' "$sim" >> ${OUTPUT_FILE}
 PROTON_MAX_ENERGY=($(head -${FINAL_LINE_NUMBER}  ${sim}/diag${DIAG_STEP_TO_BE_READ}.dat |tail -1 |awk '{print $4}'))
 printf '%s %s %s %s %s\n' "${PREPLASMA_LENGTH}" "${DENSITY}" "${RAMP_LENGTH}" "${BULK_LENGTH}" "${PROTON_MAX_ENERGY}" >> ${OUTPUT_FILE}
done

