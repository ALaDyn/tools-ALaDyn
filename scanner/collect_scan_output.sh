#! /bin/bash

SIM_HEADER="pre_"

FINAL_LINE_NUMBER=146

DIAG_STEP_TO_BE_READ=10

OUTPUT_FILE="energy_scan.txt"




SIMULATION_FOLDERS=($(find . -name "${SIM_HEADER}*" -type d))

rm -f ${OUTPUT_FILE}
touch ${OUTPUT_FILE}

#HEADER_LINE_NUMBER=45
#head -${HEADER_LINE_NUMBER} ${SIMULATION_FOLDERS[0]}/diag${DIAG_STEP_TO_BE_READ}.dat |tail -1 >> ${OUTPUT_FILE}

printf '# PreplasmaLength \t Density \t RampLength \t BulkLength \t ProtonMaxEnergy \t ProtonTotEnergy\n' >> ${OUTPUT_FILE} 

for sim in "${SIMULATION_FOLDERS[@]}"
do
 PREPLASMA_LENGTH=($(echo $sim | awk -F'_' '{print $2}'))
 RAMP_LENGTH=($(echo $sim | awk -F'_' '{print $6}'))
 DENSITY=($(echo $sim | awk -F'_' '{print $4}'))
# BULK_LENGTH=($(echo $sim | awk -F'_' '{print $8}'))
 BULK_LENGTH="0.6"
if [ -f "${sim}/diag${DIAG_STEP_TO_BE_READ}.dat" ];
then
 PROTON_MAX_ENERGY=($(head -${FINAL_LINE_NUMBER}  ${sim}/diag${DIAG_STEP_TO_BE_READ}.dat |tail -1 |awk '{print $4}'))
 PROTON_TOT_ENERGY=($(head -${FINAL_LINE_NUMBER}  ${sim}/diag${DIAG_STEP_TO_BE_READ}.dat |tail -1 |awk '{print $3}'))
else
 PROTON_MAX_ENERGY="0.0"
 PROTON_TOT_ENERGY="0.0"
fi
 printf '\t      %.1f \t%.1f\t     %.1f   \t%.1f\t%.2f\t%.3e\n' "${PREPLASMA_LENGTH}" "${DENSITY}" "${RAMP_LENGTH}" "${BULK_LENGTH}" "${PROTON_MAX_ENERGY}" "${PROTON_TOT_ENERGY}" >> ${OUTPUT_FILE}
# printf '%s %s %s %s %s\n' "${PREPLASMA_LENGTH}" "${DENSITY}" "${RAMP_LENGTH}" "${BULK_LENGTH}" "${PROTON_MAX_ENERGY}" >> ${OUTPUT_FILE}
done

