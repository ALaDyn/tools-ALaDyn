#! /bin/bash

SIM_HEADER="pre_"

PARTICLE_TYPE='Proton'
FIT_SOFTWARE=/home/HPC/sinigardi/bin/fit
SPEC_DECODER=/home/HPC/sinigardi/bin/leggi_diag

FINAL_LINE_NUMBER=146
DIAG_STEP_TO_BE_READ=10
SPEC_TIME_TO_BE_READ=100

OUTPUT_FILE="energy_scan.txt"




SIMULATION_FOLDERS=($(find . -name "${SIM_HEADER}*" -type d))

rm -f ${OUTPUT_FILE}
touch ${OUTPUT_FILE}


printf "# PreplasmaLength \t Density \t RampLength \t BulkLength \t ${PARTICLE_TYPE}MaxEnergy \t ${PARTICLE_TYPE}TotEnergy \t ${PARTICLE_TYPE}AveEnergy \t ${PARTICLE_TYPE}TotNumber \t Sel${PARTICLE_TYPE}AveEnergy \t Sel${PARTICLE_TYPE}TotNumber \n" >> ${OUTPUT_FILE} 

for sim in "${SIMULATION_FOLDERS[@]}"
do
 PREPLASMA_LENGTH=($(echo $sim | awk -F'_' '{print $2}'))
 RAMP_LENGTH=($(echo $sim | awk -F'_' '{print $6}'))
 DENSITY=($(echo $sim | awk -F'_' '{print $4}'))
# BULK_LENGTH=($(echo $sim | awk -F'_' '{print $8}'))
 BULK_LENGTH="0.6"

cd $sim
if [ -f "diag${DIAG_STEP_TO_BE_READ}.dat" ];
then
 PROTON_MAX_ENERGY=($(head -${FINAL_LINE_NUMBER} diag${DIAG_STEP_TO_BE_READ}.dat |tail -1 |awk '{print $4}'))
 PROTON_TOT_ENERGY=($(head -${FINAL_LINE_NUMBER} diag${DIAG_STEP_TO_BE_READ}.dat |tail -1 |awk '{print $3}'))
else
 PROTON_MAX_ENERGY="0.0"
 PROTON_TOT_ENERGY="0.0"
fi


if [ -f "spec${DIAG_STEP_TO_BE_READ}.dat" ];
then
 ${SPEC_DECODER} spec${DIAG_STEP_TO_BE_READ}.dat
 aveData=($( ${FIT_SOFTWARE} -scan spec${DIAG_STEP_TO_BE_READ}.dat_${PARTICLE_TYPE}_${SPEC_TIME_TO_BE_READ}.txt  ))
else
 aveData[0]="0.0"
 aveData[1]="0.0"
 aveData[2]="0.0"
 aveData[3]="0.0"
fi

PROTON_AVE_ENERGY=${aveData[0]}
PROTON_TOT_NUMBER=${aveData[1]}
SEL_PROTON_AVE_ENERGY=${aveData[2]}
SEL_PROTON_TOT_NUMBER=${aveData[3]}

cd ..

#read -p "Press [Enter] key to continue..."
printf '\t      %.1f \t%.1f\t     %.1f   \t%.1f\t%.2f\t%.3e \t%s\t%s\t%s\t%s' "${PREPLASMA_LENGTH}" "${DENSITY}" "${RAMP_LENGTH}" "${BULK_LENGTH}" "${PROTON_MAX_ENERGY}" "${PROTON_TOT_ENERGY}" "${PROTON_AVE_ENERGY}" "${PROTON_TOT_NUMBER}" "${SEL_PROTON_AVE_ENERGY}" "${SEL_PROTON_TOT_NUMBER}" >> ${OUTPUT_FILE}

printf '\n'  >> ${OUTPUT_FILE}


done

./split_energy_scan.sh

