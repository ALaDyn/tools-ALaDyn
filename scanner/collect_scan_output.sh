#! /bin/bash
module load gnu/4.9.2
SIM_HEADER="pre_"

PARTICLE_TYPE="Hion"
#PARTICLE_TYPE="Electron"
EXP_FIT_SOFTWARE=$HOME/bin/exponential_fit
SPEC_DECODER=$HOME/bin/leggi_diag

#il tool seguente e` scan-columns in tools-Propaga
#serve per splittare i risultati collezionati dallo script su diversi files,
# in funzione di un parametro, per produrre piu` plots
SCANNER=$HOME/bin/scan-columns
DO_SCAN=false

#per il seguente, copiarsi dal prepare_scan_input la riga che genera tutti i valori (in questo caso di bulk lengths scannerizzate)
columns_values=$(awk 'BEGIN{for(i=2.0;i<=10.0;i+=1.0)print i}')
column=4
# 1:{PREPLASMA_LENGTH} 2:{DENSITY} 3:{RAMP_LENGTH} 4:{BULK_LENGTH} 5:{CONT_LENGTH} 
# 6:{PROTON_MAX_ENERGY} 7:{PROTON_TOT_ENERGY} 8:{PROTON_AVE_ENERGY} 9:{PROTON_TOT_NUMBER}" 
# 10:{FRONT_PROTON_AVE_ENERGY} 11:{FRONT_PROTON_TOT_NUMBER} 12:{REAR_PROTON_AVE_ENERGY} 13:{REAR_PROTON_TOT_NUMBER}

DIAG_STEP_TO_BE_READ=10
SPEC_TIME_TO_BE_READ=100
SPEC_VERSION=4
DIAG_VERSION=3
OUTPUT_FILE="energy_scan_${PARTICLE_TYPE}.txt"
NPHYS_OVER_NMACRO=696800

SIMULATION_FOLDERS=($(find . -name "${SIM_HEADER}*" -type d))


if [ ${PARTICLE_TYPE} == "Electron" ]
then
 COLUMN_MAX_ENERGY=3
 COLUMN_TOT_ENERGY=2
elif [ ${PARTICLE_TYPE} == "Hion" ]
then
 COLUMN_MAX_ENERGY=7
 COLUMN_TOT_ENERGY=6
else
 echo "Unrecognized particle type"
 exit
fi

rm -f ${OUTPUT_FILE}
touch ${OUTPUT_FILE}

printf "#PreplasmaLength;Density;RampLength;BulkLength;ContLength;${PARTICLE_TYPE}MaxEnergy;${PARTICLE_TYPE}TotEnergy;${PARTICLE_TYPE}TotEnergy_IS;${PARTICLE_TYPE}FrontTotEnergy;${PARTICLE_TYPE}RearTotEnergy;${PARTICLE_TYPE}AveEnergy;${PARTICLE_TYPE}TotNumber;Front${PARTICLE_TYPE}AveEnergy;Front${PARTICLE_TYPE}TotNumber;Rear${PARTICLE_TYPE}AveEnergy;Rear${PARTICLE_TYPE}TotNumber\n" >> ${OUTPUT_FILE}

for sim in "${SIMULATION_FOLDERS[@]}"
do
 PREPLASMA_LENGTH=($(echo $sim | awk -F'_' '{print $2}'))
 RAMP_LENGTH=($(echo $sim | awk -F'_' '{print $6}'))
 DENSITY=($(echo $sim | awk -F'_' '{print $4}'))
 BULK_LENGTH=($(echo $sim | awk -F'_' '{print $8}'))
 CONT_LENGTH=($(echo $sim | awk -F'_' '{print $10}'))
 cd $sim

 if [ -f "diag${DIAG_STEP_TO_BE_READ}.dat" ];
 then
  ${SPEC_DECODER} diag${DIAG_STEP_TO_BE_READ}.dat v${DIAG_VERSION}
  DIAG_MAX_ENERGY=($(tail -1 diag${DIAG_STEP_TO_BE_READ}.dat.txt | awk '{print $'${COLUMN_MAX_ENERGY}'}'))
  DIAG_TOT_ENERGY=($(tail -1 diag${DIAG_STEP_TO_BE_READ}.dat.txt | awk '{print $'${COLUMN_TOT_ENERGY}'}'))
 else
  DIAG_MAX_ENERGY="-1"
  DIAG_TOT_ENERGY="-1"
 fi

 if [ -f "spec${DIAG_STEP_TO_BE_READ}.dat" ];
 then
  ${SPEC_DECODER} spec${DIAG_STEP_TO_BE_READ}.dat v${SPEC_VERSION}
  aveData=($( ${EXP_FIT_SOFTWARE} -nm ${NPHYS_OVER_NMACRO} -scan spec${DIAG_STEP_TO_BE_READ}.dat_${PARTICLE_TYPE}_${SPEC_TIME_TO_BE_READ}.* ))
 else
  aveData[0]="-1"
  aveData[1]="-1"
  aveData[2]="-1"
  aveData[3]="-1"
  aveData[4]="-1"
  aveData[5]="-1"
  aveData[6]="-1"
  aveData[7]="-1"
  aveData[8]="-1"
  aveData[9]="-1"
 fi

 TOT_ENERGY=${aveData[0]}
 TOT_NUMBER=${aveData[1]}
 FRONT_AVE_ENERGY=${aveData[2]}
 FRONT_TOT_NUMBER=${aveData[3]}
 REAR_AVE_ENERGY=${aveData[4]}
 REAR_TOT_NUMBER=${aveData[5]}
 MAX_ENERGY=${aveData[6]}
 TOT_ENERGY_INTEGRATED_SPECTRUM=${aveData[7]}
 FRONT_TOT_ENERGY_INTEGRATED_SPECTRUM=${aveData[8]}
 REAR_TOT_ENERGY_INTEGRATED_SPECTRUM=${aveData[9]}

 cd ..
 printf '%.1f;%.1f;%.1f;%.1f;%.2f;%g;%g;%g;%g;%g;%g;%g;%g;%g;%g;%g;%g;%g' "${PREPLASMA_LENGTH}" "${DENSITY}" "${RAMP_LENGTH}" "${BULK_LENGTH}" "${CONT_LENGTH}" "${MAX_ENERGY}" "${TOT_ENERGY}" "${TOT_ENERGY_INTEGRATED_SPECTRUM}" "${FRONT_TOT_ENERGY_INTEGRATED_SPECTRUM}" "${REAR_TOT_ENERGY_INTEGRATED_SPECTRUM}" "${AVE_ENERGY}" "${TOT_NUMBER}" "${FRONT_AVE_ENERGY}" "${FRONT_TOT_NUMBER}" "${REAR_AVE_ENERGY}" "${REAR_TOT_NUMBER}" "${DIAG_MAX_ENERGY}" "${DIAG_TOT_ENERGY}" >> ${OUTPUT_FILE}
 printf '\n'  >> ${OUTPUT_FILE}

done

if [ ${DO_SCAN} ] 
then
 for value in ${columns_values}
 do
  $SCANNER -in energy_scan_${PARTICLE_TYPE}.txt -out energy_scan_${PARTICLE_TYPE}_$value.txt -select $column $value
 done
else 
 echo "output on ${OUTPUT_FILE}"
fi


