#! /bin/bash

EXP_FIT_SOFTWARE=$HOME/bin/exponential_fit.exe
LOG_FIT_SOFTWARE=$HOME/bin/logaritmic_fit.exe
SPEC_DECODER=$HOME/bin/leggi_diagspec.exe
SCANNER=$HOME/bin/scan-columns.exe
BILINEAR_FILTER=$HOME/bin/interpolate_scan_results.exe
GNUPLOT=$(which gnuplot)

SIMULATION_FOLDERS=($(find . -mindepth 1 -maxdepth 1 -name "*_*_*_*_*_*_*_*_*" -type d))
#SIMULATION_FOLDERS=($(find . -mindepth 1 -maxdepth 1 -name "${SIM_HEADER}*" -type d ))
column_values=
#column_values=(0 15 30)

PARTICLE_TYPE="A2-ion"
#PARTICLE_TYPE="Electron"
X_INTERPOLATION=512
Y_INTERPOLATION=512

DIAG_STEP_TO_BE_READ=04
SPEC_TIME_TO_BE_READ=80
SPEC_VERSION=4
DIAG_VERSION=4
OUTPUT_FILE="energy_scan_${SIM_HEADER}${PARTICLE_TYPE}.txt"

# DESCRIPTION OF DIAG COLUMNS
#              1:timestep 2:Efields2 
# (Electrons)  3:Etot     4:Emax   5:Jz   6:px   7:py   8:pz   9:sigma_px  10:sigma_py  11:sigma_pz  12:mean_charge  13:charge_per_cell
# (A1-ions)   14:Etot    15:Emax  16:Jz  17:px  18:py  19:pz  20:sigma_px  21:sigma_py  22:sigma_pz  23:mean_charge  24:charge_per_cell
# (A2-ions)   25:Etot    26:Emax  27:Jz  28:px  29:py  30:pz  31:sigma_px  32:sigma_py  33:sigma_pz  34:mean_charge  35:charge_per_cell
COLUMN_TIMESTEP=1
if [ "${PARTICLE_TYPE}" = "Electron" ]
then
 COLUMN_MAX_ENERGY=4
 COLUMN_TOT_ENERGY=3
elif [ "${PARTICLE_TYPE}" = "A2-ion" ]
then
 COLUMN_MAX_ENERGY=26
 COLUMN_TOT_ENERGY=25
else
 echo "Unrecognized particle type"
 exit
fi

rm -f "${OUTPUT_FILE}"
touch "${OUTPUT_FILE}"

printf "#ANGLE;PREPLASMA_LENGTH;PREPLASMA_DENSITY;RAMP_LENGTH;BULK_LENGTH;DENSITY;CONT_LENGTH;LASER_A_NOTE;LASER_LENGTH;MAX_ENERGY;TOT_ENERGY;TOT_ENERGY_INTEGRATED_SPECTRUM;FRONT_TOT_ENERGY_INTEGRATED_SPECTRUM;REAR_TOT_ENERGY_INTEGRATED_SPECTRUM;AVE_ENERGY;TOT_NUMBER;FRONT_AVE_ENERGY;FRONT_TOT_NUMBER;REAR_AVE_ENERGY;REAR_TOT_NUMBER;DIAG_MAX_ENERGY;DIAG_TOT_ENERGY;DIAG_MAX_ENERGY_FIT\n" >> "${OUTPUT_FILE}"

for sim in "${SIMULATION_FOLDERS[@]}"
do
 ANGLE=$(basename "$sim" | awk -F'_' '{print $1}')
 PREPLASMA_LENGTH=$(basename "$sim" | awk -F'_' '{print $2}')
 PREPLASMA_DENSITY=$(basename "$sim" | awk -F'_' '{print $3}')
 RAMP_LENGTH=$(basename "$sim" | awk -F'_' '{print $4}')
 BULK_LENGTH=$(basename "$sim" | awk -F'_' '{print $5}')
 DENSITY=$(basename "$sim" | awk -F'_' '{print $6}')
 CONT_LENGTH=$(basename "$sim" | awk -F'_' '{print $7}')
 LASER_A_NOTE=$(basename "$sim" | awk -F'_' '{print $8}')
 LASER_LENGTH=$(basename "$sim" | awk -F'_' '{print $9}')

 cd "$sim/diagnostics" || exit
 rm -f ./*.txt
 declare -a logFitData
 if [ -f "diag${DIAG_STEP_TO_BE_READ}.dat" ];
 then
  ${SPEC_DECODER} diag${DIAG_STEP_TO_BE_READ}.dat v${DIAG_VERSION}
  DIAG_MAX_ENERGY=$(tail -1 diag${DIAG_STEP_TO_BE_READ}.dat.particles.txt | awk '{print $'${COLUMN_MAX_ENERGY}'}')
  DIAG_TOT_ENERGY=$(tail -1 diag${DIAG_STEP_TO_BE_READ}.dat.particles.txt | awk '{print $'${COLUMN_TOT_ENERGY}'}')
  logFitData=($( ${LOG_FIT_SOFTWARE} -x ${COLUMN_TIMESTEP} -y ${COLUMN_MAX_ENERGY} -scan diag${DIAG_STEP_TO_BE_READ}.dat.particles.txt ))
 else
  DIAG_MAX_ENERGY="-1"
  DIAG_TOT_ENERGY="-1"
  logFitData[0]="-1"
  logFitData[1]="-1"
 fi

 declare -a expFitData
 if [ -f "spec${DIAG_STEP_TO_BE_READ}.dat" ];
 then
  ${SPEC_DECODER} spec${DIAG_STEP_TO_BE_READ}.dat v${SPEC_VERSION}
  expFitData=($( ${EXP_FIT_SOFTWARE} -min 10.0 -max 40.0 -scan spec${DIAG_STEP_TO_BE_READ}.dat_${PARTICLE_TYPE}_${SPEC_TIME_TO_BE_READ}.* ))
 else
  expFitData[0]="-1"
  expFitData[1]="-1"
  expFitData[2]="-1"
  expFitData[3]="-1"
  expFitData[4]="-1"
  expFitData[5]="-1"
  expFitData[6]="-1"
  expFitData[7]="-1"
  expFitData[8]="-1"
  expFitData[9]="-1"
 fi

 TOT_ENERGY=${expFitData[0]}
 TOT_NUMBER=${expFitData[1]}
 FRONT_AVE_ENERGY=${expFitData[2]}
 FRONT_TOT_NUMBER=${expFitData[3]}
 REAR_AVE_ENERGY=${expFitData[4]}
 REAR_TOT_NUMBER=${expFitData[5]}
 MAX_ENERGY=${expFitData[6]}
 TOT_ENERGY_INTEGRATED_SPECTRUM=${expFitData[7]}
 FRONT_TOT_ENERGY_INTEGRATED_SPECTRUM=${expFitData[8]}
 REAR_TOT_ENERGY_INTEGRATED_SPECTRUM=${expFitData[9]}
 DIAG_MAX_ENERGY_FIT=${logFitData[1]}

 cd ../.. || exit
 #  1:"${ANGLE}"                                   2:"${PREPLASMA_LENGTH}"                       3:"${PREPLASMA_DENSITY}"
 #  4:"${RAMP_LENGTH}"                             5:"${BULK_LENGTH}"                            6:"${DENSITY}"
 #  7:"${CONT_LENGTH}"                             8:"${LASER_A_NOTE}"                           9:"${LASER_LENGTH}"
 # 10:"${MAX_ENERGY}"                             11:"${TOT_ENERGY}"                            12:"${TOT_ENERGY_INTEGRATED_SPECTRUM}"
 # 13:"${FRONT_TOT_ENERGY_INTEGRATED_SPECTRUM}"   14:"${REAR_TOT_ENERGY_INTEGRATED_SPECTRUM}"   15:"${AVE_ENERGY}"
 # 16:"${TOT_NUMBER}"                             17:"${FRONT_AVE_ENERGY}"                      18:"${FRONT_TOT_NUMBER}"
 # 19:"${REAR_AVE_ENERGY}"                        20:"${REAR_TOT_NUMBER}"                       21:"${DIAG_MAX_ENERGY}"
 # 22:"${DIAG_TOT_ENERGY}"                        23:"${DIAG_MAX_ENERGY_FIT}"
 printf '%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s' "${ANGLE}" "${PREPLASMA_LENGTH}" "${PREPLASMA_DENSITY}" "${RAMP_LENGTH}" "${BULK_LENGTH}" "${DENSITY}" "${CONT_LENGTH}" "${LASER_A_NOTE}" "${LASER_LENGTH}" "${MAX_ENERGY}" "${TOT_ENERGY}" "${TOT_ENERGY_INTEGRATED_SPECTRUM}" "${FRONT_TOT_ENERGY_INTEGRATED_SPECTRUM}" "${REAR_TOT_ENERGY_INTEGRATED_SPECTRUM}" "${AVE_ENERGY}" "${TOT_NUMBER}" "${FRONT_AVE_ENERGY}" "${FRONT_TOT_NUMBER}" "${REAR_AVE_ENERGY}" "${REAR_TOT_NUMBER}" "${DIAG_MAX_ENERGY}" "${DIAG_TOT_ENERGY}" "${DIAG_MAX_ENERGY_FIT}" >> "${OUTPUT_FILE}"
 printf '\n'  >> "${OUTPUT_FILE}"

done

#the following is the value used to split the collected files. It's usual to split along the bulk length
column=5

for value in "${column_values[@]}"
do
 $SCANNER -in "${OUTPUT_FILE}" -out "${OUTPUT_FILE}_${value}.txt" -select "${column}" "${value}"

 # x: Preplasma length; y: Ramp length; cb: Maximum proton energy
 ${BILINEAR_FILTER} -cx 2 -cy 4 -ce 10 -nx ${X_INTERPOLATION} -ny ${Y_INTERPOLATION} -file "${OUTPUT_FILE}_${value}.txt" -gnuplot -title "bulk length ${value} um" -xlabel "Preplasma length (um)" -ylabel "Ramp length (um)" -cblabel "Maximum proton energy (MeV)"
 $GNUPLOT plot.plt

 # x: Preplasma length; y: Ramp length; cb: Average proton energy
 ${BILINEAR_FILTER} -cx 2 -cy 4 -ce 15 -nx ${X_INTERPOLATION} -ny ${Y_INTERPOLATION} -file "${OUTPUT_FILE}_${value}.txt" -gnuplot -title "bulk length ${value} um" -xlabel "Preplasma length (um)" -ylabel "Ramp length (um)" -cblabel "Average proton temperature (MeV)"
 $GNUPLOT plot.plt

 # x: Preplasma length; y: Ramp length; cb: Total proton energy
 ${BILINEAR_FILTER} -cx 2 -cy 4 -ce 22 -nx ${X_INTERPOLATION} -ny ${Y_INTERPOLATION} -file "${OUTPUT_FILE}_${value}.txt" -gnuplot -title "bulk length ${value} um" -xlabel "Preplasma length (um)" -ylabel "Ramp length (um)" -cblabel "Total proton energy (mJ)" -cb_magn 1000
 $GNUPLOT plot.plt
done

