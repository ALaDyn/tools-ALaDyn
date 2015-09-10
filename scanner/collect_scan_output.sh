#! /bin/bash

module load compilers/gcc-4.9.0


#per il seguente, copiarsi dal prepare_scan_input la riga che genera tutti i valori (in questo caso di bulk lengths scannerizzate)
columns_values=$(awk 'BEGIN{for(i=2.0;i<=6.0;i+=2.0)print i}')
SIM_HEADER="ion_"

#columns_values=$(awk 'BEGIN{for(i=2.0;i<=10.0;i+=2.0)print i}')
#SIM_HEADER="pre_"

#columns_values=$(awk 'BEGIN{for(i=2.0;i<=6.0;i+=2.0)print i}')
#SIM_HEADER="15deg_"

#columns_values=$(awk 'BEGIN{for(i=2.0;i<=10.0;i+=2.0)print i}')
#SIM_HEADER="ch2_"

PARTICLE_TYPE="A2-ion"
#PARTICLE_TYPE="Electron"
EXP_FIT_SOFTWARE=$HOME/bin/exponential_fit
SPEC_DECODER=$HOME/bin/leggi_diag

#il tool seguente e` scan-columns in tools-Propaga
#serve per splittare i risultati collezionati dallo script su diversi files,
# in funzione di un parametro, per produrre piu` plots
SCANNER=$HOME/bin/scan-columns
BILINEAR_FILTER=$HOME/bin/interpolate
X_INTERPOLATION=1024
Y_INTERPOLATION=512
GNUPLOT=`which gnuplot`

#the following is the value used to split the collected files. It's usual to split along the bulk length
column=4
# 1:{PREPLASMA_LENGTH} 2:{DENSITY} 3:{RAMP_LENGTH} 4:{BULK_LENGTH} 5:{CONT_LENGTH} 
# 6:{MAX_ENERGY} 7:{TOT_ENERGY} 8:{TOT_ENERGY_INTEGRATED_SPECTRUM} 9:{FRONT_TOT_ENERGY_INTEGRATED_SPECTRUM}
#10:{REAR_TOT_ENERGY_INTEGRATED_SPECTRUM} 11:{AVE_ENERGY} 12:{TOT_NUMBER} 13:{FRONT_AVE_ENERGY}
#14:{FRONT_TOT_NUMBER} 15:{REAR_AVE_ENERGY} 16:{REAR_TOT_NUMBER} 17:{DIAG_MAX_ENERGY} 18:{DIAG_TOT_ENERGY}"

DIAG_STEP_TO_BE_READ=10
SPEC_TIME_TO_BE_READ=100
SPEC_VERSION=4
DIAG_VERSION=3
OUTPUT_FILE="energy_scan_${SIM_HEADER}${PARTICLE_TYPE}.txt"

SIMULATION_FOLDERS=($(find . -mindepth 1 -maxdepth 1 -type d -name "${SIM_HEADER}*"))


# DESCRIPTION OF DIAG COLUMNS
#              1:timestep 2:Efields2 
# (Electrons)  3:Etot     4:Emax   5:Jz   6:px   7:py   8:pz   9:sigma_px  10:sigma_py  11:sigma_pz  12:mean_charge  13:charge_per_cell
# (A1-ions)   14:Etot    15:Emax  16:Jz  17:px  18:py  19:pz  20:sigma_px  21:sigma_py  22:sigma_pz  23:mean_charge  24:charge_per_cell
# (A2-ions)   25:Etot    26:Emax  27:Jz  28:px  29:py  30:pz  31:sigma_px  32:sigma_py  33:sigma_pz  34:mean_charge  35:charge_per_cell
if [ ${PARTICLE_TYPE} == "Electron" ]
then
 COLUMN_MAX_ENERGY=4
 COLUMN_TOT_ENERGY=3
elif [ ${PARTICLE_TYPE} == "A2-ion" ]
then
 COLUMN_MAX_ENERGY=26
 COLUMN_TOT_ENERGY=25
else
 echo "Unrecognized particle type"
 exit
fi

rm -f ${OUTPUT_FILE}
touch ${OUTPUT_FILE}

printf "#PreplasmaLength;Density;RampLength;BulkLength;ContLength;${PARTICLE_TYPE}MaxEnergy;${PARTICLE_TYPE}TotEnergy;${PARTICLE_TYPE}TotEnergy_IS;${PARTICLE_TYPE}FrontTotEnergy;${PARTICLE_TYPE}RearTotEnergy;${PARTICLE_TYPE}AveEnergy;${PARTICLE_TYPE}TotNumber;Front${PARTICLE_TYPE}AveEnergy;Front${PARTICLE_TYPE}TotNumber;Rear${PARTICLE_TYPE}AveEnergy;Rear${PARTICLE_TYPE}TotNumber;DiagMaxEnergy;DiagTotEnergy\n" >> ${OUTPUT_FILE}

for sim in "${SIMULATION_FOLDERS[@]}"
do
 PREPLASMA_LENGTH=($(echo $sim | awk -F'_' '{print $3}'))
 RAMP_LENGTH=($(echo $sim | awk -F'_' '{print $7}'))
 DENSITY=($(echo $sim | awk -F'_' '{print $5}'))
 BULK_LENGTH=($(echo $sim | awk -F'_' '{print $9}'))
 CONT_LENGTH=($(echo $sim | awk -F'_' '{print $11}'))
 cd $sim

 if [ -f "diag${DIAG_STEP_TO_BE_READ}.dat" ];
 then
  ${SPEC_DECODER} diag${DIAG_STEP_TO_BE_READ}.dat v${DIAG_VERSION}
  DIAG_MAX_ENERGY=($(tail -1 diag${DIAG_STEP_TO_BE_READ}.dat.particles.txt | awk '{print $'${COLUMN_MAX_ENERGY}'}'))
  DIAG_TOT_ENERGY=($(tail -1 diag${DIAG_STEP_TO_BE_READ}.dat.particles.txt | awk '{print $'${COLUMN_TOT_ENERGY}'}'))
 else
  DIAG_MAX_ENERGY="-1"
  DIAG_TOT_ENERGY="-1"
 fi

 if [ -f "spec${DIAG_STEP_TO_BE_READ}.dat" ];
 then
  ${SPEC_DECODER} spec${DIAG_STEP_TO_BE_READ}.dat v${SPEC_VERSION}
  aveData=($( ${EXP_FIT_SOFTWARE} -scan spec${DIAG_STEP_TO_BE_READ}.dat_${PARTICLE_TYPE}_${SPEC_TIME_TO_BE_READ}.* ))
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

for value in ${columns_values}
do
 $SCANNER -in energy_scan_${SIM_HEADER}${PARTICLE_TYPE}.txt -out energy_scan_${SIM_HEADER}${PARTICLE_TYPE}_$value.txt -select $column $value
 ${BILINEAR_FILTER} -cx 1 -cy 3 -ce 6 -nx ${X_INTERPOLATION} -ny ${Y_INTERPOLATION} -file energy_scan_${SIM_HEADER}${PARTICLE_TYPE}_${value}.txt -gnuplot -title "bulk length ${value} {/Symbol.ttf m}m" -xlabel "Preplasma length ({/Symbol.ttf m}m)" -ylabel "Ramp length ({/Symbol.ttf m}m)" -cblabel "Maximum proton energy (MeV)"
 $GNUPLOT plot.plt
 ${BILINEAR_FILTER} -cx 1 -cy 3 -ce 15 -nx ${X_INTERPOLATION} -ny ${Y_INTERPOLATION} -file energy_scan_${SIM_HEADER}${PARTICLE_TYPE}_${value}.txt -gnuplot -title "bulk length ${value} {/Symbol.ttf m}m" -xlabel "Preplasma length ({/Symbol.ttf m}m)" -ylabel "Ramp length ({/Symbol.ttf m}m)" -cblabel "Average proton temperature (MeV)"
 $GNUPLOT plot.plt
 ${BILINEAR_FILTER} -cx 1 -cy 3 -ce 18 -nx ${X_INTERPOLATION} -ny ${Y_INTERPOLATION} -file energy_scan_${SIM_HEADER}${PARTICLE_TYPE}_${value}.txt -gnuplot -title "bulk length ${value} {/Symbol.ttf m}m" -xlabel "Preplasma length ({/Symbol.ttf m}m)" -ylabel "Ramp length ({/Symbol.ttf m}m)" -cblabel "Total proton energy ({/Symbol.ttf m}J)" -cb_magn 1000000
 $GNUPLOT plot.plt
done

