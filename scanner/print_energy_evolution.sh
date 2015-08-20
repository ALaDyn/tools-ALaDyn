#! /bin/bash
module load gnu/4.9.2
module load gnuplot/5.0.0

SIM_HEADER="pre_"
LEGGI_DIAG=$HOME/bin/leggi_diag
DIAG_VERSION=4
DIAG_TIME_TO_BE_READ="10"
ANOTHER1_DIAG_TIME_TO_BE_READ="04"
ANOTHER2_DIAG_TIME_TO_BE_READ="05"
GNUPLOT=/cineca/prod/tools/gnuplot/5.0.0/gnu--4.9.2/bin/gnuplot
IMAGE_TYPE=png
SIZEX=960
SIZEY=540
SIMULATION_FOLDERS=($(find . -name "${SIM_HEADER}*" -type d))
PRINT_ENERGY_COMPARISON=true
PRINT_JUST_P=false
PRINT_JUST_E=false

FIELD_ETOT_COLUMN=2
PROTON_ETOT_COLUMN=25
PROTON_EMAX_COLUMN=26
ELECTRON_ETOT_COLUMN=3
ELECTRON_EMAX_COLUMN=4
ION_ETOT_COLUMN=14
ION_EMAX_COLUMN=15



for sim in "${SIMULATION_FOLDERS[@]}"
do
 PREPLASMA_LENGTH=($(echo $sim | awk -F'_' '{print $2}'))
 RAMP_LENGTH=($(echo $sim | awk -F'_' '{print $6}'))
 DENSITY=($(echo $sim | awk -F'_' '{print $4}'))
 BULK_LENGTH=($(echo $sim | awk -F'_' '{print $8}'))
 CONT_LENGTH=($(echo $sim | awk -F'_' '{print $10}'))
 
 cd $sim
 filename=`ls -1 diag${DIAG_TIME_TO_BE_READ}.dat`
 another1_filename=`ls -1 diag${ANOTHER1_DIAG_TIME_TO_BE_READ}.dat`
 another2_filename=`ls -1 diag${ANOTHER2_DIAG_TIME_TO_BE_READ}.dat`

 files=$(ls ${filename} ${another1_filename} ${another2_filename} 2> /dev/null | wc -w)
 if [ $files ]
# if [ -f $filename ]
 then
  ${LEGGI_DIAG} $filename v${DIAG_VERSION}
  ${LEGGI_DIAG} $another1_filename v${DIAG_VERSION}
  ${LEGGI_DIAG} $another2_filename v${DIAG_VERSION}

  cat ${another1_filename}.particles.txt ${another2_filename}.particles.txt ${filename}.particles.txt > temp.txt
  sort -k1 -n temp.txt > temp_sort.txt
  filename=temp_sort

  GNUPLOT_FILE=diag.plt
  rm -f ${GNUPLOT_FILE}
  touch ${GNUPLOT_FILE}
  chmod 755 ${GNUPLOT_FILE}

  if ${PRINT_JUST_P}
  then
   printf "#!/gnuplot\n" > ${GNUPLOT_FILE}
   printf "FILE_IN='%s.txt'\n" "${filename}" >> ${GNUPLOT_FILE}
   printf "FILE_OUT='diag_prot_pre_%s_ramp_%s_den_%s_bulk_%s_cont_%s.png'\n" "${PREPLASMA_LENGTH}" "${RAMP_LENGTH}" "${DENSITY}" "${BULK_LENGTH}" "${CONT_LENGTH}" >> ${GNUPLOT_FILE}
   printf "c=10/3\n" >> ${GNUPLOT_FILE}
   printf "set terminal ${IMAGE_TYPE} truecolor enhanced size ${SIZEX},${SIZEY}\n" >> ${GNUPLOT_FILE}
   printf "set output FILE_OUT\n" >> ${GNUPLOT_FILE}
   printf "set title 'prot\_pre\_%s\_ramp\_%s\_den\_%s\_bulk\_%s\_cont\_%s.png'\n" "${PREPLASMA_LENGTH}" "${RAMP_LENGTH}" "${DENSITY}" "${BULK_LENGTH}" "${CONT_LENGTH}" >> ${GNUPLOT_FILE}
  # printf "set xrange[22:29]\n" >> ${GNUPLOT_FILE}
  # printf "set yrange[-nc*1.1:nc*1.1]\n" >> ${GNUPLOT_FILE}
   printf "set xlabel 't (fs)'\n" >> ${GNUPLOT_FILE}
   printf "set ylabel 'E_{max} (MeV)'\n" >> ${GNUPLOT_FILE}
   printf "set y2label 'E_{tot} (uJ)'\n" >> ${GNUPLOT_FILE}
   printf "set key top left\n" >> ${GNUPLOT_FILE}
  # printf "set format y2 '%.6f'\n" >> ${GNUPLOT_FILE}
  # printf "set format y2 '%.0s 10^{%T}'\n" >> ${GNUPLOT_FILE}
   printf "set ytics nomirror\n" >> ${GNUPLOT_FILE}
   printf "set y2tics\n" >> ${GNUPLOT_FILE}
   printf "plot FILE_IN u (\$1*c):%s w points pt 7 ps 2 lc rgb 'red' t 'E_{max}' axes x1y1,\\" "${PROTON_EMAX_COLUMN}" >> ${GNUPLOT_FILE}
  # printf "plot FILE_IN u (\$1*c):%s w lines lt 1 lw 5 lc rgb 'red' t 'E_{max}' axes x1y1,\\" "${PROTON_EMAX_COLUMN}" >> ${GNUPLOT_FILE}
   printf "\n" >> ${GNUPLOT_FILE}
   printf "FILE_IN u (\$1*c):(\$%s*1E6) w points pt 7 ps 2 lc rgb 'blue' t 'E_{tot}' axes x1y2\n" "${PROTON_ETOT_COLUMN}" >> ${GNUPLOT_FILE}
  # printf "FILE_IN u (\$1*c):(\$%s*1E6) w lines lt 1 lw 5 lc rgb 'blue' t 'E_{tot}' axes x1y2\n" "${PROTON_ETOT_COLUMN}" >> ${GNUPLOT_FILE}
  
   $GNUPLOT ${GNUPLOT_FILE}
   mv diag_prot_pre_${PREPLASMA_LENGTH}_ramp_${RAMP_LENGTH}_den_${DENSITY}_bulk_${BULK_LENGTH}_cont_${CONT_LENGTH}.png ../
  fi

  if ${PRINT_JUST_E}
  then
   printf "#!/gnuplot\n" > ${GNUPLOT_FILE}
   printf "FILE_IN='%s.txt'\n" "${filename}" >> ${GNUPLOT_FILE}
   printf "FILE_OUT='diag_el_pre_%s_ramp_%s_den_%s_bulk_%s_cont_%s.png'\n" "${PREPLASMA_LENGTH}" "${RAMP_LENGTH}" "${DENSITY}" "${BULK_LENGTH}" "${CONT_LENGTH}" >> ${GNUPLOT_FILE}
   printf "c=10/3\n" >> ${GNUPLOT_FILE}
   printf "set terminal ${IMAGE_TYPE} truecolor enhanced size ${SIZEX},${SIZEY}\n" >> ${GNUPLOT_FILE}
   printf "set output FILE_OUT\n" >> ${GNUPLOT_FILE}
   printf "set title 'el\_pre\_%s\_ramp\_%s\_den\_%s\_bulk\_%s\_cont\_%s.png'\n" "${PREPLASMA_LENGTH}" "${RAMP_LENGTH}" "${DENSITY}" "${BULK_LENGTH}" "${CONT_LENGTH}" >> ${GNUPLOT_FILE}
  # printf "set xrange[22:29]\n" >> ${GNUPLOT_FILE}
  # printf "set yrange[-nc*1.1:nc*1.1]\n" >> ${GNUPLOT_FILE}
   printf "set xlabel 't (fs)'\n" >> ${GNUPLOT_FILE}
   printf "set ylabel 'E_{max} (MeV)'\n" >> ${GNUPLOT_FILE}
   printf "set y2label 'E_{tot} (uJ)'\n" >> ${GNUPLOT_FILE}
   printf "set key top left\n" >> ${GNUPLOT_FILE}
  # printf "set format y2 '%.4f'\n" >> ${GNUPLOT_FILE}
  # printf "set format y2 '%.0s 10^{%T}'\n" >> ${GNUPLOT_FILE}
   printf "set ytics nomirror\n" >> ${GNUPLOT_FILE}
   printf "set y2tics\n" >> ${GNUPLOT_FILE}
   printf "plot FILE_IN u (\$1*c):%s w points pt 7 ps 2 lc rgb 'red' t 'E_{max}' axes x1y1,\\" "${ELECTRON_EMAX_COLUMN}" >> ${GNUPLOT_FILE}
  # printf "plot FILE_IN u (\$1*c):%s w lines lt 1 lw 5 lc rgb 'red' t 'E_{max}' axes x1y1,\\" "${ELECTRON_EMAX_COLUMN}" >> ${GNUPLOT_FILE}
   printf "\n" >> ${GNUPLOT_FILE}
   printf "FILE_IN u (\$1*c):(\$%s*1E6) w points pt 7 ps 2 lc rgb 'blue' t 'E_{tot}' axes x1y2\n" "${ELECTRON_ETOT_COLUMN}" >> ${GNUPLOT_FILE}
  # printf "FILE_IN u (\$1*c):(\$%s*1E6) w lines lt 1 lw 5 lc rgb 'blue' t 'E_{tot}' axes x1y2\n" "${ELECTRON_ETOT_COLUMN}" >> ${GNUPLOT_FILE}

   $GNUPLOT ${GNUPLOT_FILE}
   mv diag_el_pre_${PREPLASMA_LENGTH}_ramp_${RAMP_LENGTH}_den_${DENSITY}_bulk_${BULK_LENGTH}_cont_${CONT_LENGTH}.png ../
  fi

  if ${PRINT_ENERGY_COMPARISON}
  then
   printf "#!/gnuplot\n" > ${GNUPLOT_FILE}
   printf "FILE_IN='%s.txt'\n" "${filename}" >> ${GNUPLOT_FILE}
   printf "FILE_OUT='diag_rel_pre_%s_ramp_%s_den_%s_bulk_%s_cont_%s.png'\n" "${PREPLASMA_LENGTH}" "${RAMP_LENGTH}" "${DENSITY}" "${BULK_LENGTH}" "${CONT_LENGTH}" >> ${GNUPLOT_FILE}
   printf "c=10/3\n" >> ${GNUPLOT_FILE}
   printf "set terminal ${IMAGE_TYPE} truecolor enhanced size ${SIZEX},${SIZEY}\n" >> ${GNUPLOT_FILE}
   printf "set output FILE_OUT\n" >> ${GNUPLOT_FILE}
   printf "set title 'pre\_%s\_ramp\_%s\_den\_%s\_bulk\_%s\_cont\_%s.png'\n" "${PREPLASMA_LENGTH}" "${RAMP_LENGTH}" "${DENSITY}" "${BULK_LENGTH}" "${CONT_LENGTH}" >> ${GNUPLOT_FILE}
  # printf "set xrange[22:29]\n" >> ${GNUPLOT_FILE}
  # printf "set yrange[-nc*1.1:nc*1.1]\n" >> ${GNUPLOT_FILE}
   printf "set xlabel 't (fs)'\n" >> ${GNUPLOT_FILE}
   printf "set ylabel 'E_{max} (MeV)'\n" >> ${GNUPLOT_FILE}
   printf "set y2label 'E/E_{tot}'\n" >> ${GNUPLOT_FILE}
   printf "set y2range[0:100]\n" >> ${GNUPLOT_FILE}
   printf "set format y2 '%%g %%%%'\n" >> ${GNUPLOT_FILE}
   printf "set key top left\n" >> ${GNUPLOT_FILE}
  # printf "set format y2 '%.4f'\n" >> ${GNUPLOT_FILE}
  # printf "set format y2 '%.0s 10^{%T}'\n" >> ${GNUPLOT_FILE}
   printf "set ytics nomirror\n" >> ${GNUPLOT_FILE}
   printf "set y2tics\n" >> ${GNUPLOT_FILE}
   printf "plot FILE_IN u (\$1*c):%s w lines lt 1 lw 5 lc rgb 'red' t 'Emax_{el}' axes x1y1,\\"  "${ELECTRON_EMAX_COLUMN}" >> ${GNUPLOT_FILE}
   printf "\n" >> ${GNUPLOT_FILE}
   printf "FILE_IN u (\$1*c):(\$%s*1E6)*100.0/((\$%s*1E6)+(\$%s*1E6)+(\$%s*1E6)+(\$%s*1E6)) w lines lt 1 lw 5 lc rgb 'blue' t 'Etot_{el}/E_{tot}' axes x1y2,\\" "${ELECTRON_ETOT_COLUMN}" "${FIELD_ETOT_COLUMN}" "${ELECTRON_ETOT_COLUMN}" "${PROTON_ETOT_COLUMN}" "${ION_ETOT_COLUMN}" >> ${GNUPLOT_FILE}
   printf "\n" >> ${GNUPLOT_FILE}
   printf "FILE_IN u (\$1*c):%s w lines lt 1 lw 5 lc rgb 'orange' t 'Emax_{pr}' axes x1y1,\\"  "${PROTON_EMAX_COLUMN}" >> ${GNUPLOT_FILE}
   printf "\n" >> ${GNUPLOT_FILE}
   printf "FILE_IN u (\$1*c):(\$%s*1E6)*100.0/((\$%s*1E6)+(\$%s*1E6)+(\$%s*1E6)+(\$%s*1E6)) w lines lt 1 lw 5 lc rgb 'cyan' t 'Etot_{pr}/E_{tot}' axes x1y2" "${PROTON_ETOT_COLUMN}" "${FIELD_ETOT_COLUMN}" "${ELECTRON_ETOT_COLUMN}" "${PROTON_ETOT_COLUMN}" "${ION_ETOT_COLUMN}" >> ${GNUPLOT_FILE}
   printf "\n" >> ${GNUPLOT_FILE}

   $GNUPLOT ${GNUPLOT_FILE}

   printf "#!/gnuplot\n" > ${GNUPLOT_FILE}
   printf "FILE_IN='%s.txt'\n" "${filename}" >> ${GNUPLOT_FILE}
   printf "FILE_OUT='diag_abs_pre_%s_ramp_%s_den_%s_bulk_%s_cont_%s.png'\n" "${PREPLASMA_LENGTH}" "${RAMP_LENGTH}" "${DENSITY}" "${BULK_LENGTH}" "${CONT_LENGTH}" >> ${GNUPLOT_FILE}
   printf "c=10/3\n" >> ${GNUPLOT_FILE}
   printf "set terminal ${IMAGE_TYPE} truecolor enhanced size ${SIZEX},${SIZEY}\n" >> ${GNUPLOT_FILE}
   printf "set output FILE_OUT\n" >> ${GNUPLOT_FILE}
   printf "set title 'pre\_%s\_ramp\_%s\_den\_%s\_bulk\_%s\_cont\_%s.png'\n" "${PREPLASMA_LENGTH}" "${RAMP_LENGTH}" "${DENSITY}" "${BULK_LENGTH}" "${CONT_LENGTH}" >> ${GNUPLOT_FILE}
  # printf "set xrange[22:29]\n" >> ${GNUPLOT_FILE}
  # printf "set yrange[-nc*1.1:nc*1.1]\n" >> ${GNUPLOT_FILE}
   printf "set xlabel 't (fs)'\n" >> ${GNUPLOT_FILE}
   printf "set ylabel 'E_{max} (MeV)'\n" >> ${GNUPLOT_FILE}
   printf "set y2label 'E_{tot} (uJ)'\n" >> ${GNUPLOT_FILE}
   printf "set key top left\n" >> ${GNUPLOT_FILE}
  # printf "set format y2 '%.4f'\n" >> ${GNUPLOT_FILE}
  # printf "set format y2 '%.0s 10^{%T}'\n" >> ${GNUPLOT_FILE}
   printf "set ytics nomirror\n" >> ${GNUPLOT_FILE}
   printf "set y2tics\n" >> ${GNUPLOT_FILE}
   printf "plot FILE_IN u (\$1*c):((\$%s*1E6)+(\$%s*1E6)+(\$%s*1E6)) w lines lt 1 lw 5 lc rgb 'black' t 'E_{tot}' axes x1y2,\\" "${ELECTRON_ETOT_COLUMN}" "${PROTON_ETOT_COLUMN}" "${ION_ETOT_COLUMN}" >> ${GNUPLOT_FILE}
   printf "\n" >> ${GNUPLOT_FILE}
   printf "FILE_IN u (\$1*c):%s w lines lt 1 lw 5 lc rgb 'red' t 'el E_{max}' axes x1y1,\\" "${ELECTRON_EMAX_COLUMN}" >> ${GNUPLOT_FILE}
   printf "\n" >> ${GNUPLOT_FILE}
   printf "FILE_IN u (\$1*c):(\$%s*1E6) w lines lt 1 lw 5 lc rgb 'blue' t 'el E_{tot}' axes x1y2,\\" "${ELECTRON_ETOT_COLUMN}" >> ${GNUPLOT_FILE}
   printf "\n" >> ${GNUPLOT_FILE}
   printf "FILE_IN u (\$1*c):%s w lines lt 1 lw 5 lc rgb 'orange' t 'pr E_{max}' axes x1y1,\\" "${PROTON_EMAX_COLUMN}" >> ${GNUPLOT_FILE}
   printf "\n" >> ${GNUPLOT_FILE}
   printf "FILE_IN u (\$1*c):(\$%s*1E6) w lines lt 1 lw 5 lc rgb 'cyan' t 'pr E_{tot}' axes x1y2" "${PROTON_ETOT_COLUMN}" >> ${GNUPLOT_FILE}
   printf "\n" >> ${GNUPLOT_FILE}

   $GNUPLOT ${GNUPLOT_FILE}
   mv *.png ../
  fi

 fi

 cd ..

done

