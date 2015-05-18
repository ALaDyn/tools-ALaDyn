#! /bin/bash
module load gnu/4.9.2
module load gnuplot/5.0.0

SIM_HEADER="pre_"
LEGGI_DIAG=$HOME/bin/leggi_diag
DIAG_VERSION=4
DIAG_TIME_TO_BE_READ=10
GNUPLOT=/cineca/prod/tools/gnuplot/5.0.0/gnu--4.9.2/bin/gnuplot
IMAGE_TYPE=png
SIZEX=960
SIZEY=540
SIMULATION_FOLDERS=($(find . -name "${SIM_HEADER}*" -type d))


for sim in "${SIMULATION_FOLDERS[@]}"
do
 PREPLASMA_LENGTH=($(echo $sim | awk -F'_' '{print $2}'))
 RAMP_LENGTH=($(echo $sim | awk -F'_' '{print $6}'))
 DENSITY=($(echo $sim | awk -F'_' '{print $4}'))
 BULK_LENGTH=($(echo $sim | awk -F'_' '{print $8}'))
 CONT_LENGTH=($(echo $sim | awk -F'_' '{print $10}'))
 
 cd $sim
 filename=`ls -1 diag${DIAG_TIME_TO_BE_READ}.dat`
 if [ -f $filename ]
 then
  ${LEGGI_DIAG} $filename v${DIAG_VERSION}

  GNUPLOT_FILE=diag.plt
  rm -f ${GNUPLOT_FILE}
  touch ${GNUPLOT_FILE}
  chmod 775 ${GNUPLOT_FILE}

  printf "#!/gnuplot\n" >> ${GNUPLOT_FILE}
  printf "FILE_IN='%s.txt'\n" "${filename}" >> ${GNUPLOT_FILE}
  printf "FILE_OUT='diag_pre_%s_ramp_%s_den_%s_bulk_%s_cont_%s.png'\n" "${PREPLASMA_LENGTH}" "${RAMP_LENGTH}" "${DENSITY}" "${BULK_LENGTH}" "${CONT_LENGTH}" >> ${GNUPLOT_FILE}
  printf "c=10/3\n" >> ${GNUPLOT_FILE}
  printf "set terminal ${IMAGE_TYPE} truecolor enhanced size ${SIZEX},${SIZEY}\n" >> ${GNUPLOT_FILE}
  printf "set output FILE_OUT\n" >> ${GNUPLOT_FILE}
  printf "set title '%s'\n" "${sim}" >> ${GNUPLOT_FILE}
 # printf "set xrange[22:29]\n" >> ${GNUPLOT_FILE}
 # printf "set yrange[-nc*1.1:nc*1.1]\n" >> ${GNUPLOT_FILE}
  printf "set xlabel 't (fs)'\n" >> ${GNUPLOT_FILE}
  printf "set ylabel 'E_{max} (MeV)'\n" >> ${GNUPLOT_FILE}
  printf "set y2label 'E_{tot} (mJ)'\n" >> ${GNUPLOT_FILE}
  printf "set key top left\n" >> ${GNUPLOT_FILE}
  printf "set format y2 '%.4f'\n" >> ${GNUPLOT_FILE}
 # printf "set format y2 '%.0s 10^{%T}'\n" >> ${GNUPLOT_FILE}
  printf "set ytics nomirror\n" >> ${GNUPLOT_FILE}
  printf "set y2tics\n" >> ${GNUPLOT_FILE}
  printf "plot FILE_IN u (\$1*c):5 w lines lt 1 lw 5 lc rgb 'red' t 'E_{max}' axes x1y1,\\" >> ${GNUPLOT_FILE}
  printf "\n" >> ${GNUPLOT_FILE}
  printf "FILE_IN u (\$1*c):(\$4*1E3) w lines lt 1 lw 5 lc rgb 'blue' t 'E_{tot}' axes x1y2\n" >> ${GNUPLOT_FILE}
  
  $GNUPLOT ${GNUPLOT_FILE}
  mv diag_pre_${PREPLASMA_LENGTH}_ramp_${RAMP_LENGTH}_den_${DENSITY}_bulk_${BULK_LENGTH}_cont_${CONT_LENGTH}.png ../

 fi

 cd ..

done

