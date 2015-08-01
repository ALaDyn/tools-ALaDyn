#! /bin/bash
module load gnu/4.9.2
module load gnuplot/5.0.0

SIM_HEADER="pre_"
LEGGI_DIAG=$HOME/bin/leggi_diag
DIAG_VERSION=4
declare -a SIM_OUT_TO_BE_READ
SIM_OUT_TO_BE_READ={"00","01","02","03","04","05","06","07","08","09","10"}
GNUPLOT=/cineca/prod/tools/gnuplot/5.0.0/gnu--4.9.2/bin/gnuplot
BINARY_READER=/galileo/home/userexternal/ssinigar/bin/leggi_bin
IMAGE_TYPE=png
SIZEX=960
SIZEY=540
FONTSIZE=14

XMIN=30
XMAX=40
YMIN=-20
YMAX=20
nc=100
CBMIN_DEN=0
CBMAX_DEN=2*${nc}
CBMIN_FLD=-5
CBMAX_FLD=5
YMIN_FLD=-500
YMAX_FLD=500

SIMULATION_FOLDERS=($(find . -name "${SIM_HEADER}*" -type d))


#LINEOUT_CREATION_REQUIRED=true
#PRINT_EFIELDS_LINEOUT=true
#PRINT_BFIELDS_LINEOUT=true
#PRINT_EDENSITY_LINEOUT=true
#PRINT_PDENSITY_LINEOUT=true
#PRINT_IDENSITY_LINEOUT=true



for sim in "${SIMULATION_FOLDERS[@]}"
do
 PREPLASMA_LENGTH=($(echo $sim | awk -F'_' '{print $2}'))
 RAMP_LENGTH=($(echo $sim | awk -F'_' '{print $6}'))
 DENSITY=($(echo $sim | awk -F'_' '{print $4}'))
 BULK_LENGTH=($(echo $sim | awk -F'_' '{print $8}'))
 CONT_LENGTH=($(echo $sim | awk -F'_' '{print $10}'))
 
 cd $sim

 for file in `ls -1 E?fout*.dat B?fout*.dat ?denout*.dat Hidnout*.dat` ; do
  basefile=`echo ${file%.*}`
  ${BINARY_READER} $basefile -dump_gnuplot -dump_lineoutx -dontask
 done


 for num in "${SIM_OUT_TO_BE_READ[@]}" ; do
  echo 'sto eseguendo il '$num
  lineout_file_E=Edenout${num}_lineout.txt
  lineout_file_P=Pdenout${num}_lineout.txt
  lineout_file_H=Hidnout${num}_lineout.txt
  lineout_file_out=EPHdenout${num}_lineout.$IMAGE_TYPE
  GNUPLOT_FILE='lineout.plt'
  rm -f ${GNUPLOT_FILE}
  touch ${GNUPLOT_FILE}
  chmod 775 ${GNUPLOT_FILE}
  printf "#!/gnuplot\n" >> ${GNUPLOT_FILE}
  printf "FILE_IN_1='${lineout_file_E}'\n" >> ${GNUPLOT_FILE}
  printf "FILE_IN_2='${lineout_file_P}'\n" >> ${GNUPLOT_FILE}
  printf "FILE_IN_3='${lineout_file_H}'\n" >> ${GNUPLOT_FILE}
  printf "FILE_OUT='${lineout_file_out}'\n" >> ${GNUPLOT_FILE}
  printf "set terminal ${IMAGE_TYPE} truecolor enhanced size ${SIZEX},${SIZEY} ${FONTSIZE}\n" >> ${GNUPLOT_FILE}
  printf "set output FILE_OUT\n" >> ${GNUPLOT_FILE}
  printf "nc=${nc}\n" >> ${GNUPLOT_FILE}
  printf "set xlabel 'x {/Symbol.ttf m}m' \n" >> ${GNUPLOT_FILE}
  printf "set ylabel 'n (n_c)' \n" >> ${GNUPLOT_FILE}
  #printf "set format y '10^{%%L}'\n" >> ${GNUPLOT_FILE}
  #printf "set logscale y\n" >> ${GNUPLOT_FILE}
  printf "set xrange[%s:%s]\n" "${XMIN}" "${XMAX}" >> ${GNUPLOT_FILE}
  printf "set yrange[-nc*1.2:nc*1.2]\n" >> ${GNUPLOT_FILE}
  printf "plot FILE_IN_1 u 1:(\$2*nc) w lines lt 1 lw 5 lc rgb 'red' t 'electron density',\\" >> ${GNUPLOT_FILE}
  printf "\n" >> ${GNUPLOT_FILE}
  printf "FILE_IN_2 u 1:(\$2*nc) w lines lt 1 lw 5 lc rgb 'blue' t 'proton density',\\" >> ${GNUPLOT_FILE}
  printf "\n" >> ${GNUPLOT_FILE}
  printf "FILE_IN_3 u 1:(\$2*nc) w lines lt 1 lw 5 lc rgb 'dark-green' t 'aluminium density'" >> ${GNUPLOT_FILE}
  $GNUPLOT ${GNUPLOT_FILE}
 done


 for file in `ls -1 *nout??.txt` ; do
  basefilename=${file%.*}
  FILEOUT=${basefilename}.${IMAGE_TYPE}
  GNUPLOT_FILE='grid_image.plt'
  rm -f ${GNUPLOT_FILE}
  touch ${GNUPLOT_FILE}
  chmod 775 ${GNUPLOT_FILE}
  printf "#!/gnuplot\n" >> ${GNUPLOT_FILE}
  printf "FILE_IN='${file}'\n" >> ${GNUPLOT_FILE}
  printf "FILE_OUT='${FILEOUT}'\n" >> ${GNUPLOT_FILE}
  printf "set terminal ${IMAGE_TYPE} truecolor enhanced size ${SIZEX},${SIZEY} ${FONTSIZE}\n" >> ${GNUPLOT_FILE}
  printf "set output FILE_OUT\n" >> ${GNUPLOT_FILE}
  printf "nc=${nc}\n" >> ${GNUPLOT_FILE}
  printf "set title '${basefilename}' \n" >> ${GNUPLOT_FILE}
  printf "set xlabel 'x ({/Symbol.ttf m}m)' \n" >> ${GNUPLOT_FILE}
  printf "set ylabel 'y ({/Symbol.ttf m}m)' \n" >> ${GNUPLOT_FILE}
  printf "set cblabel 'n/nc' \n" >> ${GNUPLOT_FILE}
  printf "set xrange[%s:%s]\n" "${XMIN}" "${XMAX}" >> ${GNUPLOT_FILE}
  printf "set yrange[%s:%s]\n" "${YMIN}" "${YMAX}" >> ${GNUPLOT_FILE}
  #printf "set cbrange[%s:%s]\n" "${CBMIN_DEN}" "${CBMAX_DEN}" >> ${GNUPLOT_FILE}
  #printf "set xtics 45\n" >> ${GNUPLOT_FILE}
  #printf "set ytics 1\n" >> ${GNUPLOT_FILE}
  #printf "set palette model RGB defined (0 'white', 0.25 'cyan', 0.625 'blue', 1 'black', 1.375 'red', 1.75 'yellow', 2 'white')\n" >> ${GNUPLOT_FILE}
  #printf "set palette model RGB defined (0 'black', 0.15 'blue', 0.85 'cyan', 1 'white', 1.15 'yellow', 1.85 'red', 2 'black')\n" >> ${GNUPLOT_FILE}
  printf "set palette model RGB defined (0 '#023858', 0.25 '#0570B0', 0.5 '#74A9CF', 0.85 '#D0D1E6',  1 '#FFF7FB', 1 '#FFF7EC', 1.25 '#FDD49E' ,1.5 '#FC8D59', 1.75 '#D7301F', 2 '#7F0000')\n" >> ${GNUPLOT_FILE}
  #printf "set palette model RGB defined (0 'blue', 1 'white', 2 'red')\n" >> ${GNUPLOT_FILE}
  #printf "set palette rgbformulae 7,5,15\n" >> ${GNUPLOT_FILE}
  #printf "set logscale cb\n" >> ${GNUPLOT_FILE}
  #printf "set format cb '%L10^{%L}'\n" >> ${GNUPLOT_FILE}
  printf "plot FILE_IN u 1:2:(abs(\$3*nc)) w image notitle\n" >> ${GNUPLOT_FILE}
  $GNUPLOT $GNUPLOT_FILE
 done


 for num in "${SIM_OUT_TO_BE_READ[@]}" ; do
  echo 'sto eseguendo il '$num
  lineout_file_X=Exfout${num}_lineout.txt
  lineout_file_Y=Eyfout${num}_lineout.txt
  lineout_file_Z=Bzfout${num}_lineout.txt
  lineout_file_out=ExyBzfout${num}_lineout.$IMAGE_TYPE
  GNUPLOT_FILE='lineout.plt'
  rm -f ${GNUPLOT_FILE}
  touch ${GNUPLOT_FILE}
  chmod 775 ${GNUPLOT_FILE}
  printf "#!/gnuplot\n" >> ${GNUPLOT_FILE}
  printf "FILE_IN_1='${lineout_file_X}'\n" >> ${GNUPLOT_FILE}
  printf "FILE_IN_2='${lineout_file_Y}'\n" >> ${GNUPLOT_FILE}
  printf "FILE_IN_3='${lineout_file_Z}'\n" >> ${GNUPLOT_FILE}
  printf "FILE_OUT='${lineout_file_out}'\n" >> ${GNUPLOT_FILE}
  printf "set terminal ${IMAGE_TYPE} truecolor enhanced size ${SIZEX},${SIZEY} ${FONTSIZE}\n" >> ${GNUPLOT_FILE}
  printf "set output FILE_OUT\n" >> ${GNUPLOT_FILE}
  printf "set xlabel 'x {/Symbol.ttf m}m' \n" >> ${GNUPLOT_FILE}
  printf "set ylabel 'E (TV/m)' \n" >> ${GNUPLOT_FILE}
  printf "set xrange[%s:%s]\n" "${XMIN}" "${XMAX}" >> ${GNUPLOT_FILE}
  printf "set yrange[%s:%s]\n" "${YMIN_FLD}" "${YMAX_FLD}" >> ${GNUPLOT_FILE}
  printf "plot FILE_IN_1 u 1:(\$2) w lines lt 1 lw 5 lc rgb 'red' t 'Ex',\\" >> ${GNUPLOT_FILE}
  printf "\n" >> ${GNUPLOT_FILE}
  printf "FILE_IN_2 u 1:(\$2) w lines lt 1 lw 5 lc rgb 'blue' t 'Ey',\\" >> ${GNUPLOT_FILE}
  printf "\n" >> ${GNUPLOT_FILE}
  printf "FILE_IN_3 u 1:(\$2) w lines lt 1 lw 5 lc rgb 'dark-green' t 'Bz'" >> ${GNUPLOT_FILE}
  $GNUPLOT ${GNUPLOT_FILE}
 done


 for file in `ls -1 *fout??.txt` ; do
  basefilename=${file%.*}
  FILEOUT=${basefilename}.${IMAGE_TYPE}
  GNUPLOT_FILE='grid_image.plt'
  rm -f ${GNUPLOT_FILE}
  touch ${GNUPLOT_FILE}
  chmod 775 ${GNUPLOT_FILE}
  printf "#!/gnuplot\n" >> ${GNUPLOT_FILE}
  printf "FILE_IN='${file}'\n" >> ${GNUPLOT_FILE}
  printf "FILE_OUT='${FILEOUT}'\n" >> ${GNUPLOT_FILE}
  printf "set terminal ${IMAGE_TYPE} truecolor enhanced size ${SIZEX},${SIZEY} ${FONTSIZE}\n" >> ${GNUPLOT_FILE}
  printf "set output FILE_OUT\n" >> ${GNUPLOT_FILE}
  printf "set title '${basefilename}' \n" >> ${GNUPLOT_FILE}
  printf "set xlabel 'x ({/Symbol.ttf m}m)' \n" >> ${GNUPLOT_FILE}
  printf "set ylabel 'y ({/Symbol.ttf m}m)' \n" >> ${GNUPLOT_FILE}
  printf "set cblabel 'TV/m' \n" >> ${GNUPLOT_FILE}
  printf "set xrange[%s:%s]\n" "${XMIN}" "${XMAX}" >> ${GNUPLOT_FILE}
  printf "set yrange[%s:%s]\n" "${YMIN}" "${YMAX}" >> ${GNUPLOT_FILE}
  printf "set cbrange[%s:%s]\n" "${CBMIN_FLD}" "${CBMAX_FLD}" >> ${GNUPLOT_FILE}
  #printf "set xtics 45\n" >> ${GNUPLOT_FILE}
  #printf "set ytics 1\n" >> ${GNUPLOT_FILE}
  #printf "set palette model RGB defined (0 'white', 0.25 'cyan', 0.625 'blue', 1 'black', 1.375 'red', 1.75 'yellow', 2 'white')\n" >> ${GNUPLOT_FILE}
  #printf "set palette model RGB defined (0 'black', 0.15 'blue', 0.85 'cyan', 1 'white', 1.15 'yellow', 1.85 'red', 2 'black')\n" >> ${GNUPLOT_FILE}
  printf "set palette model RGB defined (0 '#023858', 0.25 '#0570B0', 0.5 '#74A9CF', 0.85 '#D0D1E6',  1 '#FFF7FB', 1 '#FFF7EC', 1.25 '#FDD49E' ,1.5 '#FC8D59', 1.75 '#D7301F', 2 '#7F0000')\n" >> ${GNUPLOT_FILE}
  #printf "set palette model RGB defined (0 'blue', 1 'white', 2 'red')\n" >> ${GNUPLOT_FILE}
  #printf "set palette rgbformulae 7,5,15\n" >> ${GNUPLOT_FILE}
  #printf "set logscale cb\n" >> ${GNUPLOT_FILE}
  #printf "set format cb '%L10^{%L}'\n" >> ${GNUPLOT_FILE}
  printf "plot FILE_IN u 1:2:3 w image notitle\n" >> ${GNUPLOT_FILE}
  $GNUPLOT $GNUPLOT_FILE
 cd ..

done

