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
SPEC_DECODER=/galileo/home/userexternal/ssinigar/bin/leggi_diag
SPEC_VERSION=3
IMAGE_TYPE=png
SIZEX=960
SIZEY=540
FONTSIZE=14

#XMIN=0
#XMAX=70
#YMIN=-25
#YMAX=25
ricampionamento='322600*2'

SIMULATION_FOLDERS=($(find . -name "${SIM_HEADER}*" -type d))




for sim in "${SIMULATION_FOLDERS[@]}"
do
 PREPLASMA_LENGTH=($(echo $sim | awk -F'_' '{print $2}'))
 RAMP_LENGTH=($(echo $sim | awk -F'_' '{print $6}'))
 DENSITY=($(echo $sim | awk -F'_' '{print $4}'))
 BULK_LENGTH=($(echo $sim | awk -F'_' '{print $8}'))
 CONT_LENGTH=($(echo $sim | awk -F'_' '{print $10}'))
 
 cd $sim

 for file in `ls -1 ??pout??.dat` ; do
  basefilename=${file%.*}
  ${BINARY_READER} ${basefilename} -find_minmax -dontask
  ${BINARY_READER} ${basefilename} -readparamsfromfile ${basefilename}.extremes -do_binning -plot_xpx -plot_ethetaT -plot_espec -plot_thetaspec -nbin 120 -dontask
 done

 for file in `ls -1 spec??.dat` ; do
  ${SPEC_DECODER} file v${SPEC_VERSION}
 done


 for file in `ls -1 Pr*_Espec.txt` ; do
  echo 'elaboro file '$file
#  how to obtain spec file related to the _Espec one in a simple way, in order to print both in gnuplot? need to think about it!
#  SPEC_FILE=
  basefilename=${file%.*}
  FILEOUT=${basefilename}.${IMAGE_TYPE}
  GNUPLOT_FILE='Espec.plt'
  rm -f ${GNUPLOT_FILE}
  touch ${GNUPLOT_FILE}
  chmod 775 ${GNUPLOT_FILE}
  printf "#!/gnuplot\n" >> ${GNUPLOT_FILE}
  printf "FILE_IN='${file}'\n" >> ${GNUPLOT_FILE}
  printf "FILE_OUT='${FILEOUT}'\n" >> ${GNUPLOT_FILE}
  printf "set terminal ${IMAGE_TYPE} truecolor enhanced size ${SIZEX},${SIZEY} ${FONTSIZE}\n" >> ${GNUPLOT_FILE}
  printf "set output FILE_OUT\n" >> ${GNUPLOT_FILE}
  printf "set xlabel 'E (MeV)' \n" >> ${GNUPLOT_FILE}
  printf "set ylabel 'dN/dE (MeV^{-1})' \n" >> ${GNUPLOT_FILE}
#  printf "set xrange[%s:%s]\n" "${XMIN}" "${XMAX}" >> ${GNUPLOT_FILE}
  printf "set format y '10^{%%L}'\n" >> ${GNUPLOT_FILE}
  printf "set logscale y\n" >> ${GNUPLOT_FILE}
  printf "plot FILE_IN u (\$1+\$2)*0.5:(\$3*${ricampionamento}) with histeps lt 1 lc rgb 'blue' notitle\n" >> ${GNUPLOT_FILE}
  $GNUPLOT $GNUPLOT_FILE
 done

 cd ..

done

