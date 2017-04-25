#! /bin/bash

# per estrarre tutti i png generati nelle varie cartelle, e raccoglierli in un albero identico
# ma privo dei files .txt, .bin, .dat inutili, questo comando e' utilissimo:
#####
# find . -type f -name '*.png' | cpio -p -d -v ../target_folder
#####

ENABLE_CNAF_JOB_SUBMISSION=true
BINARY_READER=$HOME/bin/leggi_bin
FIRST_SIM=0
LAST_SIM=99
GNUPLOT="$(which gnuplot)"
IMAGE_TYPE=png
SIZEX=1280
SIZEY=720
SIZEX_COMPOSITION=1440
SIZEY_COMPOSITION=720
FONTSIZE=16

XMIN=0
XMAX=120
YMIN=-60
YMAX=60
nc=100
CBMIN_DEN=0.1
CBMAX_DEN=100
CBMIN_FLD=-2
CBMAX_FLD=2
 
mkdir -p images

for folder in $(seq -f "%04g" ${FIRST_SIM} ${LAST_SIM}) ; do

 cd "$folder" || exit
 num=${folder:2:3}

 #for file in $(find . -name "E?fout${num}.bin" -o -name "B?fout${num}.bin" -o -name "?denout${num}.bin" -o -name "H?dnout${num}.bin") ; do
 for file in $(find . -name "Exfout${num}.bin" -o -name "Edenout${num}.bin") ; do
  basefile=$(basename "${file}" .bin)
  if [ ! -s "${basefile}.txt" ] ; then ${BINARY_READER} "$basefile" -dump_gnuplot ; fi
 done

 file_E=Edenout${num}.txt
 file_P=Pdenout${num}.txt
 file_H=H1dnout${num}.txt
 file_Ex=Exfout${num}.txt
 file_Ey=Eyfout${num}.txt
 file_Bz=Bzfout${num}.txt
 gnuplot_file_E=Edenout${num}.plt
 gnuplot_file_P=Pdenout${num}.plt
 gnuplot_file_H=H1dnout${num}.plt
 gnuplot_file_Ex=Exfout${num}.plt
 gnuplot_file_Ey=Eyfout${num}.plt
 gnuplot_file_Bz=Bzfout${num}.plt
 gnuplot_file_ExfEden=Bzfout${num}.plt
 JOB_GNUPLOT=plt${num}.cmd
 file_outE=Edenout${num}.${IMAGE_TYPE}
 file_outP=Pdenout${num}.${IMAGE_TYPE}
 file_outH=H1denout${num}.${IMAGE_TYPE}
 file_outEx=Exfout${num}.${IMAGE_TYPE}
 file_outEy=Eyfout${num}.${IMAGE_TYPE}
 file_outBz=Bzfout${num}.${IMAGE_TYPE}
 file_outExfEden=ExfEden${num}.${IMAGE_TYPE}

 GNUPLOT_FILE=${gnuplot_file_E}
 echo 'sto plottando la griglia densita` elettroni #'"$num"
 rm -f "${GNUPLOT_FILE}" ; touch "${GNUPLOT_FILE}" ; chmod 775 "${GNUPLOT_FILE}"
 {
   printf "#!%s\n" "${GNUPLOT}"
   printf "FILE_IN='%s'\n" "${file_E}"
   printf "FILE_OUT='%s'\n" "${file_outE}"
   printf "set terminal %scairo size %s,%s font ',%s'\n" "${IMAGE_TYPE}" "${SIZEX}" "${SIZEY}" "${FONTSIZE}"
   printf "set output FILE_OUT\n"
   printf "nc=%s\n" "${nc}"
   printf "set xlabel 'x {/Symbol.ttf m}m' \n"
   printf "set ylabel 'y {/Symbol.ttf m}m' \n"
   printf "set cblabel 'n (n_c)' \n"
#   printf "set xrange[%s:%s]\n" "${XMIN}" "${XMAX}"
#   printf "set yrange[%s:%s]\n" "${YMIN}" "${YMAX}"
   printf "set cbrange[%s:%s]\n" "${CBMIN_DEN}" "${CBMAX_DEN}"
   printf "set palette model RGB defined (0 '#023858', 0.25 '#0570B0', 0.5 '#74A9CF', 0.85 '#D0D1E6',  1 '#FFF7FB', 1 '#FFF7EC', 1.25 '#FDD49E' ,1.5 '#FC8D59', 1.75 '#D7301F', 2 '#7F0000')\n"
   printf "plot FILE_IN u 1:2:(\$3*nc) w image notitle"
 } >> "${GNUPLOT_FILE}"
 #if [ ! -s "${file_outE}" ] ; then $GNUPLOT "${GNUPLOT_FILE}" ; fi
  
 GNUPLOT_FILE=${gnuplot_file_P}
 echo 'sto plottando la griglia densita` protoni #'"$num"
 rm -f "${GNUPLOT_FILE}" ; touch "${GNUPLOT_FILE}" ; chmod 775 "${GNUPLOT_FILE}"
 {
   printf "#!%s\n" "${GNUPLOT}"
   printf "FILE_IN='%s'\n" "${file_P}"
   printf "FILE_OUT='%s'\n" "${file_outP}"
   printf "set terminal %scairo size %s,%s font ',%s'\n" "${IMAGE_TYPE}" "${SIZEX}" "${SIZEY}" "${FONTSIZE}"
   printf "set output FILE_OUT\n"
   printf "nc=%s\n" "${nc}"
   printf "set xlabel 'x {/Symbol.ttf m}m' \n"
   printf "set ylabel 'y {/Symbol.ttf m}m' \n"
   printf "set cblabel 'n (n_c)' \n"
#   printf "set xrange[%s:%s]\n" "${XMIN}" "${XMAX}"
#   printf "set yrange[%s:%s]\n" "${YMIN}" "${YMAX}"
   printf "set logscale cb\n"
   printf "set cbrange[0.01:nc*0.1]\n"
   printf "set palette defined (0 '#FFFFCC',1 '#FFEDA0',2 '#FED976',3 '#FEB24C',4 '#FD8D3C',5 '#FC4E2A',6 '#E31A1C',7 '#B10026')\n"
   printf "plot FILE_IN u 1:2:(\$3*nc) w image notitle"
 } >> "${GNUPLOT_FILE}"
 #if [ ! -s "${file_outP}" ] ; then $GNUPLOT "${GNUPLOT_FILE}" ; fi
  
 GNUPLOT_FILE=${gnuplot_file_H}
 echo 'sto plottando la griglia densita` ioni #'"$num"
 rm -f "${GNUPLOT_FILE}" ; touch "${GNUPLOT_FILE}" ; chmod 775 "${GNUPLOT_FILE}"
 {
   printf "#!%s\n" "${GNUPLOT}"
   printf "FILE_IN='%s'\n" "${file_H}"
   printf "FILE_OUT='%s'\n" "${file_outH}"
   printf "set terminal %scairo size %s,%s font ',%s'\n" "${IMAGE_TYPE}" "${SIZEX}" "${SIZEY}" "${FONTSIZE}"
   printf "set output FILE_OUT\n"
   printf "nc=%s\n" "${nc}"
   printf "set xlabel 'x {/Symbol.ttf m}m' \n"
   printf "set ylabel 'y {/Symbol.ttf m}m' \n"
   printf "set cblabel 'n (n_c)' \n"
#   printf "set xrange[%s:%s]\n" "${XMIN}" "${XMAX}"
#   printf "set yrange[%s:%s]\n" "${YMIN}" "${YMAX}"
   printf "set cbrange[%s:%s]\n" "${CBMIN_DEN}" "${CBMAX_DEN}"
   printf "set palette defined (0 '#352a87',1 '#0363e1',2 '#1485d4',3 '#06a7c6',4 '#38b99e',5 '#92bf73',6 '#d9ba56',7 '#fcce2e',8 '#f9fb0e')\n"
   printf "plot FILE_IN u 1:2:(\$3*nc) w image notitle"
 } >> "${GNUPLOT_FILE}"
 #if [ ! -s "${file_outH}" ] ; then $GNUPLOT "${GNUPLOT_FILE}" ; fi

 GNUPLOT_FILE=${gnuplot_file_Ex}
 echo 'sto plottando la griglia Ex #'"$num"
 rm -f "${GNUPLOT_FILE}" ; touch "${GNUPLOT_FILE}" ; chmod 775 "${GNUPLOT_FILE}"
 {
   printf "#!%s\n" "${GNUPLOT}"
   printf "FILE_IN='%s'\n" "${file_Ex}"
   printf "FILE_OUT='%s'\n" "${file_outEx}"
   printf "set terminal %scairo size %s,%s font ',%s'\n" "${IMAGE_TYPE}" "${SIZEX}" "${SIZEY}" "${FONTSIZE}"
   printf "set output FILE_OUT\n"
   printf "set xlabel 'x {/Symbol.ttf m}m' \n"
   printf "set ylabel 'y {/Symbol.ttf m}m' \n"
   printf "set cblabel 'Ex (TV/m)' \n"
#   printf "set xrange[%s:%s]\n" "${XMIN}" "${XMAX}"
#   printf "set yrange[%s:%s]\n" "${YMIN}" "${YMAX}"
   printf "set cbrange[%s:%s]\n" "${CBMIN_FLD}" "${CBMAX_FLD}"
   printf "set palette defined (0 '#352a87',1 '#0363e1',2 '#1485d4',3 '#06a7c6',4 '#38b99e',5 '#92bf73',6 '#d9ba56',7 '#fcce2e',8 '#f9fb0e')\n"
   printf "plot FILE_IN u 1:2:3 w image notitle"
 } >> "${GNUPLOT_FILE}"
 #if [ ! -s "${file_outEx}" ] ; then $GNUPLOT "${GNUPLOT_FILE}" ; fi

 GNUPLOT_FILE=${gnuplot_file_Ey}
 echo 'sto plottando la griglia Ey #'"$num"
 rm -f "${GNUPLOT_FILE}" ; touch "${GNUPLOT_FILE}" ; chmod 775 "${GNUPLOT_FILE}"
 {
   printf "#!%s\n" "${GNUPLOT}"
   printf "FILE_IN='%s'\n" "${file_Ey}"
   printf "FILE_OUT='%s'\n" "${file_outEy}"
   printf "set terminal %scairo size %s,%s font ',%s'\n" "${IMAGE_TYPE}" "${SIZEX}" "${SIZEY}" "${FONTSIZE}"
   printf "set output FILE_OUT\n"
   printf "set xlabel 'x {/Symbol.ttf m}m' \n"
   printf "set ylabel 'y {/Symbol.ttf m}m' \n"
   printf "set cblabel 'Ey (TV/m)' \n"
#   printf "set xrange[%s:%s]\n" "${XMIN}" "${XMAX}"
#   printf "set yrange[%s:%s]\n" "${YMIN}" "${YMAX}"
   printf "set cbrange[%s:%s]\n" "${CBMIN_FLD}" "${CBMAX_FLD}"
   printf "set palette defined (0 '#352a87',1 '#0363e1',2 '#1485d4',3 '#06a7c6',4 '#38b99e',5 '#92bf73',6 '#d9ba56',7 '#fcce2e',8 '#f9fb0e')\n"
   printf "plot FILE_IN u 1:2:3 w image notitle"
 } >> "${GNUPLOT_FILE}"
 #if [ ! -s "${file_outEy}" ] ; then $GNUPLOT "${GNUPLOT_FILE}" ; fi

 GNUPLOT_FILE=${gnuplot_file_Bz}
 echo 'sto plottando la griglia Bz #'"$num"
 rm -f "${GNUPLOT_FILE}" ; touch "${GNUPLOT_FILE}" ; chmod 775 "${GNUPLOT_FILE}"
 {
   printf "#!%s\n" "${GNUPLOT}"
   printf "FILE_IN='%s'\n" "${file_Bz}"
   printf "FILE_OUT='%s'\n" "${file_outBz}"
   printf "set terminal %scairo size %s,%s font ',%s'\n" "${IMAGE_TYPE}" "${SIZEX}" "${SIZEY}" "${FONTSIZE}"
   printf "set output FILE_OUT\n"
   printf "set xlabel 'x {/Symbol.ttf m}m' \n"
   printf "set ylabel 'y {/Symbol.ttf m}m' \n"
   printf "set cblabel 'Bz (TV/m)' \n"
#   printf "set xrange[%s:%s]\n" "${XMIN}" "${XMAX}"
#   printf "set yrange[%s:%s]\n" "${YMIN}" "${YMAX}"
   printf "set cbrange[%s:%s]\n" "${CBMIN_FLD}" "${CBMAX_FLD}"
   printf "set palette defined (0 '#352a87',1 '#0363e1',2 '#1485d4',3 '#06a7c6',4 '#38b99e',5 '#92bf73',6 '#d9ba56',7 '#fcce2e',8 '#f9fb0e')\n"
   printf "plot FILE_IN u 1:2:3 w image notitle"
 } >> "${GNUPLOT_FILE}"
 #if [ ! -s "${file_outBz}" ] ; then $GNUPLOT "${GNUPLOT_FILE}" ; fi

 GNUPLOT_FILE=${gnuplot_file_ExfEden}
 echo 'sto plottando la composizione Ex-Eden #'"$num"
 rm -f "${GNUPLOT_FILE}" ; touch "${GNUPLOT_FILE}" ; chmod 775 "${GNUPLOT_FILE}"
 {
   printf "#!/usr/bin/gnuplot\n"
   printf "FILE_IN_EDEN='%s'\n" "${file_E}"
   printf "FILE_IN_EXF='%s'\n" "${file_Ex}"
   printf "FILE_OUT='%s'\n" "${file_outExfEden}"
   printf "set terminal %scairo size %s,%s font ',%s'\n" "${IMAGE_TYPE}" "${SIZEX_COMPOSITION}" "${SIZEY_COMPOSITION}" "${FONTSIZE}"
   printf "set output FILE_OUT\n"
   printf "set multiplot layout 1, 2 title 'Laser-plasma interaction'\n"
   printf "set title 'E_x field'\n"
   printf "set ylabel 'x ({/Symbol.ttf m}m)'\n"
   printf "set xlabel 'y ({/Symbol.ttf m}m)'\n"
   printf "set cblabel 'E_x (TV/m)'\n"
   printf "set cbrange[%s:%s]\n" "${CBMIN_FLD}" "${CBMAX_FLD}"
   printf "set xrange[%s:%s]\n" "${YMIN}" "${YMAX}"
   printf "set yrange[%s:%s]\n" "${XMIN}" "${XMAX}"
   printf "set palette model RGB defined  (0 '#023858',\\"
   printf "\n"
   printf "                      0.25 '#0570B0',\\"
   printf "\n"
   printf "                      0.5 '#74A9CF',\\"
   printf "\n"
   printf "                      0.85 '#D0D1E6',\\"
   printf "\n"
   printf "                      1 '#FFF7FB',\\"
   printf "\n"
   printf "                      1 '#FFF7EC',\\"
   printf "\n"
   printf "                      1.25 '#FDD49E',\\"
   printf "\n"
   printf "                      1.5 '#FC8D59',\\"
   printf "\n"
   printf "                      1.75 '#D7301F',\\"
   printf "\n"
   printf "                      2 '#7F0000')\n"
   printf "plot FILE_IN_EXF u 2:1:3 w image notitle\n"
   printf "set title 'e^- density'\n"
   printf "nc=%s\n" "${nc}"
   printf "set ylabel 'x ({/Symbol.ttf m}m)'\n"
   printf "set xlabel 'y ({/Symbol.ttf m}m)'\n"
   printf "set cblabel 'n (n_c)'\n"
   printf "set cbrange[%s:%s]\n" "${CBMIN_DEN}" "${CBMAX_DEN}"
   printf "set xrange[%s:%s]\n" "${YMIN}" "${YMAX}"
   printf "set yrange[%s:%s]\n" "${XMIN}" "${XMAX}"
   printf "set logscale cb\n"
   printf "set palette defined    (0 '#352a87',\\"
   printf "\n"
   printf "                1 '#0363e1',\\"
   printf "\n"
   printf "                2 '#1485d4',\\"
   printf "\n"
   printf "                3 '#06a7c6',\\"
   printf "\n"
   printf "                4 '#38b99e',\\"
   printf "\n"
   printf "                5 '#92bf73',\\"
   printf "\n"
   printf "                6 '#d9ba56',\\"
   printf "\n"
   printf "                7 '#fcce2e',\\"
   printf "\n"
   printf "                8 '#f9fb0e')\n"
   printf "plot FILE_IN_EDEN u 2:1:(\$3*nc) w image notitle\n"
   printf "unset multiplot\n"
 } >> "${GNUPLOT_FILE}"
 if [ ! -s "${file_outExfEden}" ] ; then 
 	if [ ${ENABLE_CNAF_JOB_SUBMISSION} ] ; then 
 		rm -f "${JOB_GNUPLOT}" ; touch "${JOB_GNUPLOT}" ; chmod 775 "${JOB_GNUPLOT}"
		{
			printf "#BSUB -J impi${num}\n"
			printf "#BSUB -o %%J.out\n"
			printf "#BSUB -e %%J.err\n"
			printf "#BSUB -q hpc_short\n"
			printf "#BSUB -n 16\n"
			printf "#BSUB -R \"span[ptile=16]\"\n"
			printf "module load compilers/gcc-4.9.0\n"
			printf "module load boost_1_56_0_gcc4_9_0\n"
			printf "cd \"$folder\"\n"
			printf "${GNUPLOT} ${GNUPLOT_FILE}\n"
		} >> "${JOB_GNUPLOT}"
 		bsub < "${JOB_GNUPLOT}"
 	else
 		$GNUPLOT "${GNUPLOT_FILE}"
 	fi
 fi

 #cp ./*.png ../images/

 cd .. || exit

done


