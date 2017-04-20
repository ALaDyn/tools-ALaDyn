#! /bin/bash

# per estrarre tutti i png generati nelle varie cartelle, e raccoglierli in un albero identico
# ma privo dei files .txt, .bin, .dat inutili, questo comando e' utilissimo:
#####
# find . -type f -name '*.png' | cpio -p -d -v ../target_folder
#####

BINARY_READER=$HOME/bin/leggi_bin
FIRST_SIM=0
LAST_SIM=100
GNUPLOT="$(which gnuplot)"
IMAGE_TYPE=png
SIZEX=1280
SIZEY=720
FONTSIZE=16

XMIN=30
XMAX=50
YMIN=-20
YMAX=20
nc=100
CBMIN_DEN=0
CBMAX_DEN="1.2*nc"
CBMIN_FLD=-5
CBMAX_FLD=5
 

 for file in $(find . -name "E?fout*.dat" -o -name "B?fout*.dat" -o -name "?denout*.dat" -o -name "H?dnout*.dat") ; do
  basefile=$(basename "${file}" .dat)
  ${BINARY_READER} "$basefile" -dump_gnuplot
 done


 for num in $(seq ${FIRST_SIM} ${LAST_SIM}) ; do
  file_E=Edenout${num}.txt
  file_P=Pdenout${num}.txt
  file_H=H1dnout${num}.txt
  file_Ex=Exfout${num}.txt
  file_Ey=Eyfout${num}.txt
  file_Bz=Bzfout${num}.txt
  file_outE=Edenout${num}.$IMAGE_TYPE
  file_outP=Pdenout${num}.$IMAGE_TYPE
  file_outH=H1denout${num}.$IMAGE_TYPE
  file_outEx=Exfout${num}.$IMAGE_TYPE
  file_outEy=Eyfout${num}.$IMAGE_TYPE
  file_outBz=Bzfout${num}.$IMAGE_TYPE
  GNUPLOT_FILE='grid.plt'

  echo 'sto plottando la griglia densita` elettroni #'"$num"
  rm -f ${GNUPLOT_FILE} && touch ${GNUPLOT_FILE} && chmod 775 ${GNUPLOT_FILE}
  {
    printf "#!/gnuplot\n"
    printf "FILE_IN='%s'\n" "${file_E}"
    printf "FILE_OUT='%s'\n" "${file_outE}"
    printf "set terminal %scairo size %s,%s font ',%s'\n" "${IMAGE_TYPE}" "${SIZEX}" "${SIZEY}" "${FONTSIZE}"
    printf "set output FILE_OUT\n"
    printf "nc=%s\n" "${nc}"
    printf "set xlabel 'x {/Symbol.ttf m}m' \n"
    printf "set ylabel 'y {/Symbol.ttf m}m' \n"
    printf "set cblabel 'n (n_c)' \n"
#    printf "set xrange[%s:%s]\n" "${XMIN}" "${XMAX}"
#    printf "set yrange[%s:%s]\n" "${YMIN}" "${YMAX}"
    printf "set cbrange[%s:%s]\n" "${CBMIN_DEN}" "${CBMAX_DEN}"
    printf "set palette model RGB defined (0 '#023858', 0.25 '#0570B0', 0.5 '#74A9CF', 0.85 '#D0D1E6',  1 '#FFF7FB', 1 '#FFF7EC', 1.25 '#FDD49E' ,1.5 '#FC8D59', 1.75 '#D7301F', 2 '#7F0000')\n"
    printf "plot FILE_IN u 1:2:(\$3*nc) w image notitle"
  } >> ${GNUPLOT_FILE}
  $GNUPLOT ${GNUPLOT_FILE}
  
  echo 'sto plottando la griglia densita` protoni #'"$num"
  rm -f ${GNUPLOT_FILE} && touch ${GNUPLOT_FILE} && chmod 775 ${GNUPLOT_FILE}
  {
    printf "#!/gnuplot\n"
    printf "FILE_IN='%s'\n" "${file_P}"
    printf "FILE_OUT='%s'\n" "${file_outP}"
    printf "set terminal %scairo size %s,%s font ',%s'\n" "${IMAGE_TYPE}" "${SIZEX}" "${SIZEY}" "${FONTSIZE}"
    printf "set output FILE_OUT\n"
    printf "nc=%s\n" "${nc}"
    printf "set xlabel 'x {/Symbol.ttf m}m' \n"
    printf "set ylabel 'y {/Symbol.ttf m}m' \n"
    printf "set cblabel 'n (n_c)' \n"
#    printf "set xrange[%s:%s]\n" "${XMIN}" "${XMAX}"
#    printf "set yrange[%s:%s]\n" "${YMIN}" "${YMAX}"
    printf "set logscale cb\n"
    printf "set cbrange[0.01:nc*0.1]\n"
    printf "set palette defined (0 '#FFFFCC',1 '#FFEDA0',2 '#FED976',3 '#FEB24C',4 '#FD8D3C',5 '#FC4E2A',6 '#E31A1C',7 '#B10026')\n"
    printf "plot FILE_IN u 1:2:(\$3*nc) w image notitle"
  } >> ${GNUPLOT_FILE}
  $GNUPLOT ${GNUPLOT_FILE}
  
  echo 'sto plottando la griglia densita` ioni #'"$num"
  rm -f ${GNUPLOT_FILE} && touch ${GNUPLOT_FILE} && chmod 775 ${GNUPLOT_FILE}
  {
    printf "#!/gnuplot\n"
    printf "FILE_IN='%s'\n" "${file_H}"
    printf "FILE_OUT='%s'\n" "${file_outH}"
    printf "set terminal %scairo size %s,%s font ',%s'\n" "${IMAGE_TYPE}" "${SIZEX}" "${SIZEY}" "${FONTSIZE}"
    printf "set output FILE_OUT\n"
    printf "nc=%s\n" "${nc}"
    printf "set xlabel 'x {/Symbol.ttf m}m' \n"
    printf "set ylabel 'y {/Symbol.ttf m}m' \n"
    printf "set cblabel 'n (n_c)' \n"
#    printf "set xrange[%s:%s]\n" "${XMIN}" "${XMAX}"
#    printf "set yrange[%s:%s]\n" "${YMIN}" "${YMAX}"
    printf "set cbrange[%s:%s]\n" "${CBMIN_DEN}" "${CBMAX_DEN}"
    printf "set palette defined (0 '#352a87',1 '#0363e1',2 '#1485d4',3 '#06a7c6',4 '#38b99e',5 '#92bf73',6 '#d9ba56',7 '#fcce2e',8 '#f9fb0e')\n"
    printf "plot FILE_IN u 1:2:(\$3*nc) w image notitle"
  } >> ${GNUPLOT_FILE}
  $GNUPLOT ${GNUPLOT_FILE}

  echo 'sto plottando la griglia Ex #'"$num"
  rm -f ${GNUPLOT_FILE} && touch ${GNUPLOT_FILE} && chmod 775 ${GNUPLOT_FILE}
  {
    printf "#!/gnuplot\n"
    printf "FILE_IN='%s'\n" "${file_Ex}"
    printf "FILE_OUT='%s'\n" "${file_outEx}"
    printf "set terminal %scairo size %s,%s font ',%s'\n" "${IMAGE_TYPE}" "${SIZEX}" "${SIZEY}" "${FONTSIZE}"
    printf "set output FILE_OUT\n"
    printf "nc=%s\n" "${nc}"
    printf "set xlabel 'x {/Symbol.ttf m}m' \n"
    printf "set ylabel 'y {/Symbol.ttf m}m' \n"
    printf "set cblabel 'Ex (TV/m)' \n"
#    printf "set xrange[%s:%s]\n" "${XMIN}" "${XMAX}"
#    printf "set yrange[%s:%s]\n" "${YMIN}" "${YMAX}"
    printf "set cbrange[%s:%s]\n" "${CBMIN_FLD}" "${CBMAX_FLD}"
    printf "set palette defined (0 '#352a87',1 '#0363e1',2 '#1485d4',3 '#06a7c6',4 '#38b99e',5 '#92bf73',6 '#d9ba56',7 '#fcce2e',8 '#f9fb0e')\n"
    printf "plot FILE_IN u 1:2:3 w image notitle"
  } >> ${GNUPLOT_FILE}
  $GNUPLOT ${GNUPLOT_FILE}

  echo 'sto plottando la griglia Ey #'"$num"
  rm -f ${GNUPLOT_FILE} && touch ${GNUPLOT_FILE} && chmod 775 ${GNUPLOT_FILE}
  {
    printf "#!/gnuplot\n"
    printf "FILE_IN='%s'\n" "${file_Ey}"
    printf "FILE_OUT='%s'\n" "${file_outEy}"
    printf "set terminal %scairo size %s,%s font ',%s'\n" "${IMAGE_TYPE}" "${SIZEX}" "${SIZEY}" "${FONTSIZE}"
    printf "set output FILE_OUT\n"
    printf "nc=%s\n" "${nc}"
    printf "set xlabel 'x {/Symbol.ttf m}m' \n"
    printf "set ylabel 'y {/Symbol.ttf m}m' \n"
    printf "set cblabel 'Ey (TV/m)' \n"
#    printf "set xrange[%s:%s]\n" "${XMIN}" "${XMAX}"
#    printf "set yrange[%s:%s]\n" "${YMIN}" "${YMAX}"
    printf "set cbrange[%s:%s]\n" "${CBMIN_FLD}" "${CBMAX_FLD}"
    printf "set palette defined (0 '#352a87',1 '#0363e1',2 '#1485d4',3 '#06a7c6',4 '#38b99e',5 '#92bf73',6 '#d9ba56',7 '#fcce2e',8 '#f9fb0e')\n"
    printf "plot FILE_IN u 1:2:3 w image notitle"
  } >> ${GNUPLOT_FILE}
  $GNUPLOT ${GNUPLOT_FILE}

  echo 'sto plottando la griglia Bz #'"$num"
  rm -f ${GNUPLOT_FILE} && touch ${GNUPLOT_FILE} && chmod 775 ${GNUPLOT_FILE}
  {
    printf "#!/gnuplot\n"
    printf "FILE_IN='%s'\n" "${file_Bz}"
    printf "FILE_OUT='%s'\n" "${file_outBz}"
    printf "set terminal %scairo size %s,%s font ',%s'\n" "${IMAGE_TYPE}" "${SIZEX}" "${SIZEY}" "${FONTSIZE}"
    printf "set output FILE_OUT\n"
    printf "nc=%s\n" "${nc}"
    printf "set xlabel 'x {/Symbol.ttf m}m' \n"
    printf "set ylabel 'y {/Symbol.ttf m}m' \n"
    printf "set cblabel 'Bz (TV/m)' \n"
#    printf "set xrange[%s:%s]\n" "${XMIN}" "${XMAX}"
#    printf "set yrange[%s:%s]\n" "${YMIN}" "${YMAX}"
    printf "set cbrange[%s:%s]\n" "${CBMIN_FLD}" "${CBMAX_FLD}"
    printf "set palette defined (0 '#352a87',1 '#0363e1',2 '#1485d4',3 '#06a7c6',4 '#38b99e',5 '#92bf73',6 '#d9ba56',7 '#fcce2e',8 '#f9fb0e')\n"
    printf "plot FILE_IN u 1:2:3 w image notitle"
  } >> ${GNUPLOT_FILE}
  $GNUPLOT ${GNUPLOT_FILE}
 done


