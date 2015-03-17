#!/gnuplot

set terminal png enhanced font ",15"

set xlabel 'Preplasma length ({/Symbol.ttf m}m)'
set ylabel 'Ramp length ({/Symbol.ttf m}m)'
set cblabel 'Total proton energy (J)'
set size ratio -1
set xrange[-0.2:4.2]
set yrange[-0.2:2.2]
#set cbrange[0.18:0.32]
set palette rgbformulae 22,13,10
set datafile missing '-1'
set pm3d map
#set pm3d interpolate 0,0
#set pm3d interpolate 2,2
#splot FILE matrix

FILE_IN="energy_scan_1.txt"
FILE_OUT="tot_energy_scan_1.0.png"
set output FILE_OUT
set title 'initial density n=1.0n_c'
plot FILE_IN u 1:3:7 w points ps 6.2 pt 5 palette notitle

FILE_IN="energy_scan_2.txt"
FILE_OUT="tot_energy_scan_2.0.png"
set output FILE_OUT
set title 'initial density n=2.0n_c'
plot FILE_IN u 1:3:7 w points ps 6.2 pt 5 palette notitle
