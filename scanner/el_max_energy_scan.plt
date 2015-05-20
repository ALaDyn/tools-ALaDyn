#!/gnuplot

set terminal png enhanced font ",15"

set xlabel 'Preplasma length ({/Symbol.ttf m}m)'
set ylabel 'Ramp length ({/Symbol.ttf m}m)'
set cblabel 'Maximum electron energy (MeV)'
#set size ratio -1
set xrange[0.5:3.5]
set yrange[0:1]
#set cbrange[2.0:4.0]
set palette rgbformulae 22,13,10
set datafile missing '-1'
set pm3d map
#set pm3d interpolate 0,0
#set pm3d interpolate 2,2
#splot FILE matrix

FILE_IN="energy_scan_Electron_2.txt"
FILE_OUT="max_energy_scan_Electron_2.0.png"
set output FILE_OUT
set title 'bulk length 2.0 {/Symbol.ttf m}m'
plot FILE_IN u 1:3:6 w points ps 24 pt 5 palette notitle

FILE_IN="energy_scan_Electron_3.txt"
FILE_OUT="max_energy_scan_Electron_3.0.png"
set output FILE_OUT
set title 'bulk length 3.0 {/Symbol.ttf m}m'
plot FILE_IN u 1:3:6 w points ps 24 pt 5 palette notitle

FILE_IN="energy_scan_Electron_4.txt"
FILE_OUT="max_energy_scan_Electron_4.0.png"
set output FILE_OUT
set title 'bulk length 4.0 {/Symbol.ttf m}m'
plot FILE_IN u 1:3:6 w points ps 24 pt 5 palette notitle

FILE_IN="energy_scan_Electron_5.txt"
FILE_OUT="max_energy_scan_Electron5.0.png"
set output FILE_OUT
set title 'bulk length 5.0 {/Symbol.ttf m}m'
plot FILE_IN u 1:3:6 w points ps 24 pt 5 palette notitle

FILE_IN="energy_scan_Electron_6.txt"
FILE_OUT="max_energy_scan_Electron_6.0.png"
set output FILE_OUT
set title 'bulk length 6.0 {/Symbol.ttf m}m'
plot FILE_IN u 1:3:6 w points ps 24 pt 5 palette notitle

FILE_IN="energy_scan_Electron_7.txt"
FILE_OUT="max_energy_scan_Electron_7.0.png"
set output FILE_OUT
set title 'bulk length 7.0 {/Symbol.ttf m}m'
plot FILE_IN u 1:3:6 w points ps 24 pt 5 palette notitle

FILE_IN="energy_scan_Electron_8.txt"
FILE_OUT="max_energy_scan_Electron_8.0.png"
set output FILE_OUT
set title 'bulk length 8.0 {/Symbol.ttf m}m'
plot FILE_IN u 1:3:6 w points ps 24 pt 5 palette notitle

FILE_IN="energy_scan_Electron_9.txt"
FILE_OUT="max_energy_scan_Electron_9.0.png"
set output FILE_OUT
set title 'bulk length 9.0 {/Symbol.ttf m}m'
plot FILE_IN u 1:3:6 w points ps 24 pt 5 palette notitle

FILE_IN="energy_scan_Electron_10.txt"
FILE_OUT="max_energy_scan_Electron_10.0.png"
set output FILE_OUT
set title 'bulk length 10.0 {/Symbol.ttf m}m'
plot FILE_IN u 1:3:6 w points ps 24 pt 5 palette notitle

