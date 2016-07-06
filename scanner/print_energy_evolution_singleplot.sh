#! /bin/bash
module load gnu

SIM_HEADER="0.0_pre_"
LEGGI_DIAG=$HOME/bin/leggi_diagspec.exe
DIAG_VERSION=4
DIAG_TIME_TO_BE_READ="05"
TERMINAL_TYPE=pngcairo
SIZEX=3840
SIZEY=2160
FONTSIZE=25
SIMULATION_FOLDERS=($(find . -name "${SIM_HEADER}*" -type d))
PROTON_EMAX_COLUMN=26

for sim in "${SIMULATION_FOLDERS[@]}"
do
 ANGLE=($(echo $sim | awk -F'_' '{print $1}'))
 PREPLASMA_LENGTH=($(echo $sim | awk -F'_' '{print $3}'))
 RAMP_LENGTH=($(echo $sim | awk -F'_' '{print $7}'))
 DENSITY=($(echo $sim | awk -F'_' '{print $5}'))
 BULK_LENGTH=($(echo $sim | awk -F'_' '{print $9}'))
 CONT_LENGTH=($(echo $sim | awk -F'_' '{print $11}'))
 
 cd $sim/diagnostics

 if [ -f diag${DIAG_TIME_TO_BE_READ}.dat ]
 then
  filename=`ls -1 diag${DIAG_TIME_TO_BE_READ}.dat`
  rm *.txt
  ${LEGGI_DIAG} $filename v${DIAG_VERSION}
 fi

 cd ../..

done

GNUPLOT_FILE=diag.plt
rm -f $GNUPLOT_FILE
touch $GNUPLOT_FILE
chmod 755 ${GNUPLOT_FILE}

printf "#!/gnuplot\n" > ${GNUPLOT_FILE}
printf "FILE_OUT='diags.png'\n" "${PREPLASMA_LENGTH}" "${RAMP_LENGTH}" "${DENSITY}" "${BULK_LENGTH}" "${CONT_LENGTH}" >> ${GNUPLOT_FILE}
printf "c=10/3\n" >> ${GNUPLOT_FILE}
printf "set terminal ${TERMINAL_TYPE} size ${SIZEX},${SIZEY} font \",${FONTSIZE}\"\n" >> ${GNUPLOT_FILE}
printf "set output FILE_OUT\n" >> ${GNUPLOT_FILE}
printf "set title 'E_{max} for different preplasma configurations'\n" "${PREPLASMA_LENGTH}" "${RAMP_LENGTH}" "${DENSITY}" "${BULK_LENGTH}" "${CONT_LENGTH}" >> ${GNUPLOT_FILE}
printf "set xlabel 't (fs)'\n" >> ${GNUPLOT_FILE}
printf "set ylabel 'E_{max} (MeV)'\n" >> ${GNUPLOT_FILE}
printf "set style line 1 lt 1 lw 5 lc rgb '#000000'\n" >> ${GNUPLOT_FILE}
printf "set style line 2 lt 1 lw 5 lc rgb '#FFFF00'\n" >> ${GNUPLOT_FILE}
printf "set style line 3 lt 1 lw 5 lc rgb '#1CE6FF'\n" >> ${GNUPLOT_FILE}
printf "set style line 4 lt 1 lw 5 lc rgb '#FF34FF'\n" >> ${GNUPLOT_FILE}
printf "set style line 5 lt 1 lw 5 lc rgb '#FF4A46'\n" >> ${GNUPLOT_FILE}
printf "set style line 6 lt 1 lw 5 lc rgb '#008941'\n" >> ${GNUPLOT_FILE}
printf "set style line 7 lt 1 lw 5 lc rgb '#006FA6'\n" >> ${GNUPLOT_FILE}
printf "set style line 8 lt 1 lw 5 lc rgb '#A30059'\n" >> ${GNUPLOT_FILE}
printf "set style line 9 lt 1 lw 5 lc rgb '#FFDBE5'\n" >> ${GNUPLOT_FILE}
printf "set style line 10 lt 1 lw 5 lc rgb '#7A4900'\n" >> ${GNUPLOT_FILE}
printf "set style line 11 lt 1 lw 5 lc rgb '#0000A6'\n" >> ${GNUPLOT_FILE}
printf "set style line 12 lt 1 lw 5 lc rgb '#63FFAC'\n" >> ${GNUPLOT_FILE}
printf "set style line 13 lt 1 lw 5 lc rgb '#B79762'\n" >> ${GNUPLOT_FILE}
printf "set style line 14 lt 1 lw 5 lc rgb '#004D43'\n" >> ${GNUPLOT_FILE}
printf "set style line 15 lt 1 lw 5 lc rgb '#8FB0FF'\n" >> ${GNUPLOT_FILE}
printf "set style line 16 lt 1 lw 5 lc rgb '#997D87'\n" >> ${GNUPLOT_FILE}
printf "set style line 17 lt 1 lw 5 lc rgb '#5A0007'\n" >> ${GNUPLOT_FILE}
printf "set style line 18 lt 1 lw 5 lc rgb '#809693'\n" >> ${GNUPLOT_FILE}
printf "set style line 19 lt 1 lw 5 lc rgb '#FEFFE6'\n" >> ${GNUPLOT_FILE}
printf "set style line 20 lt 1 lw 5 lc rgb '#1B4400'\n" >> ${GNUPLOT_FILE}
printf "set style line 21 lt 1 lw 5 lc rgb '#4FC601'\n" >> ${GNUPLOT_FILE}
printf "set style line 22 lt 1 lw 5 lc rgb '#3B5DFF'\n" >> ${GNUPLOT_FILE}
printf "set style line 23 lt 1 lw 5 lc rgb '#4A3B53'\n" >> ${GNUPLOT_FILE}
printf "set style line 24 lt 1 lw 5 lc rgb '#FF2F80'\n" >> ${GNUPLOT_FILE}
printf "set style line 25 lt 1 lw 5 lc rgb '#61615A'\n" >> ${GNUPLOT_FILE}
printf "set style line 26 lt 1 lw 5 lc rgb '#BA0900'\n" >> ${GNUPLOT_FILE}
printf "set style line 27 lt 1 lw 5 lc rgb '#6B7900'\n" >> ${GNUPLOT_FILE}
printf "set style line 28 lt 1 lw 5 lc rgb '#00C2A0'\n" >> ${GNUPLOT_FILE}
printf "set style line 29 lt 1 lw 5 lc rgb '#FFAA92'\n" >> ${GNUPLOT_FILE}
printf "set style line 30 lt 1 lw 5 lc rgb '#FF90C9'\n" >> ${GNUPLOT_FILE}
printf "set style line 31 lt 1 lw 5 lc rgb '#B903AA'\n" >> ${GNUPLOT_FILE}
printf "set style line 32 lt 1 lw 5 lc rgb '#D16100'\n" >> ${GNUPLOT_FILE}
printf "set style line 33 lt 1 lw 5 lc rgb '#DDEFFF'\n" >> ${GNUPLOT_FILE}
printf "set style line 34 lt 1 lw 5 lc rgb '#000035'\n" >> ${GNUPLOT_FILE}
printf "set style line 35 lt 1 lw 5 lc rgb '#7B4F4B'\n" >> ${GNUPLOT_FILE}
printf "set style line 36 lt 1 lw 5 lc rgb '#A1C299'\n" >> ${GNUPLOT_FILE}
printf "set style line 37 lt 1 lw 5 lc rgb '#300018'\n" >> ${GNUPLOT_FILE}
printf "set style line 38 lt 1 lw 5 lc rgb '#0AA6D8'\n" >> ${GNUPLOT_FILE}
printf "set style line 39 lt 1 lw 5 lc rgb '#013349'\n" >> ${GNUPLOT_FILE}
printf "set style line 40 lt 1 lw 5 lc rgb '#00846F'\n" >> ${GNUPLOT_FILE}
printf "set style line 41 lt 1 lw 5 lc rgb '#372101'\n" >> ${GNUPLOT_FILE}
printf "set style line 42 lt 1 lw 5 lc rgb '#FFB500'\n" >> ${GNUPLOT_FILE}
printf "set style line 43 lt 1 lw 5 lc rgb '#C2FFED'\n" >> ${GNUPLOT_FILE}
printf "set style line 44 lt 1 lw 5 lc rgb '#A079BF'\n" >> ${GNUPLOT_FILE}
printf "set style line 45 lt 1 lw 5 lc rgb '#CC0744'\n" >> ${GNUPLOT_FILE}
printf "set style line 46 lt 1 lw 5 lc rgb '#C0B9B2'\n" >> ${GNUPLOT_FILE}
printf "set style line 47 lt 1 lw 5 lc rgb '#C2FF99'\n" >> ${GNUPLOT_FILE}
printf "set style line 48 lt 1 lw 5 lc rgb '#001E09'\n" >> ${GNUPLOT_FILE}
printf "set style line 49 lt 1 lw 5 lc rgb '#00489C'\n" >> ${GNUPLOT_FILE}
printf "set style line 50 lt 1 lw 5 lc rgb '#6F0062'\n" >> ${GNUPLOT_FILE}
printf "set style line 51 lt 1 lw 5 lc rgb '#0CBD66'\n" >> ${GNUPLOT_FILE}
printf "set style line 52 lt 1 lw 5 lc rgb '#EEC3FF'\n" >> ${GNUPLOT_FILE}
printf "set style line 53 lt 1 lw 5 lc rgb '#456D75'\n" >> ${GNUPLOT_FILE}
printf "set style line 54 lt 1 lw 5 lc rgb '#B77B68'\n" >> ${GNUPLOT_FILE}
printf "set style line 55 lt 1 lw 5 lc rgb '#7A87A1'\n" >> ${GNUPLOT_FILE}
printf "set style line 56 lt 1 lw 5 lc rgb '#788D66'\n" >> ${GNUPLOT_FILE}
printf "set style line 57 lt 1 lw 5 lc rgb '#885578'\n" >> ${GNUPLOT_FILE}
printf "set style line 58 lt 1 lw 5 lc rgb '#FAD09F'\n" >> ${GNUPLOT_FILE}
printf "set style line 59 lt 1 lw 5 lc rgb '#FF8A9A'\n" >> ${GNUPLOT_FILE}
printf "set style line 60 lt 1 lw 5 lc rgb '#D157A0'\n" >> ${GNUPLOT_FILE}
printf "set style line 61 lt 1 lw 5 lc rgb '#BEC459'\n" >> ${GNUPLOT_FILE}
printf "set style line 62 lt 1 lw 5 lc rgb '#456648'\n" >> ${GNUPLOT_FILE}
printf "set style line 63 lt 1 lw 5 lc rgb '#0086ED'\n" >> ${GNUPLOT_FILE}
printf "set style line 64 lt 1 lw 5 lc rgb '#886F4C'\n" >> ${GNUPLOT_FILE}
printf "set style line 65 lt 1 lw 5 lc rgb '#34362D'\n" >> ${GNUPLOT_FILE}
printf "set style line 66 lt 1 lw 5 lc rgb '#B4A8BD'\n" >> ${GNUPLOT_FILE}
printf "set style line 67 lt 1 lw 5 lc rgb '#00A6AA'\n" >> ${GNUPLOT_FILE}
printf "set style line 68 lt 1 lw 5 lc rgb '#452C2C'\n" >> ${GNUPLOT_FILE}
printf "set style line 69 lt 1 lw 5 lc rgb '#636375'\n" >> ${GNUPLOT_FILE}
printf "set style line 70 lt 1 lw 5 lc rgb '#A3C8C9'\n" >> ${GNUPLOT_FILE}
printf "set style line 71 lt 1 lw 5 lc rgb '#FF913F'\n" >> ${GNUPLOT_FILE}
printf "set style line 72 lt 1 lw 5 lc rgb '#938A81'\n" >> ${GNUPLOT_FILE}
printf "set style line 73 lt 1 lw 5 lc rgb '#575329'\n" >> ${GNUPLOT_FILE}
printf "set style line 74 lt 1 lw 5 lc rgb '#00FECF'\n" >> ${GNUPLOT_FILE}
printf "set style line 75 lt 1 lw 5 lc rgb '#B05B6F'\n" >> ${GNUPLOT_FILE}
printf "set style line 76 lt 1 lw 5 lc rgb '#8CD0FF'\n" >> ${GNUPLOT_FILE}
printf "set style line 77 lt 1 lw 5 lc rgb '#3B9700'\n" >> ${GNUPLOT_FILE}
printf "set style line 78 lt 1 lw 5 lc rgb '#04F757'\n" >> ${GNUPLOT_FILE}
printf "set style line 79 lt 1 lw 5 lc rgb '#C8A1A1'\n" >> ${GNUPLOT_FILE}
printf "set style line 80 lt 1 lw 5 lc rgb '#1E6E00'\n" >> ${GNUPLOT_FILE}
printf "set style line 81 lt 1 lw 5 lc rgb '#7900D7'\n" >> ${GNUPLOT_FILE}
printf "set style line 82 lt 1 lw 5 lc rgb '#A77500'\n" >> ${GNUPLOT_FILE}
printf "set style line 83 lt 1 lw 5 lc rgb '#6367A9'\n" >> ${GNUPLOT_FILE}
printf "set style line 84 lt 1 lw 5 lc rgb '#A05837'\n" >> ${GNUPLOT_FILE}
printf "set style line 85 lt 1 lw 5 lc rgb '#6B002C'\n" >> ${GNUPLOT_FILE}
printf "set style line 86 lt 1 lw 5 lc rgb '#772600'\n" >> ${GNUPLOT_FILE}
printf "set style line 87 lt 1 lw 5 lc rgb '#D790FF'\n" >> ${GNUPLOT_FILE}
printf "set style line 88 lt 1 lw 5 lc rgb '#9B9700'\n" >> ${GNUPLOT_FILE}
printf "set style line 89 lt 1 lw 5 lc rgb '#549E79'\n" >> ${GNUPLOT_FILE}
printf "set style line 90 lt 1 lw 5 lc rgb '#FFF69F'\n" >> ${GNUPLOT_FILE}
printf "set style line 91 lt 1 lw 5 lc rgb '#201625'\n" >> ${GNUPLOT_FILE}
printf "set style line 92 lt 1 lw 5 lc rgb '#72418F'\n" >> ${GNUPLOT_FILE}
printf "set style line 93 lt 1 lw 5 lc rgb '#BC23FF'\n" >> ${GNUPLOT_FILE}
printf "set style line 94 lt 1 lw 5 lc rgb '#99ADC0'\n" >> ${GNUPLOT_FILE}
printf "set style line 95 lt 1 lw 5 lc rgb '#3A2465'\n" >> ${GNUPLOT_FILE}
printf "set style line 96 lt 1 lw 5 lc rgb '#922329'\n" >> ${GNUPLOT_FILE}
printf "set style line 97 lt 1 lw 5 lc rgb '#5B4534'\n" >> ${GNUPLOT_FILE}
printf "set style line 98 lt 1 lw 5 lc rgb '#FDE8DC'\n" >> ${GNUPLOT_FILE}
printf "set style line 99 lt 1 lw 5 lc rgb '#404E55'\n" >> ${GNUPLOT_FILE}
printf "set style line 100 lt 1 lw 5 lc rgb '#0089A3'\n" >> ${GNUPLOT_FILE}
printf "set style line 101 lt 1 lw 5 lc rgb '#CB7E98'\n" >> ${GNUPLOT_FILE}
printf "set style line 102 lt 1 lw 5 lc rgb '#A4E804'\n" >> ${GNUPLOT_FILE}
printf "set style line 103 lt 1 lw 5 lc rgb '#324E72'\n" >> ${GNUPLOT_FILE}
printf "set style line 104 lt 1 lw 5 lc rgb '#6A3A4C'\n" >> ${GNUPLOT_FILE}
printf "set style line 105 lt 1 lw 5 lc rgb '#83AB58'\n" >> ${GNUPLOT_FILE}
printf "set style line 106 lt 1 lw 5 lc rgb '#001C1E'\n" >> ${GNUPLOT_FILE}
printf "set style line 107 lt 1 lw 5 lc rgb '#D1F7CE'\n" >> ${GNUPLOT_FILE}
printf "set style line 108 lt 1 lw 5 lc rgb '#004B28'\n" >> ${GNUPLOT_FILE}
printf "set style line 109 lt 1 lw 5 lc rgb '#C8D0F6'\n" >> ${GNUPLOT_FILE}
printf "set style line 110 lt 1 lw 5 lc rgb '#A3A489'\n" >> ${GNUPLOT_FILE}
printf "set style line 111 lt 1 lw 5 lc rgb '#806C66'\n" >> ${GNUPLOT_FILE}
printf "set style line 112 lt 1 lw 5 lc rgb '#222800'\n" >> ${GNUPLOT_FILE}
printf "set style line 113 lt 1 lw 5 lc rgb '#BF5650'\n" >> ${GNUPLOT_FILE}
printf "set style line 114 lt 1 lw 5 lc rgb '#E83000'\n" >> ${GNUPLOT_FILE}
printf "set style line 115 lt 1 lw 5 lc rgb '#66796D'\n" >> ${GNUPLOT_FILE}
printf "set style line 116 lt 1 lw 5 lc rgb '#DA007C'\n" >> ${GNUPLOT_FILE}
printf "set style line 117 lt 1 lw 5 lc rgb '#FF1A59'\n" >> ${GNUPLOT_FILE}
printf "set style line 118 lt 1 lw 5 lc rgb '#8ADBB4'\n" >> ${GNUPLOT_FILE}
printf "set style line 119 lt 1 lw 5 lc rgb '#1E0200'\n" >> ${GNUPLOT_FILE}
printf "set style line 120 lt 1 lw 5 lc rgb '#5B4E51'\n" >> ${GNUPLOT_FILE}
printf "set style line 121 lt 1 lw 5 lc rgb '#C895C5'\n" >> ${GNUPLOT_FILE}
printf "set style line 122 lt 1 lw 5 lc rgb '#320033'\n" >> ${GNUPLOT_FILE}
printf "set style line 123 lt 1 lw 5 lc rgb '#FF6832'\n" >> ${GNUPLOT_FILE}
printf "set style line 124 lt 1 lw 5 lc rgb '#66E1D3'\n" >> ${GNUPLOT_FILE}
printf "set style line 125 lt 1 lw 5 lc rgb '#CFCDAC'\n" >> ${GNUPLOT_FILE}
printf "set style line 126 lt 1 lw 5 lc rgb '#D0AC94'\n" >> ${GNUPLOT_FILE}
printf "set style line 127 lt 1 lw 5 lc rgb '#7ED379'\n" >> ${GNUPLOT_FILE}
printf "set style line 128 lt 1 lw 5 lc rgb '#012C58'\n" >> ${GNUPLOT_FILE}
printf "set style line 129 lt 1 lw 5 lc rgb '#7A7BFF'\n" >> ${GNUPLOT_FILE}
printf "set style line 130 lt 1 lw 5 lc rgb '#D68E01'\n" >> ${GNUPLOT_FILE}
printf "set style line 131 lt 1 lw 5 lc rgb '#353339'\n" >> ${GNUPLOT_FILE}
printf "set style line 132 lt 1 lw 5 lc rgb '#78AFA1'\n" >> ${GNUPLOT_FILE}
printf "set style line 133 lt 1 lw 5 lc rgb '#FEB2C6'\n" >> ${GNUPLOT_FILE}
printf "set style line 134 lt 1 lw 5 lc rgb '#75797C'\n" >> ${GNUPLOT_FILE}
printf "set style line 135 lt 1 lw 5 lc rgb '#837393'\n" >> ${GNUPLOT_FILE}
printf "set style line 136 lt 1 lw 5 lc rgb '#943A4D'\n" >> ${GNUPLOT_FILE}
printf "set style line 137 lt 1 lw 5 lc rgb '#B5F4FF'\n" >> ${GNUPLOT_FILE}
printf "set style line 138 lt 1 lw 5 lc rgb '#D2DCD5'\n" >> ${GNUPLOT_FILE}
printf "set style line 139 lt 1 lw 5 lc rgb '#9556BD'\n" >> ${GNUPLOT_FILE}
printf "set style line 140 lt 1 lw 5 lc rgb '#6A714A'\n" >> ${GNUPLOT_FILE}
printf "set style line 141 lt 1 lw 5 lc rgb '#001325'\n" >> ${GNUPLOT_FILE}
printf "set style line 142 lt 1 lw 5 lc rgb '#02525F'\n" >> ${GNUPLOT_FILE}
printf "set style line 143 lt 1 lw 5 lc rgb '#0AA3F7'\n" >> ${GNUPLOT_FILE}
printf "set style line 144 lt 1 lw 5 lc rgb '#E98176'\n" >> ${GNUPLOT_FILE}
printf "set style line 145 lt 1 lw 5 lc rgb '#DBD5DD'\n" >> ${GNUPLOT_FILE}
printf "set style line 146 lt 1 lw 5 lc rgb '#5EBCD1'\n" >> ${GNUPLOT_FILE}
printf "set style line 147 lt 1 lw 5 lc rgb '#3D4F44'\n" >> ${GNUPLOT_FILE}
printf "set style line 148 lt 1 lw 5 lc rgb '#7E6405'\n" >> ${GNUPLOT_FILE}
printf "set style line 149 lt 1 lw 5 lc rgb '#02684E'\n" >> ${GNUPLOT_FILE}
printf "set style line 150 lt 1 lw 5 lc rgb '#962B75'\n" >> ${GNUPLOT_FILE}
printf "set style line 151 lt 1 lw 5 lc rgb '#8D8546'\n" >> ${GNUPLOT_FILE}
printf "set style line 152 lt 1 lw 5 lc rgb '#9695C5'\n" >> ${GNUPLOT_FILE}
printf "set style line 153 lt 1 lw 5 lc rgb '#E773CE'\n" >> ${GNUPLOT_FILE}
printf "set style line 154 lt 1 lw 5 lc rgb '#D86A78'\n" >> ${GNUPLOT_FILE}
printf "set style line 155 lt 1 lw 5 lc rgb '#3E89BE'\n" >> ${GNUPLOT_FILE}
printf "set style line 156 lt 1 lw 5 lc rgb '#CA834E'\n" >> ${GNUPLOT_FILE}
printf "set style line 157 lt 1 lw 5 lc rgb '#518A87'\n" >> ${GNUPLOT_FILE}
printf "set style line 158 lt 1 lw 5 lc rgb '#5B113C'\n" >> ${GNUPLOT_FILE}
printf "set style line 159 lt 1 lw 5 lc rgb '#55813B'\n" >> ${GNUPLOT_FILE}
printf "set style line 160 lt 1 lw 5 lc rgb '#E704C4'\n" >> ${GNUPLOT_FILE}
printf "set style line 161 lt 1 lw 5 lc rgb '#00005F'\n" >> ${GNUPLOT_FILE}
printf "set style line 162 lt 1 lw 5 lc rgb '#A97399'\n" >> ${GNUPLOT_FILE}
printf "set style line 163 lt 1 lw 5 lc rgb '#4B8160'\n" >> ${GNUPLOT_FILE}
printf "set style line 164 lt 1 lw 5 lc rgb '#59738A'\n" >> ${GNUPLOT_FILE}
printf "set style line 165 lt 1 lw 5 lc rgb '#FF5DA7'\n" >> ${GNUPLOT_FILE}
printf "set style line 166 lt 1 lw 5 lc rgb '#F7C9BF'\n" >> ${GNUPLOT_FILE}
printf "set style line 167 lt 1 lw 5 lc rgb '#643127'\n" >> ${GNUPLOT_FILE}
printf "set style line 168 lt 1 lw 5 lc rgb '#513A01'\n" >> ${GNUPLOT_FILE}
printf "set style line 169 lt 1 lw 5 lc rgb '#6B94AA'\n" >> ${GNUPLOT_FILE}
printf "set style line 170 lt 1 lw 5 lc rgb '#51A058'\n" >> ${GNUPLOT_FILE}
printf "set style line 171 lt 1 lw 5 lc rgb '#A45B02'\n" >> ${GNUPLOT_FILE}
printf "set style line 172 lt 1 lw 5 lc rgb '#1D1702'\n" >> ${GNUPLOT_FILE}
printf "set style line 173 lt 1 lw 5 lc rgb '#E20027'\n" >> ${GNUPLOT_FILE}
printf "set style line 174 lt 1 lw 5 lc rgb '#E7AB63'\n" >> ${GNUPLOT_FILE}
printf "set style line 175 lt 1 lw 5 lc rgb '#4C6001'\n" >> ${GNUPLOT_FILE}
printf "set style line 176 lt 1 lw 5 lc rgb '#9C6966'\n" >> ${GNUPLOT_FILE}
printf "set style line 177 lt 1 lw 5 lc rgb '#64547B'\n" >> ${GNUPLOT_FILE}
printf "set style line 178 lt 1 lw 5 lc rgb '#97979E'\n" >> ${GNUPLOT_FILE}
printf "set style line 179 lt 1 lw 5 lc rgb '#006A66'\n" >> ${GNUPLOT_FILE}
printf "set style line 180 lt 1 lw 5 lc rgb '#391406'\n" >> ${GNUPLOT_FILE}
printf "set style line 181 lt 1 lw 5 lc rgb '#F4D749'\n" >> ${GNUPLOT_FILE}
printf "set style line 182 lt 1 lw 5 lc rgb '#0045D2'\n" >> ${GNUPLOT_FILE}
printf "set style line 183 lt 1 lw 5 lc rgb '#006C31'\n" >> ${GNUPLOT_FILE}
printf "set style line 184 lt 1 lw 5 lc rgb '#DDB6D0'\n" >> ${GNUPLOT_FILE}
printf "set style line 185 lt 1 lw 5 lc rgb '#7C6571'\n" >> ${GNUPLOT_FILE}
printf "set style line 186 lt 1 lw 5 lc rgb '#9FB2A4'\n" >> ${GNUPLOT_FILE}
printf "set style line 187 lt 1 lw 5 lc rgb '#00D891'\n" >> ${GNUPLOT_FILE}
printf "set style line 188 lt 1 lw 5 lc rgb '#15A08A'\n" >> ${GNUPLOT_FILE}
printf "set style line 189 lt 1 lw 5 lc rgb '#BC65E9'\n" >> ${GNUPLOT_FILE}
printf "set style line 190 lt 1 lw 5 lc rgb '#FFFFFE'\n" >> ${GNUPLOT_FILE}
printf "set style line 191 lt 1 lw 5 lc rgb '#C6DC99'\n" >> ${GNUPLOT_FILE}
printf "set style line 192 lt 1 lw 5 lc rgb '#203B3C'\n" >> ${GNUPLOT_FILE}
printf "set style line 193 lt 1 lw 5 lc rgb '#671190'\n" >> ${GNUPLOT_FILE}
printf "set style line 194 lt 1 lw 5 lc rgb '#6B3A64'\n" >> ${GNUPLOT_FILE}
printf "set style line 195 lt 1 lw 5 lc rgb '#F5E1FF'\n" >> ${GNUPLOT_FILE}
printf "set style line 196 lt 1 lw 5 lc rgb '#FFA0F2'\n" >> ${GNUPLOT_FILE}
printf "set style line 197 lt 1 lw 5 lc rgb '#CCAA35'\n" >> ${GNUPLOT_FILE}
printf "set style line 198 lt 1 lw 5 lc rgb '#374527'\n" >> ${GNUPLOT_FILE}
printf "set style line 199 lt 1 lw 5 lc rgb '#8BB400'\n" >> ${GNUPLOT_FILE}
printf "set style line 200 lt 1 lw 5 lc rgb '#797868'\n" >> ${GNUPLOT_FILE}
printf "set style line 201 lt 1 lw 5 lc rgb '#C6005A'\n" >> ${GNUPLOT_FILE}
printf "set style line 202 lt 1 lw 5 lc rgb '#3B000A'\n" >> ${GNUPLOT_FILE}
printf "set style line 203 lt 1 lw 5 lc rgb '#C86240'\n" >> ${GNUPLOT_FILE}
printf "set style line 204 lt 1 lw 5 lc rgb '#29607C'\n" >> ${GNUPLOT_FILE}
printf "set style line 205 lt 1 lw 5 lc rgb '#402334'\n" >> ${GNUPLOT_FILE}
printf "set style line 206 lt 1 lw 5 lc rgb '#7D5A44'\n" >> ${GNUPLOT_FILE}
printf "set style line 207 lt 1 lw 5 lc rgb '#CCB87C'\n" >> ${GNUPLOT_FILE}
printf "set style line 208 lt 1 lw 5 lc rgb '#B88183'\n" >> ${GNUPLOT_FILE}
printf "set style line 209 lt 1 lw 5 lc rgb '#AA5199'\n" >> ${GNUPLOT_FILE}
printf "set style line 210 lt 1 lw 5 lc rgb '#B5D6C3'\n" >> ${GNUPLOT_FILE}
printf "set style line 211 lt 1 lw 5 lc rgb '#A38469'\n" >> ${GNUPLOT_FILE}
printf "set style line 212 lt 1 lw 5 lc rgb '#9F94F0'\n" >> ${GNUPLOT_FILE}
printf "set style line 213 lt 1 lw 5 lc rgb '#A74571'\n" >> ${GNUPLOT_FILE}
printf "set style line 214 lt 1 lw 5 lc rgb '#B894A6'\n" >> ${GNUPLOT_FILE}
printf "set style line 215 lt 1 lw 5 lc rgb '#71BB8C'\n" >> ${GNUPLOT_FILE}
printf "set style line 216 lt 1 lw 5 lc rgb '#00B433'\n" >> ${GNUPLOT_FILE}
printf "set style line 217 lt 1 lw 5 lc rgb '#789EC9'\n" >> ${GNUPLOT_FILE}
printf "set style line 218 lt 1 lw 5 lc rgb '#6D80BA'\n" >> ${GNUPLOT_FILE}
printf "set style line 219 lt 1 lw 5 lc rgb '#953F00'\n" >> ${GNUPLOT_FILE}
printf "set style line 220 lt 1 lw 5 lc rgb '#5EFF03'\n" >> ${GNUPLOT_FILE}
printf "set style line 221 lt 1 lw 5 lc rgb '#E4FFFC'\n" >> ${GNUPLOT_FILE}
printf "set style line 222 lt 1 lw 5 lc rgb '#1BE177'\n" >> ${GNUPLOT_FILE}
printf "set style line 223 lt 1 lw 5 lc rgb '#BCB1E5'\n" >> ${GNUPLOT_FILE}
printf "set style line 224 lt 1 lw 5 lc rgb '#76912F'\n" >> ${GNUPLOT_FILE}
printf "set style line 225 lt 1 lw 5 lc rgb '#003109'\n" >> ${GNUPLOT_FILE}
printf "set style line 226 lt 1 lw 5 lc rgb '#0060CD'\n" >> ${GNUPLOT_FILE}
printf "set style line 227 lt 1 lw 5 lc rgb '#D20096'\n" >> ${GNUPLOT_FILE}
printf "set style line 228 lt 1 lw 5 lc rgb '#895563'\n" >> ${GNUPLOT_FILE}
printf "set style line 229 lt 1 lw 5 lc rgb '#29201D'\n" >> ${GNUPLOT_FILE}
printf "set style line 230 lt 1 lw 5 lc rgb '#5B3213'\n" >> ${GNUPLOT_FILE}
printf "set style line 231 lt 1 lw 5 lc rgb '#A76F42'\n" >> ${GNUPLOT_FILE}
printf "set style line 232 lt 1 lw 5 lc rgb '#89412E'\n" >> ${GNUPLOT_FILE}
printf "set style line 233 lt 1 lw 5 lc rgb '#1A3A2A'\n" >> ${GNUPLOT_FILE}
printf "set style line 234 lt 1 lw 5 lc rgb '#494B5A'\n" >> ${GNUPLOT_FILE}
printf "set style line 235 lt 1 lw 5 lc rgb '#A88C85'\n" >> ${GNUPLOT_FILE}
printf "set style line 236 lt 1 lw 5 lc rgb '#F4ABAA'\n" >> ${GNUPLOT_FILE}
printf "set style line 237 lt 1 lw 5 lc rgb '#A3F3AB'\n" >> ${GNUPLOT_FILE}
printf "set style line 238 lt 1 lw 5 lc rgb '#00C6C8'\n" >> ${GNUPLOT_FILE}
printf "set style line 239 lt 1 lw 5 lc rgb '#EA8B66'\n" >> ${GNUPLOT_FILE}
printf "set style line 240 lt 1 lw 5 lc rgb '#958A9F'\n" >> ${GNUPLOT_FILE}
printf "set style line 241 lt 1 lw 5 lc rgb '#BDC9D2'\n" >> ${GNUPLOT_FILE}
printf "set style line 242 lt 1 lw 5 lc rgb '#9FA064'\n" >> ${GNUPLOT_FILE}
printf "set style line 243 lt 1 lw 5 lc rgb '#BE4700'\n" >> ${GNUPLOT_FILE}
printf "set style line 244 lt 1 lw 5 lc rgb '#658188'\n" >> ${GNUPLOT_FILE}
printf "set style line 245 lt 1 lw 5 lc rgb '#83A485'\n" >> ${GNUPLOT_FILE}
printf "set style line 246 lt 1 lw 5 lc rgb '#453C23'\n" >> ${GNUPLOT_FILE}
printf "set style line 247 lt 1 lw 5 lc rgb '#47675D'\n" >> ${GNUPLOT_FILE}
printf "set style line 248 lt 1 lw 5 lc rgb '#3A3F00'\n" >> ${GNUPLOT_FILE}
printf "set style line 249 lt 1 lw 5 lc rgb '#061203'\n" >> ${GNUPLOT_FILE}
printf "set style line 250 lt 1 lw 5 lc rgb '#DFFB71'\n" >> ${GNUPLOT_FILE}
printf "set style line 251 lt 1 lw 5 lc rgb '#868E7E'\n" >> ${GNUPLOT_FILE}
printf "set style line 252 lt 1 lw 5 lc rgb '#98D058'\n" >> ${GNUPLOT_FILE}
printf "set style line 253 lt 1 lw 5 lc rgb '#6C8F7D'\n" >> ${GNUPLOT_FILE}
printf "set style line 254 lt 1 lw 5 lc rgb '#D7BFC2'\n" >> ${GNUPLOT_FILE}
printf "set style line 255 lt 1 lw 5 lc rgb '#3C3E6E'\n" >> ${GNUPLOT_FILE}
printf "set style line 256 lt 1 lw 5 lc rgb '#D83D66'\n" >> ${GNUPLOT_FILE}
printf "set style line 257 lt 1 lw 5 lc rgb '#2F5D9B'\n" >> ${GNUPLOT_FILE}
printf "set style line 258 lt 1 lw 5 lc rgb '#6C5E46'\n" >> ${GNUPLOT_FILE}
printf "set style line 259 lt 1 lw 5 lc rgb '#D25B88'\n" >> ${GNUPLOT_FILE}
printf "set style line 260 lt 1 lw 5 lc rgb '#5B656C'\n" >> ${GNUPLOT_FILE}
printf "set style line 261 lt 1 lw 5 lc rgb '#00B57F'\n" >> ${GNUPLOT_FILE}
printf "set style line 262 lt 1 lw 5 lc rgb '#545C46'\n" >> ${GNUPLOT_FILE}
printf "set style line 263 lt 1 lw 5 lc rgb '#866097'\n" >> ${GNUPLOT_FILE}
printf "set style line 264 lt 1 lw 5 lc rgb '#365D25'\n" >> ${GNUPLOT_FILE}
printf "set style line 265 lt 1 lw 5 lc rgb '#252F99'\n" >> ${GNUPLOT_FILE}
printf "set style line 266 lt 1 lw 5 lc rgb '#00CCFF'\n" >> ${GNUPLOT_FILE}
printf "set style line 267 lt 1 lw 5 lc rgb '#674E60'\n" >> ${GNUPLOT_FILE}
printf "set style line 268 lt 1 lw 5 lc rgb '#FC009C'\n" >> ${GNUPLOT_FILE}
printf "set style line 269 lt 1 lw 5 lc rgb '#92896B'\n" >> ${GNUPLOT_FILE}
printf "set key outside\n" >> ${GNUPLOT_FILE}
COUNTER=0
for sim in "${SIMULATION_FOLDERS[@]}"
do
 COUNTER=$(($COUNTER+1));
 PREPLASMA_LENGTH=($(echo $sim | awk -F'_' '{print $3}'))
 RAMP_LENGTH=($(echo $sim | awk -F'_' '{print $7}'))
 BULK_LENGTH=($(echo $sim | awk -F'_' '{print $9}'))
 
 if [ -f $sim/diagnostics/diag${DIAG_TIME_TO_BE_READ}.dat.particles.txt ]
  then
  filename=`ls -1 $sim/diagnostics/diag${DIAG_TIME_TO_BE_READ}.dat.particles.txt`

  if [ "$COUNTER" -eq 1 ]
  then
   printf "plot \"$filename\" u (\$1*c):%s w lines ls %s t 'pre_%s_ramp_%s_bulk_%s',\\" "${PROTON_EMAX_COLUMN}" "${COUNTER}" "${PREPLASMA_LENGTH}" "${RAMP_LENGTH}" "${BULK_LENGTH}" >> ${GNUPLOT_FILE}
   printf "\n" >> ${GNUPLOT_FILE}
  elif [ "$COUNTER" -eq ${#ArrayName[@]} ]
  then
   printf " \"$filename\" u (\$1*c):%s w lines ls %s t 'pre_%s_ramp_%s_bulk_%s'" "${PROTON_EMAX_COLUMN}" "${COUNTER}" "${PREPLASMA_LENGTH}" "${RAMP_LENGTH}" "${BULK_LENGTH}" >> ${GNUPLOT_FILE}
   printf "\n" >> ${GNUPLOT_FILE}
  else
   printf " \"$filename\" u (\$1*c):%s w lines ls %s t 'pre_%s_ramp_%s_bulk_%s',\\" "${PROTON_EMAX_COLUMN}" "${COUNTER}" "${PREPLASMA_LENGTH}" "${RAMP_LENGTH}" "${BULK_LENGTH}" >> ${GNUPLOT_FILE}
   printf "\n" >> ${GNUPLOT_FILE}
  fi
 fi
done

gnuplot ${GNUPLOT_FILE}

