#! /bin/bash

if [ $# != 5 ]
then
 echo "Usage: please give 5 arguments to the script:"
 echo "\$1 : put 0 for new simulations, 1 for restarts"
 echo "\$2 : initial sim id (e.g. put '4' if you want the output to reload from *out04.*)"
 echo "\$3 : tmax (e.g. 100.0)"
 echo "\$4 : number of time the program will dump outputs"
 echo "\$5 : number of diagnostic outputs produced during the simulation"
 exit
fi

MATERIAL="CH2"
#MATERIAL="AL"

HPC="cnaf"
#HPC="galileo"
#HPC="fermi"
JOBFILE=job.cmd

#use true to really submit jobs
#PRODUCTION=false
PRODUCTION=true



#preplasmas=($(awk 'BEGIN{for(i=1.0;i<=3.0;i+=1.0)print i}'))
preplasmas=0.0

#foam_densities=($(awk 'BEGIN{for(i=0.5;i<=3.0;i+=0.5)print i}'))
foam_densities=1.0

#ramps=($(awk 'BEGIN{for(i=0.25;i<=0.75;i+=0.25)print i}'))
ramps=0.0

centrals=($(awk 'BEGIN{for(i=2.5;i<=7.5;i+=2.5)print i}'))
#centrals=0.5
#centrals=(2.0 4.0 8.0)

#bulk_densities=($(awk 'BEGIN{for(i=20.0;i<=120.0;i+=20.0)print i}'))
bulk_densities=100.0

#contams=($(awk 'BEGIN{for(i=0.05;i<=0.1;i+=0.01)print i}'))
contams=0.1

angles=($(awk 'BEGIN{for(i=0.0;i<=30.0;i+=15.0)print i}'))
#angles=0.0

laser_anotes=($(awk 'BEGIN{for(i=10.0;i<=30.0;i+=10.0)print i}'))
#laser_anotes=17.0

laser_lengths=($(awk 'BEGIN{for(i=25.0;i<=30.0;i+=5.0)print i}'))
#laser_lengths=40.0







###############################
## do not touch from here :) ##
###############################


for pre in ${preplasmas[*]}
do
for f_dens in ${foam_densities[*]}
do
for ramp in ${ramps[*]}
do
for central in ${centrals[*]}
do
for c_dens in ${bulk_densities[*]}
do
for contam in ${contams[*]}
do
for angle in ${angles[*]}
do
for anote in ${laser_anotes[*]}
do
for l_length in ${laser_lengths[*]}
do

echo "${angle}_${pre}_${f_dens}_${ramp}_${central}_${c_dens}_${contam}_${anote}_${l_length}" >> sim_da_fare.txt
INPUTFILE=input.nml

if [ "$1" = "0" ]
	then
mkdir "${angle}_${pre}_${f_dens}_${ramp}_${central}_${c_dens}_${contam}_${anote}_${l_length}"
fi

cd "${angle}_${pre}_${f_dens}_${ramp}_${central}_${c_dens}_${contam}_${anote}_${l_length}" || return

if [ "$1" = "0" ] ; then
cp ../${JOBFILE} .
fi

NCPU=32
CREA_FILE_DUMP=0


##### nx, ny,nz
## nb:  il numero punti di griglia in y e z, se diversi da 1, deve essere divisibile per il numero di processori
##      se la simulazione e' di PWFA, anche il numero di punti di griglia in x deve essere divisibile per il num
##      di processori
## nb2: se si abilita lo stretching, il numero di punti in y e z deve anche essere divisibile per 6, dato che la
##      zona di stretching e' definita come 1/6 dei punti totali
## nb3: per ora e' meglio mettere in z un numero di punti uguale a quello in y, a meno che non si metta 1 per le simulazioni 2D

NUMERO_PUNTI_GRIGLIA_X=14336
NUMERO_PUNTI_GRIGLIA_Y=14336

#NUMERO_PUNTI_GRIGLIA_X=7168
#NUMERO_PUNTI_GRIGLIA_Y=7168
##NB: mettere il seguente numero ad 1 per fare simulazioni 2D
NUMERO_PUNTI_GRIGLIA_Z=1

#NUMERO_PUNTI_GRIGLIA_TRASVERSALI_OCCUPATI_DAL_PLASMA=1500
NUMERO_PUNTI_GRIGLIA_TRASVERSALI_OCCUPATI_DAL_PLASMA=14000


##nb: ricordarsi di risolvere la skin depth
NUMERO_PUNTI_GRIGLIA_PER_MICRON=100.0
RAPPORTO_DIMENSIONE_GRIGLIA_TRASVERSA_GRIGLIA_LONGITUDINALE=1.0

##### a questo punto sappiamo le dimensioni in micron della griglia 
##### Dimx = nx/k0    Dimy = ny*yx_rat/k0
##### Si potrebbero anche scrivere a video per l'utente

##### LPf,Der,str,iform
ORDINE_INTEGRAZIONE_LEAPFROG=2
TIPO_DERIVATA=2
##NB: per il seguente dato, 0=griglia standard, 1=stretching trasversale, 2=stretching anche lungo x
STRETCHING=0
## 0 open, 1 periodiche (da verificare)
BOUNDARY_X=0
BOUNDARY_Y=0
BOUNDARY_Z=0

##NB: per il seguente dato, 0=algoritmo standard, 1=conservazione carica, 2=conservazione energia
ALGORITMO_INTEGRAZIONE=0


##### mdl,dmdl,ibeam 
##NB: per il seguente dato, 1=polarizzato p, 2=polarizzato s, 3=polarizzato circolarmente
LASER_MODEL=1

##NB: per il seguente dato, ecco le chiamate effettuate nel codice
if [ "$MATERIAL" = "CH2" ] ; then
PLASMA_MODEL=1
else
PLASMA_MODEL=3
fi



#### nsp,nsb,Z_i,A_i 
##NB: se il plasma model e' 1, nsp deve essere 2 oppure 1: il significato e' che se e' 1 i protoni sono fissi (genero solo gli elettroni), 
##                                                         se e' 2 invece genero entrambe le specie
##    se il plasma model e' 2, nsp deve essere uguale a 4
##    se il plasma model e' 4, riferirsi alla guida in PLASMA_MODEL
NUMERO_SPECIE=3
##NB: per la documentazione riguardante il seguente numero, riferirsi alla spiegazione di PLASMA_MODEL
NUMERO_SPECIE_SECONDARIO=1
##NB: per i seguenti valori di solito si usa ionizzazione 9 e peso 27 (alluminio)
if [ "$MATERIAL" = "CH2" ] ; then
NUMERO_IONIZZAZIONE_Z1=4
NUMERO_ATOMICO_SPECIE_1=6
NUMERO_MASSA_SPECIE_1="12.0"
NUMERO_IONIZZAZIONE_Z1_MAX=4
else
NUMERO_IONIZZAZIONE_Z1=9
NUMERO_ATOMICO_SPECIE_1=13
NUMERO_MASSA_SPECIE_1="27.0"
NUMERO_IONIZZAZIONE_Z1_MAX=9
fi
NUMERO_IONIZZAZIONE_Z2=1
NUMERO_ATOMICO_SPECIE_2=1
NUMERO_MASSA_SPECIE_2="1.0"
NUMERO_IONIZZAZIONE_Z2_MAX=1

IONIZZAZIONE_ATTIVA=0
IONIZZAZIONE_MODELLO=1

TEMP_INIZIALE_PLASMA_1=0.0003
TEMP_INIZIALE_PLASMA_2=0.0
TEMP_INIZIALE_PLASMA_3=0.0
TEMP_INIZIALE_PLASMA_4=0.0

if [ "$MATERIAL" = "CH2" ] ; then
#### np_xc
NUMERO_ELETTRONI_LONGITUDINALMENTE_PER_CELLA=12
NUMERO_IONI_LONGITUDINALMENTE_PER_CELLA=4
NUMERO_PROTONI_LONGITUDINALMENTE_PER_CELLA=6
#### np_yc
NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA=12
NUMERO_IONI_TRASVERSALMENTE_PER_CELLA=6
NUMERO_PROTONI_TRASVERSALMENTE_PER_CELLA=8
else

#### np_xc
NUMERO_ELETTRONI_LONGITUDINALMENTE_PER_CELLA_LAYER_CENTRALE=12
NUMERO_IONI_LONGITUDINALMENTE_PER_CELLA_LAYER_CENTRALE=4
NUMERO_ELETTRONI_LONGITUDINALMENTE_PER_CELLA_LAYER_FRONTALE=4
NUMERO_IONI_LONGITUDINALMENTE_PER_CELLA_LAYER_FRONTALE=4
NUMERO_ELETTRONI_LONGITUDINALMENTE_PER_CELLA_LAYER_POSTERIORE=4
NUMERO_IONI_LONGITUDINALMENTE_PER_CELLA_LAYER_POSTERIORE=4
#### np_yc
NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA_LAYER_CENTRALE=9
NUMERO_IONI_TRASVERSALMENTE_PER_CELLA_LAYER_CENTRALE=3
NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA_LAYER_FRONTALE=4
NUMERO_IONI_TRASVERSALMENTE_PER_CELLA_LAYER_FRONTALE=4
NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA_LAYER_POSTERIORE=4
NUMERO_IONI_TRASVERSALMENTE_PER_CELLA_LAYER_POSTERIORE=4
fi

#### a questo punto e' meglio controllare che np_xc_el_l2*np_yc_el_l2 sia
#### uguale a Z_i*np_xc_ion_l2*np_yc_ion_l2
#### NON DOVREBBE ESSER NECESSARIO, MA PASQUALE INSISTE SU QUESTO PUNTO:
#### MEGLIO CHE TUTTE LE SPECIE ABBIANO LO STESSO PESO


#### t0, xc, wx, wy, a0,lam0 ## tutti in micrometri
POSIZIONE_INIZIALE_PICCO_IMPULSO_LASER=$(bc -l <<< "scale=2;(${l_length}/2)")
DISTANZA_INIZIALE_PICCO_IMPULSO_LASER_DAL_FUOCO=$(bc -l <<< "scale=2;(${l_length}/2)")
LUNGHEZZA_LASER_TOTALE=${l_length}
WAIST_LASER=3.0
PARAMETRO_ADIMENSIONALE_LASER_A0=${anote}
LUNGHEZZA_ONDA_LASER=0.8



#### lx(1:7)
SPESSORE_LAYER_FRONTALE=${pre}
SPESSORE_RAMPA_LAYER_FRONTALE_LAYER_CENTRALE=${ramp}
SPESSORE_LAYER_CENTRALE=${central}
SPESSORE_RAMPA_LAYER_CENTRALE_LAYER_POSTERIORE=0.0
SPESSORE_LAYER_POSTERIORE=${contam}
ANGOLO_ROTAZIONE_LASER=${angle}

if [ "$ANGOLO_ROTAZIONE_LASER" = "0.0" ] ; then
OFFSET_FINE_LASER_INIZIO_TARGHETTA=0.01
else
OFFSET_FINE_LASER_INIZIO_TARGHETTA=$(bc -l <<< "scale=5;(0.01+${WAIST_LASER})*s(${ANGOLO_ROTAZIONE_LASER}/57.29578)/c(${ANGOLO_ROTAZIONE_LASER}/57.29578)")
fi

WIRE_SIZE=0.0
INTERWIRE_DISTANCE=0.0
#### a questo punto sappiamo la posizione della targhetta
#### inizio_targhetta = xc + wx/2 + lx7


#### n/nc, n1/nc,n2/nc (tutte le densita' sono quindi rapportate alla densita' critica)
DENSITA_ELETTRONI_LAYER_FRONTALE=${f_dens}
DENSITA_ELETTRONI_LAYER_CENTRALE=${c_dens}
DENSITA_ELETTRONI_LAYER_POSTERIORE=10.0


#### nf, nd, npv, end_p
####
NUMERO_OUTPUT_CAMPI=1
NUMERO_OUTPUT_DENSITA_GRIGLIA=3
NUMERO_OUTPUT_SPAZIOFASI_PARTICELLE=0
NUMERO_OUTPUT_SPAZIOFASI_BUNCH=0

#### X1-X0 seguenti definiscono quindi l'altezza del parallelepipedo (o rettangolo in 2D) di spazio da dumpare
X0_TAGLIO_OUTPUT=0.0
X1_TAGLIO_OUTPUT=100.0
#### Il seguente valore invece definisce il semilato di base del parallelepipedo (o del rettangolo in 2D) dello spazio da dumpare
SEMILATO_BASE_TAGLIO_OUTPUT=20.0


#### jmp, pjmp
#### NB: Attenzione: è importante che JUMP_GRIGLIA sia un divisore del numero di punti di 
####     griglia per processore, altrimenti l’output delle griglie è corrotto.
JUMP_GRIGLIA=1
JUMP_PARTICELLE=1


#### wnd_sh, w_in, w_end, w_speed
#### numero time steps ogni quanti viene invocate la routine di moving window
MW_CALL_EVERY_N_TIMESTEPS=20
#### time step inizio movimento della window
MW_START_TIME=120.0
#### time step fine movimento della window
MW_END_TIME=120.0
#### velocità beta con cui si muove la window
#MW_SPEED=0.6
MW_SPEED=1.0



####  per il seguente dato, 1=laser driven simulation
####                        2=bunch driven simulation
SYM_TYPE=1


#### il seguente dato deve sempre essere <= 1.0
COURANT_FRIEDRICHS_LEWY_PARAMETER=0.8




########################################
########################################
## FINE PARAMETRI - NON TOCCARE OLTRE ##
########################################
########################################
nx=${NUMERO_PUNTI_GRIGLIA_X}
ny=${NUMERO_PUNTI_GRIGLIA_Y}
nz=${NUMERO_PUNTI_GRIGLIA_Z}
ny_targ=${NUMERO_PUNTI_GRIGLIA_TRASVERSALI_OCCUPATI_DAL_PLASMA}

k0=${NUMERO_PUNTI_GRIGLIA_PER_MICRON}
yx_rat=${RAPPORTO_DIMENSIONE_GRIGLIA_TRASVERSA_GRIGLIA_LONGITUDINALE}

lpf_ord=${ORDINE_INTEGRAZIONE_LEAPFROG}
der_ord=${TIPO_DERIVATA}
str_flag=${STRETCHING}
iform=${ALGORITMO_INTEGRAZIONE}

model_id=${LASER_MODEL}
dmodel_id=${PLASMA_MODEL}
ibx=${BOUNDARY_X}
iby=${BOUNDARY_Y}
ibz=${BOUNDARY_Z}
ibeam=${SYM_TYPE}

nsp=${NUMERO_SPECIE}
nsb=${NUMERO_SPECIE_SECONDARIO}
Z1_ion=${NUMERO_IONIZZAZIONE_Z1}
A1_ion=${NUMERO_ATOMICO_SPECIE_1}
M1_ion=${NUMERO_MASSA_SPECIE_1}
Z1_max=${NUMERO_IONIZZAZIONE_Z1_MAX}
Z2_ion=${NUMERO_IONIZZAZIONE_Z2}
A2_ion=${NUMERO_ATOMICO_SPECIE_2}
M2_ion=${NUMERO_MASSA_SPECIE_2}
Z2_max=${NUMERO_IONIZZAZIONE_Z2_MAX}
ionz_lev=${IONIZZAZIONE_ATTIVA}
ionz_model=${IONIZZAZIONE_MODELLO}

t0_pl_1=${TEMP_INIZIALE_PLASMA_1}
t0_pl_2=${TEMP_INIZIALE_PLASMA_2}
t0_pl_3=${TEMP_INIZIALE_PLASMA_3}
t0_pl_4=${TEMP_INIZIALE_PLASMA_4}

if [ "$MATERIAL" = "CH2" ] ; then
np_per_xc_1=${NUMERO_ELETTRONI_LONGITUDINALMENTE_PER_CELLA}
np_per_xc_2=${NUMERO_IONI_LONGITUDINALMENTE_PER_CELLA}
np_per_xc_3=${NUMERO_PROTONI_LONGITUDINALMENTE_PER_CELLA}
np_per_xc_4=1
np_per_xc_5=1
np_per_xc_6=1
np_per_yc_1=${NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA}
np_per_yc_2=${NUMERO_IONI_TRASVERSALMENTE_PER_CELLA}
np_per_yc_3=${NUMERO_PROTONI_TRASVERSALMENTE_PER_CELLA}
np_per_yc_4=1
np_per_yc_5=1
np_per_yc_6=1
else
np_per_xc_1=${NUMERO_ELETTRONI_LONGITUDINALMENTE_PER_CELLA_LAYER_CENTRALE}
np_per_xc_2=${NUMERO_IONI_LONGITUDINALMENTE_PER_CELLA_LAYER_CENTRALE}
np_per_xc_3=${NUMERO_ELETTRONI_LONGITUDINALMENTE_PER_CELLA_LAYER_FRONTALE}
np_per_xc_4=${NUMERO_IONI_LONGITUDINALMENTE_PER_CELLA_LAYER_FRONTALE}
np_per_xc_5=${NUMERO_ELETTRONI_LONGITUDINALMENTE_PER_CELLA_LAYER_POSTERIORE}
np_per_xc_6=${NUMERO_IONI_LONGITUDINALMENTE_PER_CELLA_LAYER_POSTERIORE}
np_per_yc_1=${NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA_LAYER_CENTRALE}
np_per_yc_2=${NUMERO_IONI_TRASVERSALMENTE_PER_CELLA_LAYER_CENTRALE}
np_per_yc_3=${NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA_LAYER_FRONTALE}
np_per_yc_4=${NUMERO_IONI_TRASVERSALMENTE_PER_CELLA_LAYER_FRONTALE}
np_per_yc_5=${NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA_LAYER_POSTERIORE}
np_per_yc_6=${NUMERO_IONI_TRASVERSALMENTE_PER_CELLA_LAYER_POSTERIORE}
fi

t0_lp=${DISTANZA_INIZIALE_PICCO_IMPULSO_LASER_DAL_FUOCO}
xc_lp=${POSIZIONE_INIZIALE_PICCO_IMPULSO_LASER}
w0_x=${LUNGHEZZA_LASER_TOTALE}
w0_y=${WAIST_LASER}
a0=${PARAMETRO_ADIMENSIONALE_LASER_A0}
lam0=${LUNGHEZZA_ONDA_LASER}


lpx_1=$(echo "${SPESSORE_LAYER_FRONTALE}*1.0" | bc -l)
lpx_2=$(echo "${SPESSORE_RAMPA_LAYER_FRONTALE_LAYER_CENTRALE}*1.0" | bc -l)
lpx_3=$(echo "${SPESSORE_LAYER_CENTRALE}*1.0" | bc -l)
lpx_4=$(echo "${SPESSORE_RAMPA_LAYER_CENTRALE_LAYER_POSTERIORE}*1.0" | bc -l)
lpx_5=$(echo "${SPESSORE_LAYER_POSTERIORE}*1.0" | bc -l)
lpx_6=$(echo "${ANGOLO_ROTAZIONE_LASER}*1.0" | bc -l)
lpx_7=$(echo "${OFFSET_FINE_LASER_INIZIO_TARGHETTA}*1.0" | bc -l)

lpy_1=${WIRE_SIZE}
lpy_2=${INTERWIRE_DISTANCE}

n_over_nc=${DENSITA_ELETTRONI_LAYER_CENTRALE}
n1_over_n=${DENSITA_ELETTRONI_LAYER_FRONTALE}
n2_over_n=${DENSITA_ELETTRONI_LAYER_POSTERIORE}

w_sh=${MW_CALL_EVERY_N_TIMESTEPS}
wi_time=${MW_START_TIME}
wf_time=${MW_END_TIME}
w_speed=${MW_SPEED}

nouts=$4
iene=$5
nvout=${NUMERO_OUTPUT_CAMPI}
nden=${NUMERO_OUTPUT_DENSITA_GRIGLIA}
npout=${NUMERO_OUTPUT_SPAZIOFASI_PARTICELLE}
nbout=${NUMERO_OUTPUT_SPAZIOFASI_BUNCH}

jump=${JUMP_GRIGLIA}
pjump=${JUMP_PARTICELLE}

xp0_out=${X0_TAGLIO_OUTPUT}
xp1_out=${X1_TAGLIO_OUTPUT}
yp_out=${SEMILATO_BASE_TAGLIO_OUTPUT}

tmax=$3
cfl=${COURANT_FRIEDRICHS_LEWY_PARAMETER}

new_sim=$1
id_new=$2
dump=${CREA_FILE_DUMP}
npe_yz=${NCPU}


 if [ ! -f ./dumpRestart/dumpout000001.bin ] && [ "$1" = "1" ] ; then
   echo "Missing dump files! Aborting submitting ${angle}_${pre}_${f_dens}_${ramp}_${central}_${c_dens}_${contam}_${anote}_${l_length} !"
   cd ..
   continue
 fi

 if [ "$1" = "1" ] ; then
  mv $INPUTFILE "$INPUTFILE.old$PREVIOUS_STEP"
 fi

 rm -f ${INPUTFILE}
 touch ${INPUTFILE}

{
 printf '&GRID\n'
 printf ' nx = %s,\n' "$nx"
 printf ' ny = %s,\n' "$ny"
 printf ' nz = %s,\n' "$nz"
 printf ' ny_targ = %s,\n' "${ny_targ}"
 printf ' k0 = %s,\n' "$k0"
 printf ' yx_rat = %s,\n' "${yx_rat}"
 printf ' zx_rat = %s\n' "${yx_rat}"
 printf '/'
 printf '\n\n'

 printf '&SIMULATION\n'
 printf ' lpf_ord = %s,\n' "${lpf_ord}"
 printf ' der_ord = %s,\n' "${der_ord}"
 printf ' str_flag = %s,\n' "${str_flag}"
 printf ' iform = %s,\n' "$iform"
 printf ' model_id = %s,\n' "${model_id}"
 printf ' dmodel_id = %s,\n' "${dmodel_id}"
 printf ' ibx = %s,\n' "$ibx"
 printf ' iby = %s,\n' "$iby"
 printf ' ibz = %s,\n' "$ibz"
 printf ' ibeam = %s\n' "$ibeam"
 printf '/'
 printf '\n\n'

 printf '&TARGET_DESCRIPTION\n'
 printf ' nsp = %s,\n' "$nsp"
 printf ' nsb = %s,\n' "$nsb"
 printf ' ionz_lev = %s,\n' "${ionz_lev}"
 printf ' ionz_model = %s,\n' "${ionz_model}"
 printf ' ion_min(1) = %s,\n' "${Z1_ion}"
 printf ' ion_min(2) = %s,\n' "${Z2_ion}"
 printf ' ion_min(3) = 1,\n'
 printf ' ion_max(1) = %s,\n' "${Z1_max}"
 printf ' ion_max(2) = %s,\n' "${Z2_max}"
 printf ' ion_max(3) = 1,\n'
 printf ' atomic_number(1) = %s,\n' "${A1_ion}"
 printf ' atomic_number(2) = %s,\n' "${A2_ion}"
 printf ' atomic_number(3) = 1,\n'
 printf ' mass_number(1) = %s,\n' "${M1_ion}"
 printf ' mass_number(2) = %s,\n' "${M2_ion}"
 printf ' mass_number(3) = 1.0,\n'
 printf ' t0_pl(1) = %s,\n' "${t0_pl_1}"
 printf ' t0_pl(2) = %s,\n' "${t0_pl_2}"
 printf ' t0_pl(3) = %s,\n' "${t0_pl_3}"
 printf ' t0_pl(4) = %s,\n' "${t0_pl_4}"
 printf ' np_per_xc(1) = %s,\n' "${np_per_xc_1}"
 printf ' np_per_xc(2) = %s,\n' "${np_per_xc_2}"
 printf ' np_per_xc(3) = %s,\n' "${np_per_xc_3}"
 printf ' np_per_xc(4) = %s,\n' "${np_per_xc_4}"
 printf ' np_per_xc(5) = %s,\n' "${np_per_xc_5}"
 printf ' np_per_xc(6) = %s,\n' "${np_per_xc_6}"
 printf ' np_per_yc(1) = %s,\n' "${np_per_yc_1}"
 printf ' np_per_yc(2) = %s,\n' "${np_per_yc_2}"
 printf ' np_per_yc(3) = %s,\n' "${np_per_yc_3}"
 printf ' np_per_yc(4) = %s,\n' "${np_per_yc_4}"
 printf ' np_per_yc(5) = %s,\n' "${np_per_yc_5}"
 printf ' np_per_yc(6) = %s,\n' "${np_per_yc_6}"
 printf ' lpx(1) = %s,\n' "$lpx_1"
 printf ' lpx(2) = %s,\n' "$lpx_2"
 printf ' lpx(3) = %s,\n' "$lpx_3"
 printf ' lpx(4) = %s,\n' "$lpx_4"
 printf ' lpx(5) = %s,\n' "$lpx_5"
 printf ' lpx(6) = %s,\n' "$lpx_6"
 printf ' lpx(7) = %s,\n' "$lpx_7"
 printf ' lpy(1) = %s,\n' "$lpy_1"
 printf ' lpy(2) = %s,\n' "$lpy_2"
 printf ' n_over_nc = %s,\n' "${n_over_nc}"
 printf ' n1_over_n = %s,\n' "${n1_over_n}"
 printf ' n2_over_n = %s\n' "${n2_over_n}"
 printf '/'
 printf '\n\n'

 printf '&LASER\n'
 printf ' t0_lp = %s,\n' "$t0_lp"
 printf ' xc_lp = %s,\n' "$xc_lp"
 printf ' w0_x = %s,\n' "$w0_x"
 printf ' w0_y = %s,\n' "$w0_y"
 printf ' a0 = %s,\n' "$a0"
 printf ' lam0 = %s\n' "$lam0"
 printf '/'
 printf '\n\n'

 printf '&MOVING_WINDOW\n'
 printf ' w_sh = %s,\n' "${w_sh}"
 printf ' wi_time = %s,\n' "${wi_time}"
 printf ' wf_time = %s,\n' "${wf_time}"
 printf ' w_speed = %s\n' "${w_speed}"
 printf '/'
 printf '\n\n'

 printf '&OUTPUT\n'
 printf ' nouts = %s,\n' "$nouts"
 printf ' iene = %s,\n' "$iene"
 printf ' nvout = %s,\n' "$nvout"
 printf ' nden = %s,\n' "$nden"
 printf ' npout = %s,\n' "$npout"
 printf ' nbout = %s,\n' "$nbout"
 printf ' jump = %s,\n' "$jump"
 printf ' pjump = %s,\n' "$pjump"
 printf ' xp0_out = %s,\n' "${xp0_out}"
 printf ' xp1_out = %s,\n' "${xp1_out}"
 printf ' yp_out = %s,\n' "${yp_out}"
 printf ' tmax = %s,\n' "$tmax"
 printf ' cfl = %s,\n' "$cfl"
 printf ' new_sim = %s,\n' "${new_sim}"
 printf ' id_new = %s,\n' "${id_new}"
 printf ' dump = %s\n' "$dump"
 printf '/'
 printf '\n\n'

 printf '&MPIPARAMS\n'
 printf ' nprocx = 1,\n'
 printf ' nprocy = %s,\n' "${npe_yz}"
 printf ' nprocz = 1\n'
 printf '/'
 printf '\n\n'
} >> ${INPUTFILE}

 if [ "$PRODUCTION" = "true" ] ; then
  if [ "$HPC" = "cnaf" ] ; then
   bsub < $JOBFILE
  elif [ "$HPC" = "galileo" ] ; then
   qsub $JOBFILE
  elif [ "$HPC" = "fermi" ] ; then
   llsubmit $JOBFILE
  fi
 fi

cd ..

done
done
done
done
done
done
done
done
done


