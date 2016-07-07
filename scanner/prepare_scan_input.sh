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

#MATERIAL="CH2"
#MATERIAL="CD2"
MATERIAL="AL"

HPC="cnaf"
#HPC="galileo"
#HPC="fermi"
JOBFILE=job.cmd

#use true to really submit jobs
#PRODUCTION=false
PRODUCTION=true



#preplasmas=$(awk 'BEGIN{for(i=1.0;i<=3.0;i+=1.0)print i}')
preplasmas=0.0

#densities=$(awk 'BEGIN{for(i=0.5;i<=3.0;i+=0.5)print i}')
densities=1.0

#ramps=$(awk 'BEGIN{for(i=0.25;i<=0.75;i+=0.25)print i}')
ramps=0.0

#centrals=$(awk 'BEGIN{for(i=2.0;i<=4.0;i+=1.0)print i}')
centrals=2.0

#contams=$(awk 'BEGIN{for(i=0.05;i<=0.1;i+=0.01)print i}')
contams=0.08

angles=$(awk 'BEGIN{for(i=5.0;i<=15.0;i+=5.0)print i}')
#angles=0.0








###############################
## do not touch from here :) ##
###############################


for pre in $preplasmas
do
for dens in $densities
do
for ramp in $ramps
do
for central in $centrals
do
for contam in $contams
do
for angle in $angles
do

echo "${angle}_pre_${pre}_den_${dens}_ramp_${ramp}_cent_${central}_cont_${contam}" >> sim_da_fare.txt
INPUTFILE=input.nml

if [ "$1" = "0" ] ; then
mkdir ${angle}_pre_${pre}_den_${dens}_ramp_${ramp}_cent_${central}_cont_${contam}
fi

cd ${angle}_pre_${pre}_den_${dens}_ramp_${ramp}_cent_${central}_cont_${contam}

if [ "$1" = "0" ] ; then
cp ../${JOBFILE} .
fi

NCPU=64
CREA_FILE_DUMP=1


##### nx, ny,nz
## nb:  il numero punti di griglia in y e z, se diversi da 1, deve essere divisibile per il numero di processori
##      se la simulazione e' di PWFA, anche il numero di punti di griglia in x deve essere divisibile per il num
##      di processori
## nb2: se si abilita lo stretching, il numero di punti in y e z deve anche essere divisibile per 6, dato che la
##      zona di stretching e' definita come 1/6 dei punti totali
## nb3: per ora e' meglio mettere in z un numero di punti uguale a quello in y, a meno che non si metta 1 per le simulazioni 2D
#NUMERO_PUNTI_GRIGLIA_X=3072
#NUMERO_PUNTI_GRIGLIA_Y=1536
NUMERO_PUNTI_GRIGLIA_X=7168
NUMERO_PUNTI_GRIGLIA_Y=7168
##NB: mettere il seguente numero ad 1 per fare simulazioni 2D
NUMERO_PUNTI_GRIGLIA_Z=1

#NUMERO_PUNTI_GRIGLIA_TRASVERSALI_OCCUPATI_DAL_PLASMA=1500
NUMERO_PUNTI_GRIGLIA_TRASVERSALI_OCCUPATI_DAL_PLASMA=7000


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

## DA DOCUMENTARE
IBEAM=1

##NB: per il seguente dato, 0=algoritmo standard, 1=conservazione carica, 2=conservazione energia
ALGORITMO_INTEGRAZIONE=0


##### mdl,dmdl,ibeam 
##NB: per il seguente dato, 1=polarizzato p, 2=polarizzato s, 3=polarizzato circolarmente
LASER_MODEL=1


##NB: per il seguente dato, ecco le chiamate effettuate nel codice
# case(1)
#  call one_layer_multisp(ny_targ,xf0)  !e+Z1+Z2
# case(3)
#  call preplasma_multisp(ny_targ,xf0) !foam+ high Z +H coating
# case(4)
#  call multi_layer_multisp(ny_targ,xf0) !foam+ high Z +H coating
# case(5)
#  call one_layer_nano_wires(ny_targ,xf0) !foam+ high Z +H coating mass lim
# case(6)
#  call one_layer_nano_tubes(ny_targ,xf0)  ! (e+ (Z,A) ions now nsp=2
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
POSIZIONE_INIZIALE_PICCO_IMPULSO_LASER=16.5
DISTANZA_INIZIALE_PICCO_IMPULSO_LASER_DAL_FUOCO=16.5
LUNGHEZZA_LASER_TOTALE=33.0
WAIST_LASER=6.2
PARAMETRO_ADIMENSIONALE_LASER_A0=3.0
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
DENSITA_ELETTRONI_LAYER_FRONTALE=${dens}
DENSITA_ELETTRONI_LAYER_CENTRALE=100.0
DENSITA_ELETTRONI_LAYER_POSTERIORE=10.0


#### nf, nd, npv, end_p
####
NUMERO_OUTPUT_CAMPI=1
NUMERO_OUTPUT_DENSITA_GRIGLIA=4
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


lpx_1=`echo "${SPESSORE_LAYER_FRONTALE}*1.0" | bc -l`
lpx_2=`echo "${SPESSORE_RAMPA_LAYER_FRONTALE_LAYER_CENTRALE}*1.0" | bc -l`
lpx_3=`echo "${SPESSORE_LAYER_CENTRALE}*1.0" | bc -l`
lpx_4=`echo "${SPESSORE_RAMPA_LAYER_CENTRALE_LAYER_POSTERIORE}*1.0" | bc -l`
lpx_5=`echo "${SPESSORE_LAYER_POSTERIORE}*1.0" | bc -l`
lpx_6=`echo "${ANGOLO_ROTAZIONE_LASER}*1.0" | bc -l`
lpx_7=`echo "${OFFSET_FINE_LASER_INIZIO_TARGHETTA}*1.0" | bc -l`

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
   echo "Missing dump files! Aborting submitting pre_${pre}_den_${dens}_ramp_${ramp}_cent_${central}_cont_${contam} !"
   cd ..
   continue
 fi

 if [ "$1" = "1" ] ; then
  mv $INPUTFILE $INPUTFILE.old$PREVIOUS_STEP
 fi

 rm -f ${INPUTFILE}
 touch ${INPUTFILE}


 printf '&GRID\n' >> ${INPUTFILE}
 printf ' nx = %s,\n' "$nx" >> ${INPUTFILE}
 printf ' ny = %s,\n' "$ny" >> ${INPUTFILE}
 printf ' nz = %s,\n' "$nz" >> ${INPUTFILE}
 printf ' ny_targ = %s,\n' "${ny_targ}" >> ${INPUTFILE}
 printf ' k0 = %s,\n' "$k0" >> ${INPUTFILE}
 printf ' yx_rat = %s,\n' "${yx_rat}" >> ${INPUTFILE}
 printf ' zx_rat = %s\n' "${yx_rat}" >> ${INPUTFILE}
 printf '/' >> ${INPUTFILE}
 printf '\n\n' >> ${INPUTFILE}

 printf '&SIMULATION\n' >> ${INPUTFILE}
 printf ' lpf_ord = %s,\n' "${lpf_ord}" >> ${INPUTFILE}
 printf ' der_ord = %s,\n' "${der_ord}" >> ${INPUTFILE}
 printf ' str_flag = %s,\n' "${str_flag}" >> ${INPUTFILE}
 printf ' iform = %s,\n' "$iform" >> ${INPUTFILE}
 printf ' model_id = %s,\n' "${model_id}" >> ${INPUTFILE}
 printf ' dmodel_id = %s,\n' "${dmodel_id}" >> ${INPUTFILE}
 printf ' ibx = %s,\n' "$ibx" >> ${INPUTFILE}
 printf ' iby = %s,\n' "$iby" >> ${INPUTFILE}
 printf ' ibz = %s,\n' "$ibz" >> ${INPUTFILE}
 printf ' ibeam = %s\n' "$ibeam" >> ${INPUTFILE}
 printf '/' >> ${INPUTFILE}
 printf '\n\n' >> ${INPUTFILE}

 printf '&TARGET_DESCRIPTION\n' >> ${INPUTFILE}
 printf ' nsp = %s,\n' "$nsp" >> ${INPUTFILE}
 printf ' nsb = %s,\n' "$nsb" >> ${INPUTFILE}
 printf ' ionz_lev = %s,\n' "${ionz_lev}" >> ${INPUTFILE}
 printf ' ionz_model = %s,\n' "${ionz_model}" >> ${INPUTFILE}
 printf ' ion_min(1) = %s,\n' "${Z1_ion}" >> ${INPUTFILE}
 printf ' ion_min(2) = %s,\n' "${Z2_ion}" >> ${INPUTFILE}
 printf ' ion_min(3) = 1,\n' >> ${INPUTFILE}
 printf ' ion_max(1) = %s,\n' "${Z1_max}" >> ${INPUTFILE}
 printf ' ion_max(2) = %s,\n' "${Z2_max}" >> ${INPUTFILE}
 printf ' ion_max(3) = 1,\n' >> ${INPUTFILE}
 printf ' atomic_number(1) = %s,\n' "${A1_ion}" >> ${INPUTFILE}
 printf ' atomic_number(2) = %s,\n' "${A2_ion}" >> ${INPUTFILE}
 printf ' atomic_number(3) = 1,\n' >> ${INPUTFILE}
 printf ' mass_number(1) = %s,\n' "${M1_ion}" >> ${INPUTFILE}
 printf ' mass_number(2) = %s,\n' "${M2_ion}" >> ${INPUTFILE}
 printf ' mass_number(3) = 1.0,\n' >> ${INPUTFILE}
 printf ' t0_pl(1) = %s,\n' "${t0_pl_1}" >> ${INPUTFILE}
 printf ' t0_pl(2) = %s,\n' "${t0_pl_2}" >> ${INPUTFILE}
 printf ' t0_pl(3) = %s,\n' "${t0_pl_3}" >> ${INPUTFILE}
 printf ' t0_pl(4) = %s,\n' "${t0_pl_4}" >> ${INPUTFILE}
 printf ' np_per_xc(1) = %s,\n' "${np_per_xc_1}" >> ${INPUTFILE}
 printf ' np_per_xc(2) = %s,\n' "${np_per_xc_2}" >> ${INPUTFILE}
 printf ' np_per_xc(3) = %s,\n' "${np_per_xc_3}" >> ${INPUTFILE}
 printf ' np_per_xc(4) = %s,\n' "${np_per_xc_4}" >> ${INPUTFILE}
 printf ' np_per_xc(5) = %s,\n' "${np_per_xc_5}" >> ${INPUTFILE}
 printf ' np_per_xc(6) = %s,\n' "${np_per_xc_6}" >> ${INPUTFILE}
 printf ' np_per_yc(1) = %s,\n' "${np_per_yc_1}" >> ${INPUTFILE}
 printf ' np_per_yc(2) = %s,\n' "${np_per_yc_2}" >> ${INPUTFILE}
 printf ' np_per_yc(3) = %s,\n' "${np_per_yc_3}" >> ${INPUTFILE}
 printf ' np_per_yc(4) = %s,\n' "${np_per_yc_4}" >> ${INPUTFILE}
 printf ' np_per_yc(5) = %s,\n' "${np_per_yc_5}" >> ${INPUTFILE}
 printf ' np_per_yc(6) = %s,\n' "${np_per_yc_6}" >> ${INPUTFILE}
 printf ' lpx(1) = %s,\n' "$lpx_1" >> ${INPUTFILE}
 printf ' lpx(2) = %s,\n' "$lpx_2" >> ${INPUTFILE}
 printf ' lpx(3) = %s,\n' "$lpx_3" >> ${INPUTFILE}
 printf ' lpx(4) = %s,\n' "$lpx_4" >> ${INPUTFILE}
 printf ' lpx(5) = %s,\n' "$lpx_5" >> ${INPUTFILE}
 printf ' lpx(6) = %s,\n' "$lpx_6" >> ${INPUTFILE}
 printf ' lpx(7) = %s,\n' "$lpx_7" >> ${INPUTFILE}
 printf ' lpy(1) = %s,\n' "$lpy_1" >> ${INPUTFILE}
 printf ' lpy(2) = %s,\n' "$lpy_2" >> ${INPUTFILE}
 printf ' n_over_nc = %s,\n' "${n_over_nc}" >> ${INPUTFILE}
 printf ' n1_over_n = %s,\n' "${n1_over_n}" >> ${INPUTFILE}
 printf ' n2_over_n = %s\n' "${n2_over_n}" >> ${INPUTFILE}
 printf '/' >> ${INPUTFILE}
 printf '\n\n' >> ${INPUTFILE}

 printf '&LASER\n' >> ${INPUTFILE}
 printf ' t0_lp = %s,\n' "$t0_lp" >> ${INPUTFILE}
 printf ' xc_lp = %s,\n' "$xc_lp" >> ${INPUTFILE}
 printf ' w0_x = %s,\n' "$w0_x" >> ${INPUTFILE}
 printf ' w0_y = %s,\n' "$w0_y" >> ${INPUTFILE}
 printf ' a0 = %s,\n' "$a0" >> ${INPUTFILE}
 printf ' lam0 = %s\n' "$lam0" >> ${INPUTFILE}
 printf '/' >> ${INPUTFILE}
 printf '\n\n' >> ${INPUTFILE}

 printf '&MOVING_WINDOW\n' >> ${INPUTFILE}
 printf ' w_sh = %s,\n' "${w_sh}" >> ${INPUTFILE}
 printf ' wi_time = %s,\n' "${wi_time}" >> ${INPUTFILE}
 printf ' wf_time = %s,\n' "${wf_time}" >> ${INPUTFILE}
 printf ' w_speed = %s\n' "${w_speed}" >> ${INPUTFILE}
 printf '/' >> ${INPUTFILE}
 printf '\n\n' >> ${INPUTFILE}

 printf '&OUTPUT\n' >> ${INPUTFILE}
 printf ' nouts = %s,\n' "$nouts" >> ${INPUTFILE}
 printf ' iene = %s,\n' "$iene" >> ${INPUTFILE}
 printf ' nvout = %s,\n' "$nvout" >> ${INPUTFILE}
 printf ' nden = %s,\n' "$nden" >> ${INPUTFILE}
 printf ' npout = %s,\n' "$npout" >> ${INPUTFILE}
 printf ' nbout = %s,\n' "$nbout" >> ${INPUTFILE}
 printf ' jump = %s,\n' "$jump" >> ${INPUTFILE}
 printf ' pjump = %s,\n' "$pjump" >> ${INPUTFILE}
 printf ' xp0_out = %s,\n' "${xp0_out}" >> ${INPUTFILE}
 printf ' xp1_out = %s,\n' "${xp1_out}" >> ${INPUTFILE}
 printf ' yp_out = %s,\n' "${yp_out}" >> ${INPUTFILE}
 printf ' tmax = %s,\n' "$tmax" >> ${INPUTFILE}
 printf ' cfl = %s,\n' "$cfl" >> ${INPUTFILE}
 printf ' new_sim = %s,\n' "${new_sim}" >> ${INPUTFILE}
 printf ' id_new = %s,\n' "${id_new}" >> ${INPUTFILE}
 printf ' dump = %s\n' "$dump" >> ${INPUTFILE}
 printf '/' >> ${INPUTFILE}
 printf '\n\n' >> ${INPUTFILE}

 printf '&MPIPARAMS\n' >> ${INPUTFILE}
 printf ' nprocx = 1,\n' >> ${INPUTFILE}
 printf ' nprocy = %s,\n' "${npe_yz}" >> ${INPUTFILE}
 printf ' nprocz = 1\n' >> ${INPUTFILE}
 printf '/' >> ${INPUTFILE}
 printf '\n\n' >> ${INPUTFILE}

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




