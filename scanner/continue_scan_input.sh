#! /bin/bash

if [ $# != 3 ]
then
 echo "In input devono esser passati:"
 echo "\$1 : t_finale (tipo '10.0')"
 echo "\$2 : numero output completi (tipo '2')"
 echo "\$3 : numero output ridotti (tipo '20')"
 exit
fi

JOBFILE=galileo-64.cmd
PREVIOUS_STEP=5
INPUTFILE=input.nml

preplasmas=$(awk 'BEGIN{for(i=1.0;i<=3.0;i+=1.0)print i}')
#preplasmas=0.0

#densities=$(awk 'BEGIN{for(i=0.5;i<=3.0;i+=0.5)print i}')
densities=1.0

ramps=$(awk 'BEGIN{for(i=0.25;i<=0.75;i+=0.25)print i}')
#ramps=0.0

#centrals=$(awk 'BEGIN{for(i=2.0;i<=4.0;i+=1.0)print i}')
centrals=$(awk 'BEGIN{for(i=5.0;i<=10.0;i+=1.0)print i}')
#centrals=2.4

#contams=$(awk 'BEGIN{for(i=0.05;i<=0.1;i+=0.01)print i}')
contams=0.08



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

echo "pre_${pre}_den_${dens}_ramp_${ramp}_cent_${central}_cont_${contam}" >> sim_da_fare.txt

cd pre_${pre}_den_${dens}_ramp_${ramp}_cent_${central}_cont_${contam}

ALADYN_VERSION=3
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
NUMERO_PUNTI_GRIGLIA_Y=3584
##NB: mettere il seguente numero ad 1 per fare simulazioni 2D
NUMERO_PUNTI_GRIGLIA_Z=1

#NUMERO_PUNTI_GRIGLIA_TRASVERSALI_OCCUPATI_DAL_PLASMA=1500
NUMERO_PUNTI_GRIGLIA_TRASVERSALI_OCCUPATI_DAL_PLASMA=3500


##nb: ricordarsi di risolvere la skin depth
NUMERO_PUNTI_GRIGLIA_PER_MICRON=100.0
RAPPORTO_DIMENSIONE_GRIGLIA_TRASVERSA_GRIGLIA_LONGITUDINALE=2.0

##### a questo punto sappiamo le dimensioni in micron della griglia 
##### Dimx = nx/k0    Dimy = ny*yx_rat/k0
##### Si potrebbero anche scrivere a video per l'utente

##### LPf,Der,str,iform
ORDINE_INTEGRAZIONE_LEAPFROG=2
TIPO_DERIVATA=2
##NB: per il seguente dato, 0=griglia standard, 1=stretching trasversale, 2=stretching anche lungo x, 3=PML
TIPO_BOUNDARIES=2
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
PLASMA_MODEL=3




#### nsp,nsb,Z_i,A_i 
##NB: se il plasma model e' 1, nsp deve essere 2 oppure 1: il significato e' che se e' 1 i protoni sono fissi (genero solo gli elettroni), 
##                                                         se e' 2 invece genero entrambe le specie
##    se il plasma model e' 2, nsp deve essere uguale a 4
##    se il plasma model e' 4, riferirsi alla guida in PLASMA_MODEL
NUMERO_SPECIE=3
##NB: per la documentazione riguardante il seguente numero, riferirsi alla spiegazione di PLASMA_MODEL
NUMERO_SPECIE_SECONDARIO=1
##NB: per i seguenti valori di solito si usa ionizzazione 9 e peso 27 (alluminio)
NUMERO_IONIZZAZIONE_Z1=10
NUMERO_ATOMICO_SPECIE_1=13
NUMERO_IONIZZAZIONE_Z1_MAX=10
NUMERO_IONIZZAZIONE_Z2=1
NUMERO_ATOMICO_SPECIE_2=1
NUMERO_IONIZZAZIONE_Z2_MAX=1

IONIZZAZIONE_ATTIVA=0

TEMP_INIZIALE_PLASMA_1=0.0003
TEMP_INIZIALE_PLASMA_2=0.0
TEMP_INIZIALE_PLASMA_3=0.0
TEMP_INIZIALE_PLASMA_4=0.0

#### np_xc
#NUMERO_ELETTRONI_LONGITUDINALMENTE_PER_CELLA_LAYER_CENTRALE=4
#NUMERO_IONI_LONGITUDINALMENTE_PER_CELLA_LAYER_CENTRALE=1
#NUMERO_ELETTRONI_LONGITUDINALMENTE_PER_CELLA_LAYER_FRONTALE=2
#NUMERO_IONI_LONGITUDINALMENTE_PER_CELLA_LAYER_FRONTALE=2
#NUMERO_ELETTRONI_LONGITUDINALMENTE_PER_CELLA_LAYER_POSTERIORE=2
#NUMERO_IONI_LONGITUDINALMENTE_PER_CELLA_LAYER_POSTERIORE=2
NUMERO_ELETTRONI_LONGITUDINALMENTE_PER_CELLA_LAYER_CENTRALE=12
NUMERO_IONI_LONGITUDINALMENTE_PER_CELLA_LAYER_CENTRALE=4
NUMERO_ELETTRONI_LONGITUDINALMENTE_PER_CELLA_LAYER_FRONTALE=4
NUMERO_IONI_LONGITUDINALMENTE_PER_CELLA_LAYER_FRONTALE=4
NUMERO_ELETTRONI_LONGITUDINALMENTE_PER_CELLA_LAYER_POSTERIORE=4
NUMERO_IONI_LONGITUDINALMENTE_PER_CELLA_LAYER_POSTERIORE=4


#### np_yc
#NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA_LAYER_CENTRALE=3
#NUMERO_IONI_TRASVERSALMENTE_PER_CELLA_LAYER_CENTRALE=2
#NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA_LAYER_FRONTALE=2
#NUMERO_IONI_TRASVERSALMENTE_PER_CELLA_LAYER_FRONTALE=2
#NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA_LAYER_POSTERIORE=2
#NUMERO_IONI_TRASVERSALMENTE_PER_CELLA_LAYER_POSTERIORE=2
NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA_LAYER_CENTRALE=9
NUMERO_IONI_TRASVERSALMENTE_PER_CELLA_LAYER_CENTRALE=3
NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA_LAYER_FRONTALE=4
NUMERO_IONI_TRASVERSALMENTE_PER_CELLA_LAYER_FRONTALE=4
NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA_LAYER_POSTERIORE=4
NUMERO_IONI_TRASVERSALMENTE_PER_CELLA_LAYER_POSTERIORE=4


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
ANGOLO_ROTAZIONE_LASER=0.0
OFFSET_FINE_LASER_INIZIO_TARGHETTA=0.01


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
NUMERO_OUTPUT_CAMPI=0
NUMERO_OUTPUT_DENSITA_GRIGLIA=0
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
str_flag=${TIPO_BOUNDARIES}
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
Z1_max=${NUMERO_IONIZZAZIONE_Z1_MAX}
Z2_ion=${NUMERO_IONIZZAZIONE_Z2}
A2_ion=${NUMERO_ATOMICO_SPECIE_2}
Z2_max=${NUMERO_IONIZZAZIONE_Z2_MAX}
zmod=${IONIZZAZIONE_ATTIVA}

t0_pl_1=${TEMP_INIZIALE_PLASMA_1}
t0_pl_2=${TEMP_INIZIALE_PLASMA_2}
t0_pl_3=${TEMP_INIZIALE_PLASMA_3}
t0_pl_4=${TEMP_INIZIALE_PLASMA_4}

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

t0_lp=${POSIZIONE_INIZIALE_PICCO_IMPULSO_LASER}
xc_lp=${DISTANZA_INIZIALE_PICCO_IMPULSO_LASER_DAL_FUOCO}
w0_x=${LUNGHEZZA_LASER_TOTALE}
w0_y=${WAIST_LASER}
a0=${PARAMETRO_ADIMENSIONALE_LASER_A0}
lam0=${LUNGHEZZA_ONDA_LASER}


lpx_1=${SPESSORE_LAYER_FRONTALE}
lpx_2=${SPESSORE_RAMPA_LAYER_FRONTALE_LAYER_CENTRALE}
lpx_3=${SPESSORE_LAYER_CENTRALE}
lpx_4=${SPESSORE_RAMPA_LAYER_CENTRALE_LAYER_POSTERIORE}
lpx_5=${SPESSORE_LAYER_POSTERIORE}
lpx_6=${ANGOLO_ROTAZIONE_LASER}
lpx_7=${OFFSET_FINE_LASER_INIZIO_TARGHETTA}

lpy_1=${WIRE_SIZE}
lpy_2=${INTERWIRE_DISTANCE}

n_over_nc=${DENSITA_ELETTRONI_LAYER_CENTRALE}
n1_over_n=${DENSITA_ELETTRONI_LAYER_FRONTALE}
n2_over_n=${DENSITA_ELETTRONI_LAYER_POSTERIORE}

w_sh=${MW_CALL_EVERY_N_TIMESTEPS}
wi_time=${MW_START_TIME}
wf_time=${MW_END_TIME}
w_speed=${MW_SPEED}

nouts=$2
iene=$3
nvout=${NUMERO_OUTPUT_CAMPI}
nden=${NUMERO_OUTPUT_DENSITA_GRIGLIA}
npout=${NUMERO_OUTPUT_SPAZIOFASI_PARTICELLE}
nbout=${NUMERO_OUTPUT_SPAZIOFASI_BUNCH}

jump=${JUMP_GRIGLIA}
pjump=${JUMP_PARTICELLE}

xp0_out=${X0_TAGLIO_OUTPUT}
xp1_out=${X1_TAGLIO_OUTPUT}
yp_out=${SEMILATO_BASE_TAGLIO_OUTPUT}

tmax=$1
cfl=${COURANT_FRIEDRICHS_LEWY_PARAMETER}

new_sim=1
id_new=${PREVIOUS_STEP}
dump=${CREA_FILE_DUMP}
npe_yz=${NCPU}


 if [ ! -f ./dumpout000001.bin ]
 then
   echo "Missing dump files! Aborting submitting pre_${pre}_den_${dens}_ramp_${ramp}_cent_${central}_cont_${contam} !"
   cd ..
   continue
 fi

 mv $INPUTFILE $INPUTFILE.old$PREVIOUS_STEP
 touch ${INPUTFILE}


 printf '&VERSION\n' >> ${INPUTFILE}
 printf ' aladyn_version = %s,\n' "${ALADYN_VERSION}">> ${INPUTFILE}
 printf '/' >> ${INPUTFILE}
 printf '\n\n' >> ${INPUTFILE}

 printf '&GRID\n' >> ${INPUTFILE}
 printf ' nx = %s,\n' "$nx" >> ${INPUTFILE}
 printf ' ny = %s,\n' "$ny" >> ${INPUTFILE}
 printf ' nz = %s,\n' "$nz" >> ${INPUTFILE}
 printf ' ny_targ = %s,\n' "${ny_targ}" >> ${INPUTFILE}
 printf ' k0 = %s,\n' "$k0" >> ${INPUTFILE}
 printf ' yx_rat = %s,\n' "${yx_rat}" >> ${INPUTFILE}
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
 printf ' ibeam = %s,\n' "$ibeam" >> ${INPUTFILE}
 printf ' nsp = %s,\n' "$nsp" >> ${INPUTFILE}
 printf ' nsb = %s,\n' "$nsb" >> ${INPUTFILE}
 printf ' Z1_ion = %s,\n' "${Z1_ion}" >> ${INPUTFILE}
 printf ' A1_ion = %s,\n' "${A1_ion}" >> ${INPUTFILE}
 printf ' Z1_max = %s,\n' "${Z1_max}" >> ${INPUTFILE}
 printf ' Z2_ion = %s,\n' "${Z2_ion}" >> ${INPUTFILE}
 printf ' A2_ion = %s,\n' "${A2_ion}" >> ${INPUTFILE}
 printf ' Z2_max = %s,\n' "${Z2_max}" >> ${INPUTFILE}
 printf ' zmod = %s,\n' "$zmod" >> ${INPUTFILE}
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
 printf ' t0_lp = %s,\n' "$t0_lp" >> ${INPUTFILE}
 printf ' xc_lp = %s,\n' "$xc_lp" >> ${INPUTFILE}
 printf ' w0_x = %s,\n' "$w0_x" >> ${INPUTFILE}
 printf ' w0_y = %s,\n' "$w0_y" >> ${INPUTFILE}
 printf ' a0 = %s,\n' "$a0" >> ${INPUTFILE}
 printf ' lam0 = %s,\n' "$lam0" >> ${INPUTFILE}
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
 printf ' n2_over_n = %s,\n' "${n2_over_n}" >> ${INPUTFILE}
 printf ' w_sh = %s,\n' "${w_sh}" >> ${INPUTFILE}
 printf ' wi_time = %s,\n' "${wi_time}" >> ${INPUTFILE}
 printf ' wf_time = %s,\n' "${wf_time}" >> ${INPUTFILE}
 printf ' w_speed = %s,\n' "${w_speed}" >> ${INPUTFILE}
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
 printf ' dump = %s,\n' "$dump" >> ${INPUTFILE}
 printf ' npe_yz = %s,\n' "${npe_yz}" >> ${INPUTFILE}
 printf '/' >> ${INPUTFILE}
 printf '\n\n' >> ${INPUTFILE}

qsub $JOBFILE

cd ..

done
done
done
done
done

