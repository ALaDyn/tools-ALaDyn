#! /bin/bash

if [ $# != 5 ]
then
 echo "In input devono esser passati:"
 echo "\$1 : new o restart ('0' oppure '1')"
 echo "\$2 : id file iniziale (tipo '4' per partire con *out04.*, ad esempio)"
 echo "\$3 : t_finale (tipo '10.0')"
 echo "\$4 : numero output completi (tipo '2')"
 echo "\$5 : numero output ridotti (tipo '20')"
 exit
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
NUMERO_PUNTI_GRIGLIA_X=7168
NUMERO_PUNTI_GRIGLIA_Y=3584
##NB: mettere il seguente numero ad 1 per fare simulazioni 2D
NUMERO_PUNTI_GRIGLIA_Z=1

NUMERO_PUNTI_GRIGLIA_TRASVERSALI_OCCUPATI_DAL_PLASMA=3500


##nb: ricordarsi di risolvere la skin depth
NUMERO_PUNTI_GRIGLIA_PER_MICRON="100."
RAPPORTO_DIMENSIONE_GRIGLIA_TRASVERSA_GRIGLIA_LONGITUDINALE="2.0"

##### a questo punto sappiamo le dimensioni in micron della griglia 
##### Dimx = nx/k0    Dimy = ny*yx_rat/k0
##### Si potrebbero anche scrivere a video per l'utente

##### LPf,Der,str,iform
ORDINE_INTEGRAZIONE_LEAPFROG=2
TIPO_DERIVATA=2
##NB: per il seguente dato, 0=griglia standard, 1=stretching trasversale, 2=stretching anche lungo x, 3=PML
TIPO_BOUNDARIES=2
##NB: per il seguente dato, 1=conservazione carica, 2=conservazione energia
ALGORITMO_INTEGRAZIONE=1


##### mdl,dmdl,ibeam 
##NB: per il seguente dato, 1=polarizzato p, 2=polarizzato s, 3=polarizzato circolarmente
LASER_MODEL=1
##NB: per il seguente dato, 1=plasma di soli elettroni e protoni, in tutti e 5 i layers
##                          2=plasma a 4 specie, con la targhetta centrale (piu' rampa eventuale) fatta di Heavy Ions ed elettroni (ascolta solo x2,x3), 
##                            permette di definire Z_i ed A_i, piu' layer posteriore di contaminanti fatto di Light Ions+H+elettroni
##                          3=plasma di 3 specie: C, H ed elettroni, con il rapporto CHn definito dal BEAM_MODEL seguente
##                          4=plasma a 3 o 4 specie, con foam anteriore fatta di elettroni e H (con BEAM_MODEL=0 e NUMERO_SPECIE=3), 
##                            oppure di carbonio "light ions" (con BEAM_MODEL=0 e NUMERO_SPECIE=4) oppure infine dello stesso materiale
##                            del bulk della targhetta "heavy ions" (con BEAM_MODEL=1 e NUMERO_SPECIE=3);
##                            la targhetta centrale e' sempre fatta di "Heavy Ions" ed elettroni, e i contaminanti di solo idrogeno+elettroni
PLASMA_MODEL=4


#### nsp,nsb,Z_i,A_i 
##NB: se il plasma model e' 1, nsp deve essere 2 oppure 1: il significato e' che se e' 1 i protoni sono fissi (genero solo gli elettroni), 
##                                                         se e' 2 invece genero entrambe le specie
##    se il plasma model e' 2, nsp deve essere uguale a 4
##    se il plasma model e' 4, riferirsi alla guida in PLASMA_MODEL
NUMERO_SPECIE=3
##NB: per la documentazione riguardante il seguente numero, riferirsi alla spiegazione di PLASMA_MODEL
NUMERO_SPECIE_SECONDARIO=0
##NB: per i seguenti valori di solito si usa ionizzazione 9 e peso 27 (alluminio)
NUMERO_IONIZZAZIONE_Z=9
NUMERO_DI_NUCLEONI_A=27


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


#### a questo punto e' meglio controllare che np_xc_el_l2*np_yc_el_l2 sia
#### uguale a Z_i*np_xc_ion_l2*np_yc_ion_l2
#### NON DOVREBBE ESSER NECESSARIO, MA PASQUALE INSISTE SU QUESTO PUNTO:
#### MEGLIO CHE TUTTE LE SPECIE ABBIANO LO STESSO PESO


#### t0, xc, wx, wy, a0,lam0 ## tutti in micrometri
POSIZIONE_INIZIALE_PICCO_IMPULSO_LASER="12.1"
DISTANZA_INIZIALE_PICCO_IMPULSO_LASER_DAL_FUOCO="12.0"
LUNGHEZZA_LASER_TOTALE="24.0"
WAIST_LASER="12.0"
PARAMETRO_ADIMENSIONALE_LASER_A0="5.0"
LUNGHEZZA_ONDA_LASER="0.8"



#### lx(1:7)
SPESSORE_LAYER_FRONTALE="0.08"
SPESSORE_RAMPA_LAYER_FRONTALE_LAYER_CENTRALE="0.0"
SPESSORE_LAYER_CENTRALE="0.6"
SPESSORE_RAMPA_LAYER_CENTRALE_LAYER_POSTERIORE="0.0"
SPESSORE_LAYER_POSTERIORE="0.08"
ANGOLO_ROTAZIONE_LASER="0.0"
OFFSET_FINE_LASER_INIZIO_TARGHETTA="0.01"


#### a questo punto sappiamo la posizione della targhetta
#### inizio_targhetta = xc + wx/2 + lx7


#### n/nc, n1/nc,n2/nc (tutte le densita' sono quindi rapportate alla densita' critica)
DENSITA_ELETTRONI_LAYER_FRONTALE="10.0"
DENSITA_ELETTRONI_LAYER_CENTRALE="100."
DENSITA_ELETTRONI_LAYER_POSTERIORE="10.0"


#### nf, nd, npv, end_p
####
NUMERO_OUTPUT_CAMPI=3
NUMERO_OUTPUT_DENSITA_GRIGLIA=1
NUMERO_OUTPUT_SPAZIOFASI_PARTICELLE=4
#### il flag seguente, se vale 1, impone un taglio nello spazio xyz secondo le dimensioni descritte dai
#### parametri imposti nelle tre seguenti righe, al fine di ridurre le dimensioni dei file.
FLAG_TAGLIO_OUTPUT_SPAZIOFASI_PARTICELLE=0
#### X1-X0 seguenti definiscono quindi l'altezza del parallelepipedo (o rettangolo in 2D) di spazio da dumpare
X0_TAGLIO_OUTPUT="0.0"
X1_TAGLIO_OUTPUT="100.0"
#### Il seguente valore invece definisce il semilato di base del parallelepipedo (o del rettangolo in 2D) dello spazio da dumpare
SEMILATO_BASE_TAGLIO_OUTPUT="20.0"


#### jmp, pjmp
#### NB: Attenzione: è importante che JUMP_GRIGLIA sia un divisore del numero di punti di 
####     griglia per processore, altrimenti l’output delle griglie è corrotto.
JUMP_GRIGLIA=2
JUMP_PARTICELLE=2


#### wnd_sh, w_in, w_end, w_speed
#### numero time steps ogni quanti viene invocate la routine di moving window
MW_CALL_EVERY_N_TIMESTEPS=20
#### time step inizio movimento della window
MW_START_TIME="120."
#### time step fine movimento della window
MW_END_TIME="120."
#### velocità beta con cui si muove la window
MW_SPEED="1.0"



####  per il seguente dato, 1=laser driven simulation
####                        2=bunch driven simulation
SYM_TYPE=1


#### il seguente dato deve sempre essere <= 1.0
COURANT_FRIEDRICHS_LEWY_PARAMETER="0.85"




########################################
########################################
## FINE PARAMETRI - NON TOCCARE OLTRE ##
########################################
########################################
nx=${NUMERO_PUNTI_GRIGLIA_X}
ny=${NUMERO_PUNTI_GRIGLIA_Y}
nz=${NUMERO_PUNTI_GRIGLIA_Z}
nplasma=${NUMERO_PUNTI_GRIGLIA_TRASVERSALI_OCCUPATI_DAL_PLASMA}

k0=${NUMERO_PUNTI_GRIGLIA_PER_MICRON}
yx_rat=${RAPPORTO_DIMENSIONE_GRIGLIA_TRASVERSA_GRIGLIA_LONGITUDINALE}

LPf=${ORDINE_INTEGRAZIONE_LEAPFROG}
Der=${TIPO_DERIVATA}
str=${TIPO_BOUNDARIES}
iform=${ALGORITMO_INTEGRAZIONE}

mdl=${LASER_MODEL}
dmdl=${PLASMA_MODEL}
ibeam=${SYM_TYPE}

nsp=${NUMERO_SPECIE}
nsb=${NUMERO_SPECIE_SECONDARIO}
Z_i=${NUMERO_IONIZZAZIONE_Z}
A_i=${NUMERO_DI_NUCLEONI_A}

np_xc_el_l2=${NUMERO_ELETTRONI_LONGITUDINALMENTE_PER_CELLA_LAYER_CENTRALE}
np_xc_ion_l2=${NUMERO_IONI_LONGITUDINALMENTE_PER_CELLA_LAYER_CENTRALE}
np_xc_el_l1=${NUMERO_ELETTRONI_LONGITUDINALMENTE_PER_CELLA_LAYER_FRONTALE}
np_xc_ion_l1=${NUMERO_IONI_LONGITUDINALMENTE_PER_CELLA_LAYER_FRONTALE}
np_xc_el_l3=${NUMERO_ELETTRONI_LONGITUDINALMENTE_PER_CELLA_LAYER_POSTERIORE}
np_xc_ion_l3=${NUMERO_IONI_LONGITUDINALMENTE_PER_CELLA_LAYER_POSTERIORE}

np_yc_el_l2=${NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA_LAYER_CENTRALE}
np_yc_ion_l2=${NUMERO_IONI_TRASVERSALMENTE_PER_CELLA_LAYER_CENTRALE}
np_yc_el_l1=${NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA_LAYER_FRONTALE}
np_yc_ion_l1=${NUMERO_IONI_TRASVERSALMENTE_PER_CELLA_LAYER_FRONTALE}
np_yc_el_l3=${NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA_LAYER_POSTERIORE}
np_yc_ion_l3=${NUMERO_IONI_TRASVERSALMENTE_PER_CELLA_LAYER_POSTERIORE}

t0=${POSIZIONE_INIZIALE_PICCO_IMPULSO_LASER}
xc=${DISTANZA_INIZIALE_PICCO_IMPULSO_LASER_DAL_FUOCO}
wx=${LUNGHEZZA_LASER_TOTALE}
wy=${WAIST_LASER}
a0=${PARAMETRO_ADIMENSIONALE_LASER_A0}
lam0=${LUNGHEZZA_ONDA_LASER}


lx1=${SPESSORE_LAYER_FRONTALE}
lx2=${SPESSORE_RAMPA_LAYER_FRONTALE_LAYER_CENTRALE}
lx3=${SPESSORE_LAYER_CENTRALE}
lx4=${SPESSORE_RAMPA_LAYER_CENTRALE_LAYER_POSTERIORE}
lx5=${SPESSORE_LAYER_POSTERIORE}
lx6=${ANGOLO_ROTAZIONE_LASER}
lx7=${OFFSET_FINE_LASER_INIZIO_TARGHETTA}

n_nc=${DENSITA_ELETTRONI_LAYER_CENTRALE}
n1_nc=${DENSITA_ELETTRONI_LAYER_FRONTALE}
n2_nc=${DENSITA_ELETTRONI_LAYER_POSTERIORE}

wnd_sh=${MW_CALL_EVERY_N_TIMESTEPS}
w_in=${MW_START_TIME}
w_end=${MW_END_TIME}
w_speed=${MW_SPEED}

nout=$4
iene=$5
nf=${NUMERO_OUTPUT_CAMPI}
nd=${NUMERO_OUTPUT_DENSITA_GRIGLIA}
npv=${NUMERO_OUTPUT_SPAZIOFASI_PARTICELLE}
end_p=${FLAG_TAGLIO_OUTPUT_SPAZIOFASI_PARTICELLE}

jmp=${JUMP_GRIGLIA}
pjmp=${JUMP_PARTICELLE}

xp0=${X0_TAGLIO_OUTPUT}
xp1=${X1_TAGLIO_OUTPUT}
ypmax=${SEMILATO_BASE_TAGLIO_OUTPUT}

tmax=$3
cfl=${COURANT_FRIEDRICHS_LEWY_PARAMETER}

new=$1
id_ew=$2
dump=${CREA_FILE_DUMP}
pey=${NCPU}



 INPUTFILE=input.nml

 rm  ${INPUTFILE}
 touch ${INPUTFILE}

 printf '&GRID\n' >> ${INPUTFILE}
 printf ' nx = %s,\n' "$nx" >> ${INPUTFILE}
 printf ' ny = %s,\n' "$ny" >> ${INPUTFILE}
 printf ' nz = %s,\n' "$nz" >> ${INPUTFILE}
 printf ' ny_targ = %s,\n' "$nplasma" >> ${INPUTFILE}
 printf ' k0 = %s,\n' "$k0" >> ${INPUTFILE}
 printf ' yx_rat = %s,\n' "$yx_rat" >> ${INPUTFILE}
 printf '/' >> ${INPUTFILE}
 printf '\n\n' >> ${INPUTFILE}

 printf '&SIMULATION\n' >> ${INPUTFILE}
 printf ' LPf_ord = %s,\n' "$LPf" >> ${INPUTFILE}
 printf ' Der_ord = %s,\n' "$Der" >> ${INPUTFILE}
 printf ' Str_flag = %s,\n' "$str" >> ${INPUTFILE}
 printf ' iform = %s,\n' "$iform" >> ${INPUTFILE}
 printf ' model_id = %s,\n' "$mdl" >> ${INPUTFILE}
 printf ' dmodel_id = %s,\n' "$dmdl" >> ${INPUTFILE}
 printf ' ibeam = %s,\n' "$ibeam" >> ${INPUTFILE}
 printf ' nsp = %s,\n' "$nsp" >> ${INPUTFILE}
 printf ' nsb = %s,\n' "$nsb" >> ${INPUTFILE}
 printf ' Z_ion = %s,\n' "${Z_i}" >> ${INPUTFILE}
 printf ' A_ion = %s,\n' "${A_i}" >> ${INPUTFILE}
 printf ' np_per_xc(1) = %s,\n' "${np_xc_el_l2}" >> ${INPUTFILE}
 printf ' np_per_xc(2) = %s,\n' "${np_xc_ion_l2}" >> ${INPUTFILE}
 printf ' np_per_xc(3) = %s,\n' "${np_xc_el_l1}" >> ${INPUTFILE}
 printf ' np_per_xc(4) = %s,\n' "${np_xc_ion_l1}" >> ${INPUTFILE}
 printf ' np_per_xc(5) = %s,\n' "${np_xc_el_l3}" >> ${INPUTFILE}
 printf ' np_per_xc(6) = %s,\n' "${np_xc_ion_l3}" >> ${INPUTFILE}
 printf ' np_per_yc(1) = %s,\n' "${np_yc_el_l2}" >> ${INPUTFILE}
 printf ' np_per_yc(2) = %s,\n' "${np_yc_ion_l2}" >> ${INPUTFILE}
 printf ' np_per_yc(3) = %s,\n' "${np_yc_el_l1}" >> ${INPUTFILE}
 printf ' np_per_yc(4) = %s,\n' "${np_yc_ion_l1}" >> ${INPUTFILE}
 printf ' np_per_yc(5) = %s,\n' "${np_yc_el_l3}" >> ${INPUTFILE}
 printf ' np_per_yc(6) = %s,\n' "${np_yc_ion_l3}" >> ${INPUTFILE}
 printf ' t0_lp = %s,\n' "$t0" >> ${INPUTFILE}
 printf ' xc_lp = %s,\n' "$xc" >> ${INPUTFILE}
 printf ' w0_x = %s,\n' "$wx" >> ${INPUTFILE}
 printf ' w0_y = %s,\n' "$wy" >> ${INPUTFILE}
 printf ' a0 = %s,\n' "$a0" >> ${INPUTFILE}
 printf ' lam0 = %s,\n' "$lam0" >> ${INPUTFILE}
 printf ' lpx(1) = %s,\n' "$lx1" >> ${INPUTFILE}
 printf ' lpx(2) = %s,\n' "$lx2" >> ${INPUTFILE}
 printf ' lpx(3) = %s,\n' "$lx3" >> ${INPUTFILE}
 printf ' lpx(4) = %s,\n' "$lx4" >> ${INPUTFILE}
 printf ' lpx(5) = %s,\n' "$lx5" >> ${INPUTFILE}
 printf ' lpx(6) = %s,\n' "$lx6" >> ${INPUTFILE}
 printf ' lpx(7) = %s,\n' "$lx7" >> ${INPUTFILE}
 printf ' n_over_nc = %s,\n' "${n_nc}" >> ${INPUTFILE}
 printf ' n1_over_n = %s,\n' "${n1_nc}" >> ${INPUTFILE}
 printf ' n2_over_n = %s,\n' "${n2_nc}" >> ${INPUTFILE}
 printf ' w_sh = %s,\n' "${wnd_sh}" >> ${INPUTFILE}
 printf ' wi_time = %s,\n' "${w_in}" >> ${INPUTFILE}
 printf ' wf_time = %s,\n' "${w_end}" >> ${INPUTFILE}
 printf ' w_speed = %s,\n' "${w_speed}" >> ${INPUTFILE}
 printf '/' >> ${INPUTFILE}
 printf '\n\n' >> ${INPUTFILE}

 printf '&OUTPUT\n' >> ${INPUTFILE}
 printf ' nouts = %s,\n' "$nout" >> ${INPUTFILE}
 printf ' iene = %s,\n' "$iene" >> ${INPUTFILE}
 printf ' nvout = %s,\n' "$nf" >> ${INPUTFILE}
 printf ' nden = %s,\n' "$nd" >> ${INPUTFILE}
 printf ' npout = %s,\n' "$npv" >> ${INPUTFILE}
 printf ' nbout = %s,\n' "${end_p}" >> ${INPUTFILE}
 printf ' jump = %s,\n' "$jmp" >> ${INPUTFILE}
 printf ' pjump = %s,\n' "$pjmp" >> ${INPUTFILE}
 printf ' xp0_out = %s,\n' "$xp0" >> ${INPUTFILE}
 printf ' xp1_out = %s,\n' "$xp1" >> ${INPUTFILE}
 printf ' yp_out = %s,\n' "$ypmax" >> ${INPUTFILE}
 printf ' tmax = %s,\n' "$tmax" >> ${INPUTFILE}
 printf ' cfl = %s,\n' "$cfl" >> ${INPUTFILE}
 printf ' new_sim = %s,\n' "$new" >> ${INPUTFILE}
 printf ' id_new = %s,\n' "${id_ew}" >> ${INPUTFILE}
 printf ' dump = %s,\n' "$dump" >> ${INPUTFILE}
 printf ' npe_yz = %s,\n' "$pey" >> ${INPUTFILE}
 printf '/' >> ${INPUTFILE}
 printf '\n\n' >> ${INPUTFILE}


