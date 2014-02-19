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


NCPU=12
ENABLE_FINAL_DUMPS=1


##### nx, ny,nz
## nb: il numero punti di griglia in y e z, se diversi da 1, deve essere divisibile per il numero di processori
## nb2: se si abilita lo stretching, il numero di punti in y e z deve anche essere divisibile per 6, dato che la
##      zona di stretching e' definita come 1/6 dei punti totali
## nb3: per ora e' meglio mettere in z un numero di punti uguale a quello in y, a meno che non si metta 1 per le simulazioni 2D
NUMERO_PUNTI_GRIGLIA_X=2640
NUMERO_PUNTI_GRIGLIA_Y=2400
##NB: mettere il seguente numero ad 1 per fare simulazioni 2D
NUMERO_PUNTI_GRIGLIA_Z=1
##NB: a causa di un baco, mettere il seguente valore superiore a (NUMERO_PUNTI_GRIGLIA_Y - NUMERO_PUNTI_GRIGLIA_Y / NCPU)
NUMERO_PUNTI_GRIGLIA_TRASVERSALI_OCCUPATI_DAL_PLASMA=2280


##nb: ricordarsi di risolvere la skin depth con k0
##### k0,yx_rat
NUMERO_PUNTI_GRIGLIA_PER_MICRON="80.0"
RAPPORTO_DIMENSIONE_GRIGLIA_TRASVERSA_GRIGLIA_LONGITUDINALE="2.0"

##### a questo punto sappiamo le dimensioni in micron della griglia 
#####  Dimx = nx/k0    Dimy = ny*yx_rat/k0

##### LPf,Der,str,iform
ORDINE_INTEGRAZIONE_LEAPFROG=2
TIPO_DERIVATA=3
##NB: per il seguente dato, 0=niente, 1=stretching trasversale, 2=stretching anche lungo x, 3=PML
TIPO_BOUNDARIES=1
##NB: per il seguente dato, 1=conservazione carica, 2=conservazione energia
ALGORITMO_INTEGRAZIONE=2


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
PLASMA_MODEL=2


#### nsp,nsb,Z_i,A_i 
##NB: se il plasma model e' 1, nsp deve essere 2 oppure 1: il significato e' che se e' 1 i protoni sono fissi (genero solo gli elettroni), 
##                                                         se e' 2 invece genero entrambe le specie
##    se il plasma model e' 2, nsp deve essere uguale a 4
##    se il plasma model e' 4, riferirsi alla guida in PLASMA_MODEL
NUMERO_SPECIE=4
##NB: per la documentazione riguardante il seguente numero, riferirsi alla spiegazione di PLASMA_MODEL
NUMERO_SPECIE_SECONDARIO=0
##NB: per i seguenti valori di solito si usa ionizzazione 9 e peso 27 (alluminio)
NUMERO_IONIZZAZIONE_Z=9
NUMERO_DI_NUCLEONI_A=27


## Ricordarsi di risolvere la densita' in unita' di densita' critiche con il numero di particelle per cella
#### np_xc
NUMERO_ELETTRONI_LONGITUDINALMENTE_PER_CELLA_LAYER_CENTRALE=6
NUMERO_IONI_LONGITUDINALMENTE_PER_CELLA_LAYER_CENTRALE=2
NUMERO_ELETTRONI_LONGITUDINALMENTE_PER_CELLA_LAYER_FRONTALE=2
NUMERO_IONI_LONGITUDINALMENTE_PER_CELLA_LAYER_FRONTALE=2
NUMERO_ELETTRONI_LONGITUDINALMENTE_PER_CELLA_LAYER_POSTERIORE=2
NUMERO_IONI_LONGITUDINALMENTE_PER_CELLA_LAYER_POSTERIORE=2


#### np_yc
NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA_LAYER_CENTRALE=6
NUMERO_IONI_TRASVERSALMENTE_PER_CELLA_LAYER_CENTRALE=2
NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA_LAYER_FRONTALE=2
NUMERO_IONI_TRASVERSALMENTE_PER_CELLA_LAYER_FRONTALE=2
NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA_LAYER_POSTERIORE=2
NUMERO_IONI_TRASVERSALMENTE_PER_CELLA_LAYER_POSTERIORE=2


#### a questo punto e' meglio controllare che np_xc_el_l2*np_yc_el_l2 sia
#### uguale a Z_i*np_xc_ion_l2*np_yc_ion_l2
#### NON DOVREBBE ESSER NECESSARIO, MA PASQUALE INSISTE SU QUESTO PUNTO:
#### E' PIU' TRANQUILLO SE TUTTE LE SPECIE HANNO LO STESSO PESO


#### t0, xc, wx, wy, a0,lam0 ## tutti in micrometri
POSIZIONE_INIZIALE_PICCO_IMPULSO_LASER="10.7"
DISTANZA_INIZIALE_PICCO_IMPULSO_LASER_DAL_FUOCO="11."
LUNGHEZZA_LASER_TOTALE="20.6"
WAIST_LASER="3."
PARAMETRO_ADIMENSIONALE_LASER_A0="5."
LUNGHEZZA_ONDA_LASER="0.8"



#### lx(1:7)
SPESSORE_LAYER_FRONTALE="0.0"
SPESSORE_RAMPA_LAYER_FRONTALE_LAYER_CENTRALE="0.0"
SPESSORE_LAYER_CENTRALE="0.8"
SPESSORE_RAMPA_LAYER_CENTRALE_LAYER_POSTERIORE="0.0"
SPESSORE_LAYER_POSTERIORE="0.08"
ANGOLO_ROTAZIONE_LASER="0.0"
OFFSET_FINE_LASER_INIZIO_TARGHETTA="0.1"


#### a questo punto sappiamo la posizione della targhetta
#### inizio_targhetta = xc + wx/2 + lx7


#### n/nc, n1/nc,n2/nc (tutte le densita' sono quindi rapportate alla densita' critica)
DENSITA_ELETTRONI_LAYER_FRONTALE="0.0"
DENSITA_ELETTRONI_LAYER_CENTRALE="50."
DENSITA_ELETTRONI_LAYER_POSTERIORE="9."


#### nf, nd, npv, end_p
####
TIPO_OUTPUT_CAMPI=3
TIPO_OUTPUT_DENSITA_GRIGLIA=1
TIPO_OUTPUT_SPAZIOFASI_PARTICELLE=2
#### il flag seguente, se vale 1, impone un taglio nello spazio xyz secondo le dimensioni descritte dai
#### parametri imposti nelle tre seguenti righe, al fine di ridurre le dimensioni dei file.
FLAG_TAGLIO_OUTPUT_SPAZIOFASI_PARTICELLE=1
#### X1-X0 seguenti definiscono quindi l'altezza del parallelepipedo (o rettangolo in 2D) di spazio da dumpare
X0_TAGLIO_OUTPUT="24.0"
X1_TAGLIO_OUTPUT="60.0"
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
MW_START_TIME="10."
#### time step fine movimento della window
MW_END_TIME="50."
#### velocità beta con cui si muove la window
MW_SPEED="0.5"



#### Non sono anticipati qui i parametri del bunch di elettroni per l'iniezione
#### Saranno aggiunti in futuro



#### I seguenti valori non sono documentati (ibeam, cfl)
BEAM_MODEL=1
CFL="0.85"


########################################
########################################
## FINE PARAMETRI - NON TOCCARE OLTRE ##
########################################
########################################
nx=${NUMERO_PUNTI_GRIGLIA_X}
ny=${NUMERO_PUNTI_GRIGLIA_Y}
nz=${NUMERO_PUNTI_GRIGLIA_Z}
nplasma=${NUMERO_PUNTI_GRIGLIA_TRASVERSALI_OCCUPATI_DAL_PLASMA}
Line01Comments=" nx, ny,nz              !the grid dimensions"
k0=${NUMERO_PUNTI_GRIGLIA_PER_MICRON}
yx_rat=${RAPPORTO_DIMENSIONE_GRIGLIA_TRASVERSA_GRIGLIA_LONGITUDINALE}
Line02Comments=" k0,yx_rat             !number of grid points/micron Dx=k0, Dy/Dx=yx_rat "
LPf=${ORDINE_INTEGRAZIONE_LEAPFROG}
Der=${TIPO_DERIVATA}
str=${TIPO_BOUNDARIES}
iform=${ALGORITMO_INTEGRAZIONE}
Line03Comments=" LPf,Der,str,iform    !scheme order in time and space, type of boundaries"
mdl=${LASER_MODEL}
dmdl=${PLASMA_MODEL}
ibeam=${BEAM_MODEL}
Line04Comments=" mdl,dmdl,ibeam       ! laser model, initial plasma model"
nsp=${NUMERO_SPECIE}
nsb=${NUMERO_SPECIE_SECONDARIO}
Z_i=${NUMERO_IONIZZAZIONE_Z}
A_i=${NUMERO_DI_NUCLEONI_A}
Line05Comments=" nsp,nsb,Z_i,A_i            ! number of species, (Z,A) ions (if nsp=3)"
np_xc_el_l2=${NUMERO_ELETTRONI_LONGITUDINALMENTE_PER_CELLA_LAYER_CENTRALE}
np_xc_ion_l2=${NUMERO_IONI_LONGITUDINALMENTE_PER_CELLA_LAYER_CENTRALE}
np_xc_el_l1=${NUMERO_ELETTRONI_LONGITUDINALMENTE_PER_CELLA_LAYER_FRONTALE}
np_xc_ion_l1=${NUMERO_IONI_LONGITUDINALMENTE_PER_CELLA_LAYER_FRONTALE}
np_xc_el_l3=${NUMERO_ELETTRONI_LONGITUDINALMENTE_PER_CELLA_LAYER_POSTERIORE}
np_xc_ion_l3=${NUMERO_IONI_LONGITUDINALMENTE_PER_CELLA_LAYER_POSTERIORE}
Line06Comments=" np_xc                  !np_per_xcell (el,ion) in layers (2,1,3)"
np_yc_el_l2=${NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA_LAYER_CENTRALE}
np_yc_ion_l2=${NUMERO_IONI_TRASVERSALMENTE_PER_CELLA_LAYER_CENTRALE}
np_yc_el_l1=${NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA_LAYER_FRONTALE}
np_yc_ion_l1=${NUMERO_IONI_TRASVERSALMENTE_PER_CELLA_LAYER_FRONTALE}
np_yc_el_l3=${NUMERO_ELETTRONI_TRASVERSALMENTE_PER_CELLA_LAYER_POSTERIORE}
np_yc_ion_l3=${NUMERO_IONI_TRASVERSALMENTE_PER_CELLA_LAYER_POSTERIORE}
Line07Comments=" np_yc                  !np_per_ycell (el,ion) in layers (2,1,3)"
t0=${POSIZIONE_INIZIALE_PICCO_IMPULSO_LASER}
xc=${DISTANZA_INIZIALE_PICCO_IMPULSO_LASER_DAL_FUOCO}
wx=${LUNGHEZZA_LASER_TOTALE}
wy=${WAIST_LASER}
a0=${PARAMETRO_ADIMENSIONALE_LASER_A0}
lam0=${LUNGHEZZA_ONDA_LASER}
Line08Comments=" t0, xc, wx, wy, a0,lam0         ! laser: xf=xc+t0 focus, wx=length,wy=waist"
bch="0.0"
xb="50"
gam0="50."
sx="30."
sy="50."
epsx="0."
epsy="0."
dg="0."
Line09Comments=" bch,xb,gam0,sx,sy,epsx,epsy,dg        ! el bunch param "
wch="0."
xw="50."
sxw="10."
syw="10."
epsxw="0."
epsyw="0."
dgw="0."
Line10Comments=" wch,xw,sxw,syw,epsxw,epsyw,dgw        ! el bunch param "
lx1=${SPESSORE_LAYER_FRONTALE}
lx2=${SPESSORE_RAMPA_LAYER_FRONTALE_LAYER_CENTRALE}
lx3=${SPESSORE_LAYER_CENTRALE}
lx4=${SPESSORE_RAMPA_LAYER_CENTRALE_LAYER_POSTERIORE}
lx5=${SPESSORE_LAYER_POSTERIORE}
lx6=${ANGOLO_ROTAZIONE_LASER}
lx7=${OFFSET_FINE_LASER_INIZIO_TARGHETTA}
Line11Comments=" lx(1:7)                  !1 layer 2 ramp  3 central layer 4 ramp  5 end layer  6 angle  7 offset"
n_nc=${DENSITA_ELETTRONI_LAYER_CENTRALE}
n1_nc=${DENSITA_ELETTRONI_LAYER_FRONTALE}
n2_nc=${DENSITA_ELETTRONI_LAYER_POSTERIORE}
Line12Comments=" n/nc, n1/nc,n2/nc         !the el density in layer 3,1,2"
wnd_sh=${MW_CALL_EVERY_N_TIMESTEPS}
w_in=${MW_START_TIME}
w_end=${MW_END_TIME}
w_speed=${MW_SPEED}
Line13Comments=" wnd_sh,w_in,w_end,w_speed !mov.wind.:points, init. time,final time,speed"
nout=$4
iene=$5
nf=${TIPO_OUTPUT_CAMPI}
nd=${TIPO_OUTPUT_DENSITA_GRIGLIA}
npv=${TIPO_OUTPUT_SPAZIOFASI_PARTICELLE}
end_p=${FLAG_TAGLIO_OUTPUT_SPAZIOFASI_PARTICELLE}
Line14Comments=" nout,iene,nf,nd,npv,end_p !out numb nf fields, nd den fields npv part"
jmp=${JUMP_GRIGLIA}
pjmp=${JUMP_PARTICELLE}
Line15Comments=" jmp, pjmp          ! jump rate in field grid data pjmp= rate for particles "
xp0=${X0_TAGLIO_OUTPUT}
xp1=${X1_TAGLIO_OUTPUT}
ypmax=${SEMILATO_BASE_TAGLIO_OUTPUT}
Line16Comments=" xp0, xp1, ypmax        !volume for out particles "
tmax=$3
cfl=${CFL}
Line17Comments=" tmax, cfl               !"
new=$1
id_ew=$2
dump=${ENABLE_FINAL_DUMPS}
pey=${NCPU}
Line18Comments=" new,id_ew,dump,pey"


 INPUTFILE=input.data

 rm  ${INPUTFILE}
 touch ${INPUTFILE}

 printf '%s, %s, %s, %s\n' "$nx" "$ny" "$nz" "$nplasma" >> ${INPUTFILE}
 printf '%s\n' "${Line01Comments}" >> ${INPUTFILE}
 printf '%s, %s\n' "$k0" "$yx_rat" >> ${INPUTFILE}
 printf '%s\n' "${Line02Comments}" >> ${INPUTFILE}
 printf '%s, %s, %s, %s\n' "$LPf" "$Der" "$str" "$iform" >> ${INPUTFILE}
 printf '%s\n' "${Line03Comments}" >> ${INPUTFILE}
 printf '%s, %s, %s\n' "$mdl" "$dmdl" "$ibeam" >> ${INPUTFILE}
 printf '%s\n' "${Line04Comments}" >> ${INPUTFILE}
 printf '%s, %s, %s, %s\n' "$nsp" "$nsb" "${Z_i}" "${A_i}" >> ${INPUTFILE}
 printf '%s\n' "${Line05Comments}" >> ${INPUTFILE}
 printf '%s, %s, %s, %s, %s, %s\n' "${np_xc_el_l2}" "${np_xc_ion_l2}" "${np_xc_el_l1}" "${np_xc_ion_l1}" "${np_xc_el_l3}" "${np_xc_ion_l3}" >> ${INPUTFILE}
 printf '%s\n' "${Line06Comments}" >> ${INPUTFILE}
 printf '%s, %s, %s, %s, %s, %s\n' "${np_yc_el_l2}" "${np_yc_ion_l2}" "${np_yc_el_l1}" "${np_yc_ion_l1}" "${np_yc_el_l3}" "${np_yc_ion_l3}" >> ${INPUTFILE}
 printf '%s\n' "${Line07Comments}" >> ${INPUTFILE}
 printf '%s, %s, %s, %s, %s, %s\n' "$t0" "$xc" "$wx" "$wy" "$a0" "$lam0" >> ${INPUTFILE}
 printf '%s\n' "${Line08Comments}" >> ${INPUTFILE}
 printf '%s, %s, %s, %s, %s, %s, %s, %s\n' "$bch" "$xb" "$gam0" "$sx" "$sy" "$epsx" "$epsy" "$dg" >> ${INPUTFILE}
 printf '%s\n' "${Line09Comments}" >> ${INPUTFILE}
 printf '%s, %s, %s, %s, %s, %s, %s\n' "$wch" "$xw" "$sxw" "$syw" "$epsxw" "$epsyw" "$dgw" >> ${INPUTFILE}
 printf '%s\n' "${Line10Comments}" >> ${INPUTFILE}
 printf '%s, %s, %s, %s, %s, %s, %s\n' "$lx1" "$lx2" "$lx3" "$lx4" "$lx5" "$lx6" "$lx7" >> ${INPUTFILE}
 printf '%s\n' "${Line11Comments}" >> ${INPUTFILE}
 printf '%s, %s, %s\n' "${n_nc}" "${n1_nc}" "${n2_nc}" >> ${INPUTFILE}
 printf '%s\n' "${Line12Comments}" >> ${INPUTFILE}
 printf '%s, %s, %s, %s\n' "${wnd_sh}" "${w_in}" "${w_end}" "${w_speed}" >> ${INPUTFILE}
 printf '%s\n' "${Line13Comments}" >> ${INPUTFILE}
 printf '%s, %s, %s, %s, %s, %s\n' "$nout" "$iene" "$nf" "$nd" "$npv" "${end_p}" >> ${INPUTFILE}
 printf '%s\n' "${Line14Comments}" >> ${INPUTFILE}
 printf '%s, %s\n' "$jmp" "$pjmp" >> ${INPUTFILE}
 printf '%s\n' "${Line15Comments}" >> ${INPUTFILE}
 printf '%s, %s, %s\n' "$xp0" "$xp1" "$ypmax" >> ${INPUTFILE}
 printf '%s\n' "${Line16Comments}" >> ${INPUTFILE}
 printf '%s, %s\n' "$tmax" "$cfl" >> ${INPUTFILE}
 printf '%s\n' "${Line17Comments}" >> ${INPUTFILE}
 printf '%s, %s, %s, %s\n' "$new" "${id_ew}" "$dump" "$pey" >> ${INPUTFILE}
 printf '%s\n' "${Line18Comments}" >> ${INPUTFILE}

