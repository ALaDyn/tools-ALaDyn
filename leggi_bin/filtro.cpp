# ifndef __FILTRO_CPP
# define __FILTRO_CPP
# include "leggi_binario_ALaDyn_fortran.h"

/*
_Filtro

è una struct orientata al filtraggio dei dati provenienti da simulazioni ALaDyn

Al momento sono contemplati filtri di soglia minima e massima per

CIASCUNA COORDINATA DI SPAZIO DI FASE e per ENERGIA

Quest'ultima è, al presente, valutata da _Filtro stessa, ma si può molto
facilmente introdurre modifiche per evitarne il calcolo, così come introdurre
altri filtraggi su altre variabili e secondo altri criteri (Ad Majora).

Questo file fornisce anche nomi mnemonici #define-iti per le costanti intere
con un solo bit acceso -- INDIPENDENTEMENTE QUINDI DA MSB o LSB ----

STRUTTURA DEL CODICE

Il namespace cost

definisce al suo interno nomi mnemonici per i filtri da applicare: ad esempio

cost :: xmin

serve ad attivare il filtro di soglia minima per la coordinata di fase X e

cost :: ymin | cost :: pzmax

attiva SIMULTANEAMENTE i filtri di soglia minima e massima per le coordinate di
fase X e Pz. Ciò viene fatto utilizzando alcune (non tutte) delle citate costanti
#define-ite.

cost :: tutte

raggruppa in un array ordinato tutte le precedenti costanti ed è utile per poter
gestire i loop.

DESCRIZIONE DELLA struct _Filtro

La enum class _Nomi serve a "battezzare" con nomi mnemonici di valori
CONSECUTIVI le stesse costanti del namespace cost; lo stesso fa l'array static
_Filtro :: cost ma usando degli unsigned int.

I metodi static

float * _Filtro :: costruisci_filtro
void _Filtro :: individua_filtro

lavorando in coppia, risparmiano al programmatore la grana di doversi ricordare
checchessia (nomi dei filtri, ordine di presenza negli array e quant'altro).

La variabile static maschera_interna e la struct flag_filtri sono trasparenti
all'utente finale.

L'array di stringhe _Filtro :: descr è utile per trasmettere i richiesti
argomenti alla funzione _Filtro :: costruisci_filtro.

Il costruttore di _Filtro fa tutto il lavoro: può essere eseguito fornendogli
quattro parametri oppure solo tre, lasciandogli l'arbitrio del quarto (ma solo
se viene usato come terzo parametro il valore reso da costruisci_filtro).

Il primo argomento da fornire al costruttore è, ovviamente, il puntatore ai
dati da filtrare; il secondo è un puntatore a due interi che siano,
nell'ordine, il numero di "particelle" e il numero di dati per particella.
Il terzo argomento è un array ordinato di valori float che contiene le diverse
soglie da utilizzare per la filtratura. Il quarto, se fornito, rappresenta i
filtri da applicare; se non fornito viene costruito da "costruisci_filtro"
assieme all'array da usare come terzo argomento.

*/


#ifndef _WIN32
#include <strings.h>
#else
int    strcasecmp(const char* s1, const char* s2)
{
  for (;;)
  {
    int c1 = tolower(*((unsigned char*)s1++));
    int c2 = tolower(*((unsigned char*)s2++));

    if ((c1 != c2) || (c1 == '\0'))
    {
      return(c1 - c2);
    }
  }
}
#endif

namespace cost
{
  unsigned int xmin = __0X00;
  unsigned int ymin = __0X01;
  unsigned int zmin = __0X02;
  unsigned int pxmin = __0X03;
  unsigned int pymin = __0X04;
  unsigned int pzmin = __0X05;
  unsigned int xmax = __0X06;
  unsigned int ymax = __0X07;
  unsigned int zmax = __0X08;
  unsigned int pxmax = __0X09;
  unsigned int pymax = __0X10;
  unsigned int pzmax = __0X11;
  unsigned int emin = __0X12;
  unsigned int emax = __0X13;
  unsigned int thetamin = __0X14;
  unsigned int thetamax = __0X15;
  unsigned int thetaTmin = __0X16;
  unsigned int thetaTmax = __0X17;
  unsigned int tymin = __0X18;
  unsigned int tymax = __0X19;
  unsigned int tzmin = __0X20;
  unsigned int tzmax = __0X21;
  unsigned int wmin = __0X22;
  unsigned int wmax = __0X23;
  unsigned int chmin = __0X24;
  unsigned int chmax = __0X25;

  unsigned int tutte[] =
  {
    xmin, ymin, zmin,
    pxmin, pymin, pzmin,
    xmax, ymax, zmax,
    pxmax, pymax, pzmax,
    emin, emax,
    thetamin, thetamax,
    thetaTmin, thetaTmax,
    tymin, tymax,
    tzmin, tzmax,
    wmin, wmax,
    chmin, chmax
  };
}


_Filtro::_Filtro(Parametri * parametri, float *dati, unsigned int n_dati[], float *val, unsigned int maschera)
{
  float * pntt_loc, p[] = { 0, 0, 0 }, E = 0., theta = 0., thetaT = 0., ty = 0., tz = 0., w = 0., ch = 0.;
  unsigned int corrente = 0, tests[32];
  bool flag;
  unsigned char tot_test = 0;

  flag_filtri = 0;
  if (!maschera) maschera = maschera_interna;
  if (!maschera)
  {
    return;
  }
  for (unsigned char c = 0; c < 32; ++c)
  {
    unsigned int r = maschera & cost[c];
    if (!r) continue;
    tests[tot_test++] = cost[c];
  }
  for (unsigned int i = 0; i < n_dati[0]; ++i)
  {
    pntt_loc = dati + i*n_dati[1];
    flag = true;

    if (((parametri->p[NCOLONNE] == 6 || parametri->p[NCOLONNE] == 7) && parametri->aladyn_version < 3) || (parametri->p[NCOLONNE] == 8 && parametri->aladyn_version == 3))
    {
      p[0] = pntt_loc[3], p[1] = pntt_loc[4], p[2] = pntt_loc[5];
      if (parametri->p[WEIGHT] && !parametri->overwrite_weight)
        w = pntt_loc[6];
      else
        w = parametri->overwrite_weight_value;
      if (parametri->aladyn_version == 3 && !parametri->overwrite_charge)
        ch = pntt_loc[7];
      else
        ch = parametri->overwrite_charge_value;
    }

    else if (((parametri->p[NCOLONNE] == 4 || parametri->p[NCOLONNE] == 5) && parametri->aladyn_version < 3) || (parametri->p[NCOLONNE] == 6 && parametri->aladyn_version == 3))
    {
      p[0] = pntt_loc[2], p[1] = pntt_loc[3], p[2] = 0.0;
      if (parametri->p[WEIGHT] && !parametri->overwrite_weight)
        w = pntt_loc[4];
      else
        w = parametri->overwrite_weight_value;
      if (parametri->aladyn_version == 3 && !parametri->overwrite_charge)
        ch = pntt_loc[5];
      else
        ch = parametri->overwrite_charge_value;
    }


    for (unsigned char c = 0; c < tot_test; ++c)
    {
      if (!flag) break;

      if (tests[c] == __0X12 || tests[c] == __0X13)
        E = (float)(parametri->massa_particella_MeV * (sqrtf(1.0f + p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) - 1.f));

      else if (tests[c] == __0X14 || tests[c] == __0X15)
        theta = (float)(atan2(sqrt(p[1] * p[1] + p[2] * p[2]), p[0])*180. / M_PI);

      else if (tests[c] == __0X16 || tests[c] == __0X17)
        thetaT = (float)atan(sqrt((p[1] * p[1] / (p[0] * p[0])) + (p[2] * p[2] / (p[0] * p[0]))));

      else if (tests[c] == __0X18 || tests[c] == __0X19)
        ty = p[1] / p[0];

      else if (tests[c] == __0X20 || tests[c] == __0X21)
        tz = p[2] / p[0];

      switch (tests[c])
      {
      case __0X00: // cost::xmin
        nomi = xmin;
        flag_filtri.meno_xmin = pntt_loc[(int)nomi] >= val[(int)nomi];
        flag = flag && flag_filtri.meno_xmin;
        break;
      case __0X01: // cost::ymin
        nomi = ymin;
        flag_filtri.meno_ymin = pntt_loc[(int)nomi] >= val[(int)nomi];
        flag = flag && flag_filtri.meno_ymin;
        break;
      case __0X02: // cost::zmin
        nomi = zmin;
        flag_filtri.meno_zmin = pntt_loc[(int)nomi] >= val[(int)nomi];
        flag = flag && flag_filtri.meno_zmin;
        break;
      case __0X03: // cost::pxmin
        nomi = pxmin;
        flag_filtri.meno_pxmin = pntt_loc[(int)nomi] >= val[(int)nomi];
        flag = flag && flag_filtri.meno_pxmin;
        break;
      case __0X04: // cost::pymin
        nomi = pymin;
        flag_filtri.meno_pymin = pntt_loc[(int)nomi] >= val[(int)nomi];
        flag = flag && flag_filtri.meno_pymin;
        break;
      case __0X05: // cost::pzmin
        nomi = pzmin;
        flag_filtri.meno_pzmin = pntt_loc[(int)nomi] >= val[(int)nomi];
        flag = flag && flag_filtri.meno_pzmin;
        break;
      case __0X06: // cost::xmax
        nomi = xmax;
        flag_filtri.piu_xmax = pntt_loc[(int)nomi - 6] <= val[(int)nomi];
        flag = flag && flag_filtri.piu_xmax;
        break;
      case __0X07: // cost::ymax
        nomi = ymax;
        flag_filtri.piu_ymax = pntt_loc[(int)nomi - 6] <= val[(int)nomi];
        flag = flag && flag_filtri.piu_ymax;
        break;
      case __0X08: // cost::zmax
        nomi = zmax;
        flag_filtri.piu_zmax = pntt_loc[(int)nomi - 6] <= val[(int)nomi];
        flag = flag && flag_filtri.piu_zmax;
        break;
      case __0X09: // cost::pxmax
        nomi = pxmax;
        flag_filtri.piu_pxmax = pntt_loc[(int)nomi - 6] <= val[(int)nomi];
        flag = flag && flag_filtri.piu_pxmax;
        break;
      case __0X10: // cost::pymax
        nomi = pymax;
        flag_filtri.piu_pymax = pntt_loc[(int)nomi - 6] <= val[(int)nomi];
        flag = flag && flag_filtri.piu_pymax;
        break;
      case __0X11: // cost::pzmax
        nomi = pzmax;
        flag_filtri.piu_pzmax = pntt_loc[(int)nomi - 6] <= val[(int)nomi];
        flag = flag && flag_filtri.piu_pzmax;
        break;
      case __0X12: // cost::emin
        nomi = emin;
        flag_filtri.meno_Emin = E >= val[12];
        flag = flag && flag_filtri.meno_Emin;
        break;
      case __0X13:  // cost::emax
        nomi = emax;
        flag_filtri.piu_Emax = E <= val[13];
        flag = flag && flag_filtri.piu_Emax;
        break;
      case __0X14: // cost::thetamin
        nomi = thetamin;
        flag_filtri.meno_thetamin = theta >= val[14];
        flag = flag && flag_filtri.meno_thetamin;
        break;
      case __0X15: // cost::thetamax
        nomi = thetamax;
        flag_filtri.piu_thetamax = theta <= val[15];
        flag = flag && flag_filtri.piu_thetamax;
        break;
      case __0X16: // cost::thetaTmin
        nomi = thetaTmin;
        flag_filtri.meno_thetaTmin = thetaT >= val[16];
        flag = flag && flag_filtri.meno_thetaTmin;
        break;
      case __0X17: // cost::thetaTmax
        nomi = thetaTmax;
        flag_filtri.piu_thetaTmax = thetaT <= val[17];
        flag = flag && flag_filtri.piu_thetaTmax;
        break;
      case __0X18: // cost::tymin
        nomi = tymin;
        flag_filtri.meno_tymin = ty >= val[18];
        flag = flag && flag_filtri.meno_tymin;
        break;
      case __0X19: // cost::tymax
        nomi = tymax;
        flag_filtri.piu_tymax = ty <= val[19];
        flag = flag && flag_filtri.piu_tymax;
        break;
      case __0X20: // cost::tzmin
        nomi = tzmin;
        flag_filtri.meno_tzmin = tz >= val[20];
        flag = flag && flag_filtri.meno_tzmin;
        break;
      case __0X21: // cost::tzmax
        nomi = tzmax;
        flag_filtri.piu_tzmax = tz <= val[21];
        flag = flag && flag_filtri.piu_tzmax;
        break;
      case __0X22: // cost::wmin
        nomi = wmin;
        flag_filtri.meno_wmin = w >= val[22];
        flag = flag && flag_filtri.meno_wmin;
        break;
      case __0X23: // cost::wmax
        nomi = wmax;
        flag_filtri.piu_wmax = w <= val[23];
        flag = flag && flag_filtri.piu_wmax;
        break;
      case __0X24: // cost::chmin
        nomi = chmin;
        flag_filtri.meno_chmin = ch >= val[24];
        flag = flag && flag_filtri.meno_chmin;
        break;
      case __0X25: // cost::chmax
        nomi = chmax;
        flag_filtri.piu_chmax = ch <= val[25];
        flag = flag && flag_filtri.piu_chmax;
        break;
      }

    }
    if (!flag) continue;
    for (unsigned int k = 0; k < n_dati[1]; ++k) dati[n_dati[1] * corrente + k] = pntt_loc[k];
    corrente++;
  }
  n_dati[0] = corrente;
  maschera_interna = 0;
}

const char * _Filtro::descr[] =
{
  "+xmin",
  "+ymin",
  "+zmin",
  "+xmax",
  "+ymax",
  "+zmax",
  "+pxmin",
  "+pymin",
  "+pzmin",
  "+pxmax",
  "+pymax",
  "+pzmax",
  "+Emin",
  "+Emax",
  "+thetamin",
  "+thetamax",
  "+thetaTmin",
  "+thetaTmax",
  "+tymin",
  "+tymax",
  "+tzmin",
  "+tzmax",
  "+wmin",
  "+wmax",
  "+chmin",
  "+chmax"
  // varie ed eventuali
};

const unsigned int _Filtro::cost[] =
{
  __0X00, __0X01, __0X02, __0X03, __0X04, __0X05, __0X06, __0X07,
  __0X08, __0X09, __0X10, __0X11, __0X12, __0X13, __0X14, __0X15,
  __0X16, __0X17, __0X18, __0X19, __0X20, __0X21, __0X22, __0X23,
  __0X24, __0X25, __0X26, __0X27, __0X28, __0X29, __0X30, __0X31
};

float * _Filtro::costruisci_filtro(Parametri * parametri)
{
  char ** miei_args;
  float * miei_val;
  int indices[NUM_FILTRI], quanti = 0;
  for (int i = 1; i < parametri->argc; ++i)
  {
    if (parametri->argv[i][0] == '+')
      indices[quanti++] = i;
  }
  indices[quanti] = -1;
  if (!quanti) return (float *)NULL;
  miei_args = new char *[NUM_FILTRI + 1], miei_args[NUM_FILTRI] = 0;
  miei_val = new float[NUM_FILTRI + 1];
  for (int i = 0; i < quanti; ++i)
    miei_args[i] = const_cast<char*>(parametri->argv[indices[i]].c_str()),
    miei_val[i] = (float)atof(parametri->argv[indices[i] + 1].c_str());
  miei_args[quanti] = 0;
  return costruisci_filtro(
    miei_args[0], miei_val[0],
    miei_args[1], miei_val[1],
    miei_args[2], miei_val[2],
    miei_args[3], miei_val[3],
    miei_args[4], miei_val[4],
    miei_args[5], miei_val[5],
    miei_args[6], miei_val[6],
    miei_args[7], miei_val[7],
    miei_args[8], miei_val[8],
    miei_args[9], miei_val[9],
    miei_args[10], miei_val[10],
    miei_args[11], miei_val[11],
    miei_args[12], miei_val[12],
    miei_args[13], miei_val[13],
    miei_args[14], miei_val[14],
    miei_args[15], miei_val[15],
    miei_args[16], miei_val[16],
    miei_args[17], miei_val[17],
    miei_args[18], miei_val[18],
    miei_args[19], miei_val[19],
    miei_args[20], miei_val[20],
    miei_args[21], miei_val[21],
    miei_args[22], miei_val[22],
    miei_args[23], miei_val[23],
    miei_args[24], miei_val[24],
    miei_args[25], miei_val[25],
    miei_args[NUM_FILTRI]);
}


float * _Filtro::costruisci_filtro(const char *p, ...)
{
  va_list app;
  char * buff = new char[MAX_LENGTH_FILENAME], *z = buff;
  float val, *tutti_val = new float[NUM_FILTRI];
  va_start(app, p);
  strcpy(buff, p);
  do
  {
    val = (float)va_arg(app, double);
    individua_filtro(buff, val, tutti_val);
    if (!tutti_val) return (float*)NULL;
    z = va_arg(app, char *);
    if (!z) break;
    strcpy(buff, z);
  } while (1);
  va_end(app);
  return tutti_val;
}

void _Filtro::individua_filtro(char *b, float v, float *& V)
{
  int i;
  for (i = 0; i < NUM_FILTRI; ++i) if (!strcasecmp(b, descr[i])) break;
  if (i >= NUM_FILTRI)
  {
    V = (float*)NULL; // non è detto che sia la cosa migliore da farsi
    return;
  }
  V[i] = v;
  maschera_interna |= cost::tutte[i];
}

unsigned int _Filtro::maschera_interna = 0;


# endif
