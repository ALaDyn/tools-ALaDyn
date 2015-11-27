#pragma once

#include "leggi_binario_ALaDyn_fortran.h"
#include "binnaggio.h"
#include "filtro.h"
#include "scrittura.h"
#include "swap_tools.h"


class Particella_v1 {
  double x, y, z;
  double px, py, pz;
};
class Particella_v2 {
  double x, y, z;
  double px, py, pz;
  double weight;
};
class Particella_v3 {
  double x, y, z;
  double px, py, pz;
  float weight;
  int charge;
};
union Particella_vX {
  Particella_v1 particella_v1;
  Particella_v2 particella_v2;
  Particella_v3 particella_v3;
};
union double_as_two_float {
  double d;
  float f[2]; //f[0] = peso, f[1] = carica
};




int leggi_particelle(Parametri * );


