#ifndef __LEGGI_PARTICELLE_H
#define __LEGGI_PARTICELLE_H

#pragma warning(disable : 593)

#define _USE_MATH_DEFINES
#include "leggi_binario_ALaDyn_fortran.h"


#define MAX_LENGTH_FILENAME 200
#define P_MASS 938.272

int leggi_particelle(char* , int , int , int , int , parametri_binnaggio);

#endif

