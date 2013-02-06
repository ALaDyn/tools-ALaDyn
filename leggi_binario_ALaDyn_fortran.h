#ifndef __LEGGI_ALADYN_FORTRAN
#define __LEGGI_ALADYN_FORTRAN

#define _USE_MATH_DEFINES

#include<cstdio>
#include<iostream>
#include<cstdlib>
#include<cmath>
#include<cstring>
#include<string>
#include<fstream>
#include<sstream>
#include<iomanip>

#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))
#define TRUE 1
#define FALSE 0
#define P_MASS 938.272

class parametri;

#define NUMERO_MASSIMO	1.0e30
#define MAX_LENGTH_FILENAME 200

#pragma warning(disable : 593)

#define NPARAMETRI	8
#define WEIGHT		0
#define FUNZIONE	1
#define SWAP		2
#define OUT_BINARY	3
#define OUT_ASCII	4
#define	FIND_MINMAX	5
#define DO_BINNING	6
#define OUT_PARAMS	7


#include "classeparametri.cpp"

#if defined(_MSC_VER)
#include "swap_tools.h"
#include "leggi_campi.h"
#include "leggi_particelle.h"
#else
#include "swap_tools.cpp"
#include "leggi_campi.cpp"
#include "leggi_particelle.cpp"
#endif



#endif
