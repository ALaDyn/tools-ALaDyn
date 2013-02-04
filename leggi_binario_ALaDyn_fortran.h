#ifndef __LEGGI_ALADYN_FORTRAN
#define __LEGGI_ALADYN_FORTRAN

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

class parametri;

#define NPARAMETRI	5
#define WEIGHT		0
#define FUNZIONE	1
#define SWAP		2
#define OUT_BINARY	3
#define OUT_ASCII	4


#if defined(_MSC_VER)
#include "leggi_campi.h"
#include "leggi_particelle.h"
#include "swap_tools.h"
#else
#include "leggi_campi.cpp"
#include "leggi_particelle.cpp"
#include "swap_tools.cpp"
#endif


#include "classeparametri.cpp"



#endif
