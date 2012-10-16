#ifndef __LEGGI_ALADYN_FORTRAN
#define __LEGGI_ALADYN_FORTRAN

#include<cstdio>
#include<iostream>
#include<cstdlib>
#include<cmath>
#include<cstring>

#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))
#define TRUE 1
#define FALSE 0

#if defined(_MSC_VER)
#include "leggi_campi.h"
#include "leggi_particelle.h"
#include "swap_tools.h"
#else
#include "leggi_campi.cpp"
#include "leggi_particelle.cpp"
#include "swap_tools.cpp"
#endif

#endif
