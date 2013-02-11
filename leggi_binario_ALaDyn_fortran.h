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
#include<cstdarg>
#ifndef _WIN32
#include <strings.h>
#else
int    strcasecmp(const char* s1, const char* s2)
{
    for (;;) {
        int c1 = tolower( *((unsigned char*) s1++));
        int c2 = tolower( *((unsigned char*) s2++));

        if ((c1 != c2) || (c1 == '\0')) {
            return( c1 - c2);
        }
    }
}
#endif


// #define ENABLE_DEBUG


#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))
#define TRUE 1
#define FALSE 0

#define C						2.99792458e+10		// cm / s
#define ME_G					9.10938291E-28		// electron mass [g]
#define MP_G					1.6726231e-24		// proton mass [g]
#define MP_MEV					938.272013			// proton mass [MeV/c^2]
#define ME_MEV					0.510998928			// electron mass [MeV/c^2]
// #define CHARGE				4.80320425e-10		// statC	commentata perche' Turchetti la usa un po' diversa
#define CHARGE					4.803262e-10		// statC    valore usato da Turchetti; nb: e' impreciso negli ultimi due decimali
#define FROM_TESLA_TO_GAUSS		1.0e+4
// #define DA_ERG_A_MEV			6.2415097523028e+5	// conversione Servizi
#define DA_ERG_A_MEV			6.241509744512e+5	// conversione Sinigardi
#define FROM_VOLT_TO_STATVOLT	3.335640951982e-3	// 1 statvolt = 299.792458 volts.

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
