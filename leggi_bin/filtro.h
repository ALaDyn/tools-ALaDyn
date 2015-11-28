#pragma once

#include "leggi_binario_ALaDyn_fortran.h"


// definizione numero filtri "abilitati"
#ifndef NUM_FILTRI
#define NUM_FILTRI               26
#endif

#define __0X00                   0x1
#define __0X01                   0x2
#define __0X02                   0x4
#define __0X03                   0x8
#define __0X04                   0x10
#define __0X05                   0x20
#define __0X06                   0x40
#define __0X07                   0x80
#define __0X08                   0x100
#define __0X09                   0x200
#define __0X10                   0x400
#define __0X11                   0x800
#define __0X12                   0x1000
#define __0X13                   0x2000
#define __0X14                   0x4000
#define __0X15                   0x8000
#define __0X16                   0x10000
#define __0X17                   0x20000
#define __0X18                   0x40000
#define __0X19                   0x80000
#define __0X20                   0x100000
#define __0X21                   0x200000
#define __0X22                   0x400000
#define __0X23                   0x800000
#define __0X24                   0x1000000
#define __0X25                   0x2000000
// fine filtri in uso, i prossimi sono codici liberi
#define __0X26                   0x4000000
#define __0X27                   0x8000000
#define __0X28                   0x10000000
#define __0X29                   0x20000000
#define __0X30                   0x40000000
#define __0X31                   0x80000000



struct _Filtro
{
  enum _Nomi
  {
    xmin, ymin, zmin, xmax, ymax, zmax,
    pxmin, pymin, pzmin, pxmax, pymax, pzmax,
    emin, emax, thetamin, thetamax, thetaTmin, thetaTmax,
    tymin, tymax, tzmin, tzmax, wmin, wmax, chmin, chmax
  } nomi;
  static float * costruisci_filtro(const char *, ...);
  static float * costruisci_filtro(Parametri *);
  static void individua_filtro(char *, float, float *&);
  static const unsigned int cost[];
  static unsigned int maschera_interna;
  struct _flag_filtri
  {
    unsigned meno_xmin : 1;
    unsigned meno_ymin : 1;
    unsigned meno_zmin : 1;
    unsigned piu_xmax : 1;
    unsigned piu_ymax : 1;
    unsigned piu_zmax : 1;
    unsigned meno_pxmin : 1;
    unsigned meno_pymin : 1;
    unsigned meno_pzmin : 1;
    unsigned piu_pxmax : 1;
    unsigned piu_pymax : 1;
    unsigned piu_pzmax : 1;
    unsigned meno_Emin : 1;
    unsigned piu_Emax : 1;
    unsigned meno_thetamin : 1;
    unsigned piu_thetamax : 1;
    unsigned meno_thetaTmin : 1;
    unsigned piu_thetaTmax : 1;
    unsigned meno_tymin : 1;
    unsigned meno_tzmin : 1;
    unsigned piu_tymax : 1;
    unsigned piu_tzmax : 1;
    unsigned meno_wmin : 1;
    unsigned piu_wmax : 1;
    unsigned meno_chmin : 1;
    unsigned piu_chmax : 1;
    _flag_filtri operator=(int)
    {
      meno_xmin = meno_ymin = meno_zmin =
        meno_pxmin = meno_pymin = meno_pzmin =
        piu_xmax = piu_ymax = piu_zmax =
        piu_pxmax = piu_pymax = piu_pzmax =
        meno_Emin = piu_Emax = meno_thetamin = piu_thetamax =
        meno_thetaTmin = piu_thetaTmax = meno_tymin =
        piu_tymax = meno_tzmin = piu_tzmax =
        meno_wmin = piu_wmax = meno_chmin = piu_chmax = 0;
      return *this;
    }
    // varie ed eventuali
  } flag_filtri;
  static const char * descr[];
  _Filtro(Parametri*, float *, unsigned int[], float *, unsigned int = 0);
};

