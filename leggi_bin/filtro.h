#pragma once

#include "leggi_binario_ALaDyn_fortran.h"

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

