#pragma once

#include "parameters.h"


//number of enabled filters
#ifndef ENABLED_FILTERS
#define ENABLED_FILTERS               26
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
//from here going on we have spare filter codes
#define __0X26                   0x4000000
#define __0X27                   0x8000000
#define __0X28                   0x10000000
#define __0X29                   0x20000000
#define __0X30                   0x40000000
#define __0X31                   0x80000000



struct _Filter
{
  enum _Nomi
  {
    xmin, ymin, zmin, xmax, ymax, zmax,
    pxmin, pymin, pzmin, pxmax, pymax, pzmax,
    emin, emax, thetamin, thetamax, thetaTmin, thetaTmax,
    tymin, tymax, tzmin, tzmax, wmin, wmax, chmin, chmax
  } name;
  static aladyn_float * build_filter(const char *, ...);
  static aladyn_float * build_filter(Parameters *);
  static void find_filter(char *, aladyn_float, aladyn_float *&);
  static const unsigned int cost[];
  static unsigned int internal_mask;
  struct _flag_filters
  {
    unsigned minus_xmin : 1;
    unsigned minus_ymin : 1;
    unsigned minus_zmin : 1;
    unsigned plus_xmax : 1;
    unsigned plus_ymax : 1;
    unsigned plus_zmax : 1;
    unsigned minus_pxmin : 1;
    unsigned minus_pymin : 1;
    unsigned minus_pzmin : 1;
    unsigned plus_pxmax : 1;
    unsigned plus_pymax : 1;
    unsigned plus_pzmax : 1;
    unsigned minus_Emin : 1;
    unsigned plus_Emax : 1;
    unsigned minus_thetamin : 1;
    unsigned plus_thetamax : 1;
    unsigned minus_thetaTmin : 1;
    unsigned plus_thetaTmax : 1;
    unsigned minus_tymin : 1;
    unsigned minus_tzmin : 1;
    unsigned plus_tymax : 1;
    unsigned plus_tzmax : 1;
    unsigned minus_wmin : 1;
    unsigned plus_wmax : 1;
    unsigned minus_chmin : 1;
    unsigned plus_chmax : 1;
    _flag_filters operator=(int)
    {
      minus_xmin = minus_ymin = minus_zmin =
        minus_pxmin = minus_pymin = minus_pzmin =
        plus_xmax = plus_ymax = plus_zmax =
        plus_pxmax = plus_pymax = plus_pzmax =
        minus_Emin = plus_Emax = minus_thetamin = plus_thetamax =
        minus_thetaTmin = plus_thetaTmax = minus_tymin =
        plus_tymax = minus_tzmin = plus_tzmax =
        minus_wmin = plus_wmax = minus_chmin = plus_chmax = 0;
      return *this;
    }
    // varie ed eventuali
  } flag_filters;
  static const char * descr[];
  _Filter(Parameters*, aladyn_float *, unsigned int[], aladyn_float *, unsigned int = 0);
};

