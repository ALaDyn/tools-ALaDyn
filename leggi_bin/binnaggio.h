#pragma once

#include "leggi_binario_ALaDyn_fortran.h"


struct _Binning
{
  _Binning(aladyn_float *, Parameters *, aladyn_float **, std::string, std::string);
  _Binning(aladyn_float *, Parameters *, aladyn_float *, std::string);
};


