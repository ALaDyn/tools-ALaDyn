#pragma once

#include "leggi_binario_ALaDyn_fortran.h"


struct _Binnaggio
{
  _Binnaggio(float *, int, int, Parametri *, float **, std::string, std::string);
  _Binnaggio(float *, int, int, Parametri *, float *, std::string);
};


