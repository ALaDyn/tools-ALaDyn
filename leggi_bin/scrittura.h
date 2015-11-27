#pragma once

#include "leggi_binario_ALaDyn_fortran.h"

struct _Scrittura
{
  _Scrittura(Parametri *, float **, std::string, std::string, std::string);
  _Scrittura(Parametri *, float *, std::string, std::string);
};


