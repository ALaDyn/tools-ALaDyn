#pragma once

#include "leggi_binario_ALaDyn_fortran.h"

struct _Write
{
  _Write(Parameters *, aladyn_float **, std::string, std::string, std::string);
  _Write(Parameters *, aladyn_float *, std::string, std::string);
};


