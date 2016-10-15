#pragma once

#include "binary_decoder.h"


struct _Binning
{
  _Binning(aladyn_float *, Parameters *, aladyn_float **, std::string, std::string);
  _Binning(aladyn_float *, Parameters *, aladyn_float *, std::string);
};


