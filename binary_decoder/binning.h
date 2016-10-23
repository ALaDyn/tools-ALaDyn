#pragma once

#include "parameters.h"


struct _Binning
{
  _Binning(aladyn_float *, size_t, Parameters *, densityplot *);
  _Binning(aladyn_float * , size_t, Parameters * , histo * );
};


