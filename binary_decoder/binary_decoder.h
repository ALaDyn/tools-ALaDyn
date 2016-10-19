#pragma once

//#define ENABLE_DEBUG

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS

#define MAJOR_RELEASE  7
#define MINOR_RELEASE  1
#define BUGFIX_RELEASE 1

typedef float aladyn_float;

#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <functional>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <limits>
#include <cfloat>
#include <ios>
#include <cstdarg>

#if (defined CINECA)
#include <inttypes.h>
#include <stdint.h>
#endif

#if (!defined CINECA) && (defined _MSC_VER)
#include<cstdint>
#endif

#if (!defined CINECA) && (defined __GNUC__)
/* Test for GCC > 4.6.0 */
#if __GNUC__ > 4 || (__GNUC__ == 4 && (__GNUC_MINOR__ > 6))
#include<cstdint>
#else
#include <inttypes.h>
#include <stdint.h>
#endif
#endif

#if defined (_MSC_VER)
#include<wchar.h>
#define fseeko _fseeki64
#define ftello _ftelli64
#endif

#if defined (__MINGW32__)
#define fseeko fseeko64
#define ftello ftello64
#endif

#if defined (__CYGWIN__)
#define fseeko fseek
#define ftello ftell
#endif


#define MAX_NUM_OF_PARTICLES_IN_MEMORY 10000000            // in reality we store double this number -1
#define UMA_G                          1.660538921E-24     // from uma to grams
#define SPEED_OF_LIGHT                 2.99792458E+10      // cm / s
#define ME_G                           9.10938291E-28      // electron mass [g]
#define MP_G                           1.6726231E-24       // proton mass [g]
#define MP_MEV                         938.272013          // proton mass [MeV/c^2]
#define ME_MEV                         0.510998928         // electron mass [MeV/c^2]
#define MHI_UMA                        26.981538           // atomic weight of Aluminum in atomic mass units
#define MLI_UMA                        12.0107             // atomic weight of Carbon in atomic mass units
//#define CHARGE                       4.80320425e-10      // statC  - official value
#define CHARGE                         4.803262e-10        // statC  - Turchetti's value
#define FORCE_PRINTF_BUFFER_SIZE       1024
#define INVALID_COLUMN                -1
#define COLUMN_X                       0
#define COLUMN_Y                       1
#define COLUMN_Z                       2
#define COLUMN_PX                      3
#define COLUMN_PY                      4
#define COLUMN_PZ                      5
#define COLUMN_GAMMA                   6
#define COLUMN_THETA                   7
#define COLUMN_E                       8
#define COLUMN_THETAT                  9
#define COLUMN_TY                     10
#define COLUMN_TZ                     11
#define COLUMN_W                      12
#define COLUMN_CH                     13



