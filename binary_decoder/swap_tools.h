#pragma once

#include <vector>
#include <cstddef>
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

int is_big_endian(void);
void swap_endian_s(short*, int);
void swap_endian_i(int*, int);
void swap_endian_i(std::vector<int> &);
void swap_endian_f(float*, size_t);
void swap_endian_f(float*, int);
void swap_endian_f(float***, size_t, size_t, size_t);
void swap_endian_f(std::vector<float> &);

