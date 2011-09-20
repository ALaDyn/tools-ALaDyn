#ifndef __LEGGI_ALADYN_FORTRAN
#define __LEGGI_ALADYN_FORTRAN

#include<cstdio>
#include<iostream>
#include<cstdlib>
#include<cmath>
#include<cstring>

#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))
#define TRUE 1
#define FALSE 0



void swap_endian_s(short* , int);
void swap_endian_i(int* , int );
void swap_endian_f(float* , int );
int leggi_particelle(char* , int , int , int , int);
int leggi_campi(char* , int , int , int , int);


#endif
