
#include "leggi_binario_ALaDyn_fortran.h"
#include "swap_tools.h"


void swap_endian_s(short* in_s,int n)
{
	int i;
	union {short smio; char arr[2];}x;
	char buff;

	for(i=0;i<n;i++)
	{
		x.smio=in_s[i];
		buff=x.arr[0];
		x.arr[0]=x.arr[1];
		x.arr[1]=buff;
		in_s[i]=x.smio;
	}
}


void swap_endian_i(int* in_i,int n)
{
	int i;
	union {int imio; float fmio; char arr[4];}x;
	char buff;
	for(i=0;i<n;i++)
	{
		x.imio=in_i[i];
		buff=x.arr[0];
		x.arr[0]=x.arr[3];
		x.arr[3]=buff;
		buff=x.arr[1];
		x.arr[1]=x.arr[2];
		x.arr[2]=buff;
		in_i[i]=x.imio;
	}
}


void swap_endian_f(float* in_f, int n)
{
	int i;
	union {int imio; float fmio; char arr[4];}x;
	char buff;
	for(i=0;i<n;i++)
	{
		x.fmio=in_f[i];
		buff=x.arr[0];
		x.arr[0]=x.arr[3];
		x.arr[3]=buff;
		buff=x.arr[1];
		x.arr[1]=x.arr[2];
		x.arr[2]=buff;
		in_f[i]=x.fmio;
	}
}


