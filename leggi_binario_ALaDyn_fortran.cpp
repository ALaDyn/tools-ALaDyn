
#include<cstdio>
#include<iostream>
#include<cstdlib>
#include<cmath>
#include<cstring>

#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))
#define TRUE 1
#define FALSE 0

#include "leggi_binario_ALaDyn_fortran.h"
//#include "leggi_campi.cpp"
//#include "leggi_particelle.cpp"
//#include "swap_tools.cpp"


int main (int argc, char *argv[])
{
	int WEIGHT, funzione, FLAG_ENDIAN, out_swap, out_file;
	std::cout << "Programma di conversione da output binario Fortran ad output \"umano\"" << std::endl;
	if (argc == 1)
	{
		std::cout << "Usasi:  ./a.out  nomefiledaconvertire.bin  \n ";
		std::cout << "Per scriptare inserire nell'ordine WEIGHT funzione FLAG_ENDIAN out_swap out_file" << std::endl;
	}
	if (argc != 7)
	{
		std::cout << "Inserisci il WEIGHT: 0 per output vecchio, 1 per nuovo: ";
		std::cin >> WEIGHT;
		std::cout << "Che tipo di file e'? 1 campi, 2 particelle: ";
		std::cin >> funzione;
		std::cout << "Su che endian lavori? 1 big, 2 little: ";
		std::cin >> FLAG_ENDIAN;
		std::cout << "Vuoi fare lo swap dell'endian? 1 si', 0 no: ";
		std::cin >> out_swap;
		std::cout << "Vuoi un output completo? 1 si', 0 no: ";
		std::cin >> out_file;
	}
	else
	{
		WEIGHT = atoi(argv[2]);
		funzione = atoi(argv[3]);
		FLAG_ENDIAN = atoi(argv[4]);
		out_swap = atoi(argv[5]);
		out_file = atoi(argv[6]);
	}
	if (funzione == 1) leggi_campi(argv[1], WEIGHT, FLAG_ENDIAN, out_swap, out_file);
	else if (funzione == 2) leggi_particelle(argv[1], WEIGHT, FLAG_ENDIAN, out_swap, out_file);
	return 0;
}
