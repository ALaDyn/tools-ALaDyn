

#include "leggi_binario_ALaDyn_fortran.h"


int main (int argc, char *argv[])
{
	int WEIGHT, funzione, out_swap, out_binary, out_ascii;
//	int FLAG_ENDIAN;
	std::cout << "Programma di conversione da output binario Fortran ad output \"umano\"" << std::endl;
	std::cout << "Inserisci il WEIGHT: 0 per output vecchio, 1 per nuovo: ";
	std::cin >> WEIGHT;
	std::cout << "Che tipo di file e'? 1 campi, 2 particelle: ";
	std::cin >> funzione;
//	std::cout << "Su che endian lavori? 1 big, 2 little: ";
//	std::cin >> FLAG_ENDIAN;
	std::cout << "E' necessario fare lo swap dell'endian? 1 si', 0 no: ";
	std::cin >> out_swap;
	if (funzione == 2)
	{
		std::cout << "Vuoi l'output binario pulito dal fortran? 1 si', 0 no: ";
		std::cin >> out_binary;
		std::cout << "Vuoi l'output dello spazio delle fasi ascii? 1 si', 0 no: ";
		std::cin >> out_ascii;
	}
	if (funzione == 1) leggi_campi(argv[1], out_swap);
	else if (funzione == 2) leggi_particelle(argv[1], WEIGHT, out_swap, out_binary, out_ascii);
	return 0;
}
