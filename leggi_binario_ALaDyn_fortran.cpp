

#include "leggi_binario_ALaDyn_fortran.h"


int main (int argc, char *argv[])
{
	int WEIGHT, funzione, out_swap, out_binary, out_ascii;
	if (argc == 1)
	{
		std::cout << "Si usa:\n-manuale: ./reader Nomefile.bin" << std::endl;
		std::cout << "-batch: ./reader Nomefile.bin type(1|2) swap(0|1) [WEIGHT(0|1) cleanBin(0|1) asciiOut(0|1)]" << std::endl;
		std::cout << "parameters inside square brackets are required only for particle-type binary files" << std::endl;
	}
	else if (argc == 2)
	{
		std::cout << "Programma di conversione da output binario Fortran ad output \"umano\"" << std::endl;
		std::cout << "Che tipo di file e'? 1 campi, 2 particelle: ";
		std::cin >> funzione;
		std::cout << "E' necessario fare lo swap dell'endian? 1 si', 0 no: ";
		std::cin >> out_swap;
		if (funzione == 2)
		{
			std::cout << "Inserisci il WEIGHT: 0 per output vecchio, 1 per nuovo: ";
			std::cin >> WEIGHT;
			std::cout << "Vuoi l'output binario pulito dal fortran? 1 si', 0 no: ";
			std::cin >> out_binary;
			std::cout << "Vuoi l'output dello spazio delle fasi ascii? 1 si', 0 no: ";
			std::cin >> out_ascii;
		}
	}
	else if(argc == 7)
	{
		funzione = std::atoi(argv[2]);
		out_swap = std::atoi(argv[3]);
		if (funzione == 2)
		{
			WEIGHT = std::atoi(argv[4]);
			out_binary = std::atoi(argv[5]);
			out_ascii = std::atoi(argv[6]);
		}
		std::cout << "Tipo file (1 campi, 2 particelle): " << funzione << std::endl;
		std::cout << "Swapping ordine byte (1 si', 0 no): " << out_swap << std::endl;
		if (funzione == 2)
		{
			std::cout << "WEIGHT (0 per output vecchio, 1 per nuovo): " << WEIGHT << std::endl;
			std::cout << "Output binario pulito dal fortran (1 si', 0 no): " << out_binary << std::endl;
			std::cout << "Output dello spazio delle fasi ascii (1 si', 0 no): " << out_ascii << std::endl;
		}
	}
	else
	{
		std::cout << "Linea di comando non valida" << std::endl;
		return -1;
	}

	if ((funzione != 1 && funzione != 2) || (out_swap != 0 && out_swap != 1) || (funzione == 2 && WEIGHT != 0 && WEIGHT != 1)  
		|| (funzione==2 && out_ascii != 0 && out_ascii !=1) || (funzione==2 && out_binary != 0 && out_binary != 1))
	{
		std::cout << "Input sbagliato!" << std::endl;
		return -2;
	}
	
	if (funzione == 1) leggi_campi(argv[1], out_swap);
	else if (funzione == 2) leggi_particelle(argv[1], WEIGHT, out_swap, out_binary, out_ascii);
	return 0;
}
