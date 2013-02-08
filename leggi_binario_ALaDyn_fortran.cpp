

#include "leggi_binario_ALaDyn_fortran.h"



int main (int argc, char *argv[])
{
	parametri binning;
	bool testParametri;
	bool fallita_lettura_inputfile = true;
	std::ifstream inputfile;

	if (argc == 1)
	{
		std::cout << "Si usa:\n-manuale: ./reader Nomefile.bin" << std::endl;
		std::cout << "-batch: ./reader Nomefile.bin --parameters (see source code for details, for now)" << std::endl;
	}
	else if (argc == 2)
	{
		std::cout << "Programma di conversione da output binario Fortran ad output \"umano\"" << std::endl;
		binning.leggi_interattivo();
	}

	else if (std::string(argv[2]) == "-readParamfromFile")
	{
		binning.leggi_da_file(argv[3]);
	}

	else
	{
		binning.leggi_da_shell(argc, argv);
	}

	inputfile.open(argv[1]);
	fallita_lettura_inputfile = inputfile.fail();
	inputfile.close();
	testParametri = binning.check_parametri();

	if (
		    ( binning.p[FUNZIONE] != 1 && binning.p[FUNZIONE]   != 2								) 
		 || ( binning.p[SWAP]     != 0 && binning.p[SWAP]       != 1								) 
		 || ( binning.p[FUNZIONE] == 2 && binning.p[WEIGHT]     != 0 && binning.p[WEIGHT]     != 1	)  
		 || ( binning.p[FUNZIONE] == 2 && binning.p[OUT_ASCII]  != 0 && binning.p[OUT_ASCII]  != 1	) 
		 || ( binning.p[FUNZIONE] == 2 && binning.p[OUT_BINARY] != 0 && binning.p[OUT_BINARY] != 1  ) 
		 || ( testParametri == false																) 
		 || ( fallita_lettura_inputfile																)
		 || ( binning.p[FIND_MINMAX] && binning.p[DO_BINNING]										)
	   )
	{
		std::cout << "Input sbagliato! (maggiori dettagli in futuro...)" << std::endl;
		return -2;
	}

	if (binning.p[FUNZIONE] == 1) leggi_campi(argv[1], binning);
	else if (binning.p[FUNZIONE] == 2) leggi_particelle(argv[1], binning);
	return 0;
}


