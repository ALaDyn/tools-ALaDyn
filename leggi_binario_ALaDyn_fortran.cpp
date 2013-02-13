

#include "leggi_binario_ALaDyn_fortran.h"



int main (int argc, char *argv[])
{
	parametri binning;
	bool testParametri = true;;
	bool fallita_lettura_inputfile = true;
	std::ifstream inputfile;

	std::cout << "Binary file reader v2.0" << std::endl;

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

#ifdef ENABLE_DEBUG
	for (int i = 0; i < NPARAMETRI; i++) std::cout << "p[" << i << "] = " << binning.p[i] << std::endl;
#endif

	inputfile.open(argv[1]);
	fallita_lettura_inputfile = inputfile.fail();
	inputfile.close();
	testParametri = binning.check_parametri();

	if ( testParametri == false	)
	{
		std::cout << "Parametri non coerenti" << std::endl;
		return -2;
	}
	if ( fallita_lettura_inputfile )
	{
		std::cout << "Input file non trovato" << std::endl;
		return -3;
	}
	if ( binning.p[FUNZIONE] == 2 && binning.p[FIND_MINMAX] && binning.p[DO_BINNING] )
	{
		std::cout << "Non riesco a fare altro che cercare minimi e massimi" << std::endl;
		return -4;
	}

	if (binning.p[FUNZIONE] == 1) leggi_campi(argv[1], binning);
	else if (binning.p[FUNZIONE] == 2) leggi_particelle(argv[1], binning);
	return 0;
}


