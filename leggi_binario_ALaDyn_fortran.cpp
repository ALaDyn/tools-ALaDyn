

#include "leggi_binario_ALaDyn_fortran.h"



int main (const int argc, const char *argv[])
{
	parametri parametri;
	bool testParametri = true;;
	bool fallita_lettura_inputfile = true;

	std::cout << "Binary file reader v2.0" << std::endl;

	if (argc == 1)
	{
		std::cout << "Si usa:\n-manuale: ./reader Nomefile (senza estensione!)" << std::endl;
		std::cout << "-batch: ./reader Nomefile  --parameters (see source code for details, for now)" << std::endl;
		return -1;
	}

	std::ostringstream nomefile_bin, nomefile_dat;
	nomefile_bin << std::string(argv[1]) << ".bin";
	nomefile_dat << std::string(argv[1]) << ".dat";
	std::string riga_persa, endianness, columns;
	std::ifstream file_dat, file_bin;
	file_dat.open(nomefile_dat.str().c_str());
	if (file_dat.fail()) parametri.old_fortran_bin = true;
	else
	{
		parametri.old_fortran_bin = false;
		file_dat >> endianness >> columns;
		std::getline(file_dat,riga_persa); // nb: serve per pulire dallo stdin gli spazi residui fino al newline
		if (endianness == "BIG-ENDIAN") parametri.endian_file = 1;
		else if (endianness == "LITTLE-ENDIAN") parametri.endian_file = 0;
		if (parametri.endian_file == parametri.endian_machine) 
		{
			parametri.p[SWAP] = 0;
			parametri.p_b[SWAP] = false;
		}
		else 
		{
			parametri.p[SWAP] = 1;
			parametri.p_b[SWAP] = false;
		}
		parametri.p[NCOLUMNS] = std::atoi(columns.c_str());
		parametri.p_b[NCOLUMNS] = false;
	}
	file_bin.open(nomefile_bin.str().c_str());
	fallita_lettura_inputfile = file_bin.fail();
	file_bin.close();
	file_dat.close();


	if (argc == 2)
	{
		std::cout << "Programma di conversione da output binario Fortran ad output \"umano\"" << std::endl;
		parametri.leggi_interattivo();
	}

	else if (std::string(argv[2]) == "-readParamfromFile")
	{
		parametri.leggi_da_file(argv[3]);
		if (argc > 3) 
		{
			parametri.leggi_da_shell(argc,argv);
		}
		if (parametri.incompleto()) parametri.leggi_interattivo();
	}

	else
	{
		parametri.leggi_da_shell(argc, argv);
		if (parametri.incompleto()) parametri.leggi_interattivo();
	}

#ifdef ENABLE_DEBUG
	for (int i = 0; i < NPARAMETRI; i++) std::cout << "p[" << i << "] = " << parametri.p[i] << std::endl;
#endif


	testParametri = parametri.check_parametri();

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
	if ( parametri.p[FUNZIONE] == 2 && parametri.p[FIND_MINMAX] && parametri.p[DO_BINNING] )
	{
		std::cout << "Non riesco a fare altro che cercare minimi e massimi" << std::endl;
		return -4;
	}

	if (parametri.p[FUNZIONE] == 1) leggi_campi(argc, argv, &parametri);
	else if (parametri.p[FUNZIONE] == 2) leggi_particelle(argc, argv, &parametri);
	return 0;
}


