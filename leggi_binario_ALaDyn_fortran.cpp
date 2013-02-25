

#include "leggi_binario_ALaDyn_fortran.h"



int main (const int argc, const char *argv[])
{
	parametri parametri;
	bool testParametri = true;;
	bool fallita_lettura_inputfile = true;

	std::cout << "Binary file reader v3.0" << std::endl;

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
	file_bin.open(nomefile_bin.str().c_str());
	fallita_lettura_inputfile = file_bin.fail();
	file_bin.close();
	if ( fallita_lettura_inputfile )
	{
		std::cout << "Input file non trovato" << std::endl;
		return -3;
	}
	else
	{
		parametri.check_filename(argv[1]);
	}
	file_dat.open(nomefile_dat.str().c_str());
	if (file_dat.fail()) 
	{
		parametri.old_fortran_bin = true;
		parametri.chiedi_endian_file();
		parametri.chiedi_numero_colonne();
	}
	else
	{
		parametri.old_fortran_bin = false;
		file_dat >> endianness >> columns;
		std::getline(file_dat,riga_persa); // nb: serve per pulire dallo stdin gli spazi residui fino al newline
		if (endianness == "BIG-ENDIAN")
		{
			parametri.endian_file = 1;
			parametri.p[NCOLUMNS] = std::atoi(columns.c_str());
			parametri.p_b[NCOLUMNS] = false;
		}
		else if (endianness == "LITTLE-ENDIAN")
		{
			parametri.endian_file = 0;
			parametri.p[NCOLUMNS] = std::atoi(columns.c_str());
			parametri.p_b[NCOLUMNS] = false;
		}
		else
		{
			parametri.chiedi_endian_file();
			parametri.chiedi_numero_colonne();
		}
	}
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
	file_dat.close();

	parametri.check_filename(argv[1]);

	if (argc == 2)
	{
		parametri.leggi_interattivo();
	}

	else if (std::string(argv[2]) == "-readParamsfromFile")
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

	/*
	if ( parametri.p[FUNZIONE] == 2 && parametri.p[FIND_MINMAX] && parametri.p[DO_BINNING] )
	{
	std::cout << "Non riesco a fare altro che cercare minimi e massimi" << std::endl;
	return -4;
	}
	*/

	if (parametri.p[DO_BINNING]) parametri.organizza_minimi_massimi();

	if (parametri.p[FUNZIONE] == 1  || parametri.file_campi_Ex || parametri.file_campi_Ey || parametri.file_campi_Ez 
									|| parametri.file_campi_Bx || parametri.file_campi_By || parametri.file_campi_Bz 
									|| parametri.file_densita_elettroni || parametri.file_densita_protoni 
									|| parametri.file_densita_HI || parametri.file_densita_LI) leggi_campi(argc, argv, &parametri);
	else if (parametri.p[FUNZIONE] == 2 || parametri.file_particelle_P || parametri.file_particelle_E 
										|| parametri.file_particelle_HI || parametri.file_particelle_LI) leggi_particelle(argc, argv, &parametri);
	return 0;
}


