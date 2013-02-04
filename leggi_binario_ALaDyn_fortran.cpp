

#include "leggi_binario_ALaDyn_fortran.h"



int main (int argc, char *argv[])
{
	int WEIGHT = 0, funzione = -1, out_swap = -1, out_binary = -1, out_ascii = -1;
	parametri_binnaggio binning;
	bool test;

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
			std::cout << "C'e' la colonna dei weight? 0 per output vecchio (no), 1 per nuovo (si'): ";
			std::cin >> WEIGHT;
			std::cout << "Vuoi l'output binario pulito dal fortran? 1 si', 0 no: ";
			std::cin >> out_binary;
			std::cout << "Vuoi l'output dello spazio delle fasi ascii? 1 si', 0 no: ";
			std::cin >> out_ascii;
		}
	}

	else
	{
		for (int i = 1; i < argc; i++)	// * We will iterate over argv[] to get the parameters stored inside.
		{								// * Note that we're starting on 1 because we don't need to know the
										// * path of the program, which is stored in argv[0]
			if (std::string(argv[i]) == "-swap")
			{
				out_swap = 1;
			}
			else if (std::string(argv[i]) == "-field")
			{
				funzione = 1;
			}
			else if (std::string(argv[i]) == "-particles")
			{
				funzione = 2;
			}
			else if (std::string(argv[i]) == "-dump_binary")
			{
				out_binary = 1;
			}
			else if (std::string(argv[i]) == "-dump_ascii")
			{
				out_ascii = 1;
			}
			else if (std::string(argv[i]) == "-xmin")
			{
				binning.xmin = atof(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-xmax")
			{
				binning.xmax = atof(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-ymin")
			{
				binning.ymin = atof(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-ymax")
			{
				binning.ymax = atof(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-zmin")
			{
				binning.zmin = atof(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-zmax")
			{
				binning.zmax = atof(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-pxmin")
			{
				binning.pxmin = atof(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-pxmax")
			{
				binning.pxmax = atof(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-pymin")
			{
				binning.pymin = atof(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-pymax")
			{
				binning.pymax = atof(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-pzmin")
			{
				binning.pzmin = atof(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-pzmax")
			{
				binning.pzmax = atof(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-nbinx")
			{
				binning.nbin_x = atoi(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-nbiny")
			{
				binning.nbin_y = atoi(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-nbinz")
			{
				binning.nbin_z = atoi(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-nbinpx")
			{
				binning.nbin_px = atoi(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-nbinpy")
			{
				binning.nbin_py = atoi(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-nbinpz")
			{
				binning.nbin_pz = atoi(argv[i+1]);
				i++;
			}
		}
	}

	test = check_parametri(binning);

	if ((funzione != 1 && funzione != 2) || (out_swap != 0 && out_swap != 1) || (funzione == 2 && WEIGHT != 0 && WEIGHT != 1)  
		|| (funzione==2 && out_ascii != 0 && out_ascii !=1) || (funzione==2 && out_binary != 0 && out_binary != 1) || (test == false) )
	{
		std::cout << "Input sbagliato!" << std::endl;
		return -2;
	}

	if (funzione == 1) leggi_campi(argv[1], out_swap);
	else if (funzione == 2) leggi_particelle(argv[1], WEIGHT, out_swap, out_binary, out_ascii, binning);
	return 0;
}



bool check_parametri(parametri_binnaggio analizza)
{
	bool test = true;
//	( analizza.xmin < analizza.xmax ) ? test = true : test = false;
//	( analizza.ymin < analizza.ymax ) ? test = true : test = false;
//	( analizza.zmin < analizza.zmax ) ? test = true : test = false;
//	(analizza.pxmin < analizza.pxmax) ? test = true : test = false;
//	(analizza.pymin < analizza.pymax) ? test = true : test = false;
//	(analizza.pzmin < analizza.pzmax) ? test = true : test = false;
//	( analizza.nbin_x > 0) ? test = true : test = false;
//	( analizza.nbin_y > 0) ? test = true : test = false;
//	( analizza.nbin_z > 0) ? test = true : test = false;
//	(analizza.nbin_px > 0) ? test = true : test = false;
//	(analizza.nbin_py > 0) ? test = true : test = false;
//	(analizza.nbin_pz > 0) ? test = true : test = false;
	return test;
}

