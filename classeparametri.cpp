#ifndef __BINNING_C
#define __BINNING_C

#include "leggi_binario_ALaDyn_fortran.h"


class parametri
{
public:
	int nbin_x, nbin_y, nbin_z, nbin_px, nbin_py, nbin_pz, nbin_E, nbin_theta;
	int p[NPARAMETRI];
	char support_label[MAX_LENGTH_FILENAME];
	float xmin, xmax, pxmin, pxmax, ymin, ymax, pymin, pymax, zmin, zmax, pzmin, pzmax, Emin, Emax, thetamin, thetamax;
//	float dim_x, dim_y, dim_z, dim_px, dim_py, dim_pz, dim_E, dim_theta;

	/* costruttore default - inizializza a zero per ora. Siccome tutto e' inizializzato a zero 
	non puo' nemmeno calcolare le le dimensioni dei bin!
	Verificare che non ci siano cose intelligenti da poter fare! */
	parametri()
	{
		nbin_x = nbin_px = nbin_y = nbin_z = nbin_py = nbin_pz = nbin_E = nbin_theta = 0;
		xmin = xmax = pxmin = pxmax = ymin = ymax = pymin = pymax = zmin = zmax = pzmin = pzmax = thetamin = thetamax = 0.0;
//		dim_x = dim_y = dim_z = dim_px = dim_py = dim_pz = dim_E = dim_theta = 0.0;
	}

	/* costruttore parametrico 1D */
	parametri(float xm, float xM, float pxm, float pxM, int nx, int npx)
	{
		nbin_x = nx;
		nbin_px = npx;
		nbin_y = nbin_z = nbin_py = nbin_pz = nbin_E = nbin_theta = 0;
		xmin = xm;
		xmax = xM;
		pxmin = pxm;
		pxmax = pxM;
		ymin = ymax = pymin = pymax = zmin = zmax = pzmin = pzmax = thetamin = thetamax = 0.0;
		Emin = Emax = 0.0;	// E non calcolata in 1D
		nbin_E = 0;
//		dim_x = (xmax - xmin) / nbin_x;
//		dim_px = (pxmax - pxmin) / nbin_px;
//		dim_y = dim_z = dim_py = dim_pz = dim_E = dim_theta = 0.0;
	}

	float dimmi_dimx()
	{
		return (xmax - xmin) / nbin_x;
	}
	float dimmi_dimy()
	{
		return (ymax - ymin) / nbin_y;
	}
	float dimmi_dimz()
	{
		return (zmax - zmin) / nbin_z;
	}
	float dimmi_dimpx()
	{
		return (pxmax - pxmin) / nbin_px;
	}
	float dimmi_dimpy()
	{
		return (pymax - pymin) / nbin_py;
	}
	float dimmi_dimpz()
	{
		return (pzmax - pzmin) / nbin_pz;
	}
	float dimmi_dimtheta()
	{
		return (thetamax - thetamin) / nbin_theta;
	}
	float dimmi_dimE()
	{
		return (Emax - Emin) / nbin_E;
	}

	void leggi_da_file(const char *nomefile)
	{
		std::ifstream fileParametri;
		fileParametri.open(nomefile);
		std::string evita, leggi;

		std::cout << "Che tipo di file e'? 1 campi, 2 particelle: ";
		std::cin >> p[FUNZIONE];
		std::cout << "E' necessario fare lo swap dell'endian? 1 si', 0 no: ";
		std::cin >> p[SWAP];
		if (p[FUNZIONE] == 2)
		{
			std::cout << "C'e' la colonna dei weight? 0 per output vecchio (no), 1 per nuovo (si'): ";
			std::cin >> p[WEIGHT];
			std::cout << "Vuoi cercare massimi e minimi? 0 per no, 1 per si': ";
			std::cin >> p[FIND_MINMAX];
			if (p[FIND_MINMAX] != 1)
			{
				std::cout << "Vuoi fare il binnaggio dei dati? 1 si', 0 no: ";
				std::cin >> p[DO_BINNING];
				if (p[DO_BINNING] == 1)
				{
					fileParametri >> evita >> evita;
					fileParametri >> leggi; 
					xmin = (float) atof(leggi.c_str());
					fileParametri >> evita >> evita;
					fileParametri >> leggi;
					xmax = (float) atof(leggi.c_str());
					fileParametri >> evita >> evita;
					fileParametri >> leggi;
					ymin = (float) atof(leggi.c_str());
					fileParametri >> evita >> evita;
					fileParametri >> leggi;
					ymax = (float) atof(leggi.c_str());
					fileParametri >> evita >> evita;
					fileParametri >> leggi;
					zmin = (float) atof(leggi.c_str());
					fileParametri >> evita >> evita;
					fileParametri >> leggi;
					zmax = (float) atof(leggi.c_str());
					fileParametri >> evita >> evita;
					fileParametri >> leggi;
					pxmin = (float) atof(leggi.c_str());
					fileParametri >> evita >> evita;
					fileParametri >> leggi;
					pxmax = (float) atof(leggi.c_str());
					fileParametri >> evita >> evita;
					fileParametri >> leggi;
					pymin = (float) atof(leggi.c_str());
					fileParametri >> evita >> evita;
					fileParametri >> leggi;
					pymax = (float) atof(leggi.c_str());
					fileParametri >> evita >> evita;
					fileParametri >> pzmin;
					pzmin = (float) atof(leggi.c_str());
					fileParametri >> evita >> evita;
					fileParametri >> leggi;
					pzmax = (float) atof(leggi.c_str());
					fileParametri >> evita >> evita;
					fileParametri >> leggi;
					Emin = (float) atof(leggi.c_str());
					fileParametri >> evita >> evita;
					fileParametri >> leggi;
					Emax = (float) atof(leggi.c_str());
					fileParametri >> evita >> evita;
					fileParametri >> leggi;
					thetamin = (float) atof(leggi.c_str());
					fileParametri >> evita >> evita;
					fileParametri >> leggi;
					thetamax = (float) atof(leggi.c_str());
				}
			}
			else p[DO_BINNING] = 0;
			std::cout << "Vuoi l'output completo binario? 1 si', 0 no: ";
			std::cin >> p[OUT_BINARY];
			std::cout << "Vuoi l'output completo ascii? 1 si', 0 no: ";
			std::cin >> p[OUT_ASCII];
		}
		else if(p[FUNZIONE] == 1)
		{
			std::cout << "Inserisci label per il file: (i.e. Ex) ";
			std::cin >> support_label;
		}
	}


	void leggi_interattivo()
	{
		std::cout << "Che tipo di file e'? 1 campi, 2 particelle: ";
		std::cin >> p[FUNZIONE];
		std::cout << "E' necessario fare lo swap dell'endian? 1 si', 0 no: ";
		std::cin >> p[SWAP];
		if (p[FUNZIONE] == 2)
		{
			std::cout << "C'e' la colonna dei weight? 0 per output vecchio (no), 1 per nuovo (si'): ";
			std::cin >> p[WEIGHT];
			std::cout << "Vuoi cercare massimi e minimi? 0 per no, 1 per si': ";
			std::cin >> p[FIND_MINMAX];
			if (p[FIND_MINMAX] != 1)
			{
				std::cout << "Vuoi fare il binnaggio dei dati? 1 si', 0 no: ";
				std::cin >> p[DO_BINNING];
				if (p[DO_BINNING] == 1)
				{
					std::cout << "xmin = ";
					std::cin >> xmin;
					std::cout << "xmax = ";
					std::cin >> xmax;
					std::cout << "ymin = ";
					std::cin >> ymin;
					std::cout << "ymax = ";
					std::cin >> ymax;
					std::cout << "zmin = ";
					std::cin >> zmin;
					std::cout << "zmax = ";
					std::cin >> zmax;
					std::cout << "pxmin = ";
					std::cin >> pxmin;
					std::cout << "pxmax = ";
					std::cin >> pxmax;
					std::cout << "pymin = ";
					std::cin >> pymin;
					std::cout << "pymax = ";
					std::cin >> pymax;
					std::cout << "pzmin = ";
					std::cin >> pzmin;
					std::cout << "pzmax = ";
					std::cin >> pzmax;
					std::cout << "Emin = ";
					std::cin >> Emin;
					std::cout << "Emax = ";
					std::cin >> Emax;
					std::cout << "thetamin = ";
					std::cin >> thetamin;
					std::cout << "thetamax = ";
					std::cin >> thetamax;
				}
			}
			else p[DO_BINNING] = 0;
			std::cout << "Vuoi l'output completo binario? 1 si', 0 no: ";
			std::cin >> p[OUT_BINARY];
			std::cout << "Vuoi l'output completo ascii? 1 si', 0 no: ";
			std::cin >> p[OUT_ASCII];
			std::cout << "Vuoi l'output dei parametri contenuti nel file? 1 si', 0 no: ";
			std::cin >> p[OUT_PARAMS];
		}
		else if(p[FUNZIONE] == 1)
		{
			std::cout << "Inserisci label per il file: (i.e. Ex) ";
			std::cin >> support_label;
		}
	}

	void leggi_da_shell(int argc, char *argv[])
	{
		for (int i = 2; i < argc; i++)	// * We will iterate over argv[] to get the parameters stored inside.
		{								// * Note that we're starting on 1 because we don't need to know the
										// * path of the program, which is stored in argv[0], and the input file,
										// * which is supposed to be given as the first argument and so is in argv[1]
			if (std::string(argv[i]) == "-swap")
			{
				p[SWAP] = 1;
			}
			else if (std::string(argv[i]) == "-field")
			{
				p[FUNZIONE] = 1;
			}
			else if (std::string(argv[i]) == "-particles")
			{
				p[FUNZIONE] = 2;
			}
			else if (std::string(argv[i]) == "-dump_binary")
			{
				p[OUT_BINARY] = 1;
			}
			else if (std::string(argv[i]) == "-dump_ascii")
			{
				p[OUT_ASCII] = 1;
			}
			else if (std::string(argv[i]) == "-parameters")
			{
				p[OUT_PARAMS] = 1;
			}
			else if (std::string(argv[i]) == "-find_minmax")
			{
				p[FIND_MINMAX] = 1;
			}
			else if (std::string(argv[i]) == "-do_binning")
			{
				p[DO_BINNING] = 1;
			}
			else if (std::string(argv[i]) == "-xmin")
			{
				xmin = (float) atof(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-xmax")
			{
				xmax = (float) atof(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-ymin")
			{
				ymin = (float) atof(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-ymax")
			{
				ymax = (float) atof(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-zmin")
			{
				zmin = (float) atof(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-zmax")
			{
				zmax = (float) atof(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-pxmin")
			{
				pxmin = (float) atof(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-pxmax")
			{
				pxmax = (float) atof(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-pymin")
			{
				pymin = (float) atof(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-pymax")
			{
				pymax = (float) atof(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-pzmin")
			{
				pzmin = (float) atof(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-pzmax")
			{
				pzmax = (float) atof(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-nbinx")
			{
				nbin_x = atoi(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-nbiny")
			{
				nbin_y = atoi(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-nbinz")
			{
				nbin_z = atoi(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-nbinpx")
			{
				nbin_px = atoi(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-nbinpy")
			{
				nbin_py = atoi(argv[i+1]);
				i++;
			}
			else if (std::string(argv[i]) == "-nbinpz")
			{
				nbin_pz = atoi(argv[i+1]);
				i++;
			}
		}
	}


	
	bool check_parametri()
	{
		bool test = true;
//		 || ( binning.p[FUNZIONE] != 1 && binning.p[FUNZIONE]   != 2								) 
//		 || ( binning.p[SWAP]     != 0 && binning.p[SWAP]       != 1								) 
//		 || ( binning.p[FUNZIONE] == 2 && binning.p[WEIGHT]     != 0 && binning.p[WEIGHT]     != 1	)  
//		 || ( binning.p[FUNZIONE] == 2 && binning.p[OUT_ASCII]  != 0 && binning.p[OUT_ASCII]  != 1	) 
//		 || ( binning.p[FUNZIONE] == 2 && binning.p[OUT_BINARY] != 0 && binning.p[OUT_BINARY] != 1  ) 
	//	( xmin < xmax ) ? test = true : test = false;
	//	( ymin < ymax ) ? test = true : test = false;
	//	( zmin < zmax ) ? test = true : test = false;
	//	(pxmin < pxmax) ? test = true : test = false;
	//	(pymin < pymax) ? test = true : test = false;
	//	(pzmin < pzmax) ? test = true : test = false;
	//	( nbin_x > 0) ? test = true : test = false;
	//	( nbin_y > 0) ? test = true : test = false;
	//	( nbin_z > 0) ? test = true : test = false;
	//	(nbin_px > 0) ? test = true : test = false;
	//	(nbin_py > 0) ? test = true : test = false;
	//	(nbin_pz > 0) ? test = true : test = false;
		return test;
	}


};


#endif