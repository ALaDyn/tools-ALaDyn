#ifndef __BINNING_C
#define __BINNING_C

#include "leggi_binario_ALaDyn_fortran.h"



/* costruttore default - inizializza a zero per ora. Siccome tutto e' inizializzato a zero 
non puo' nemmeno calcolare le le dimensioni dei bin!
Verificare che non ci siano cose intelligenti da poter fare! */
parametri :: parametri()
{
	nbin_x = nbin_px = nbin_y = nbin_z = nbin_py = nbin_pz = nbin_E = nbin_theta = 0;
	xmin = xmax = pxmin = pxmax = ymin = ymax = pymin = pymax = zmin = zmax = pzmin = pzmax = thetamin = thetamax = Emin = Emax = gammamin = gammamax = 0.0;
	//		dim_x = dim_y = dim_z = dim_px = dim_py = dim_pz = dim_E = dim_theta = 0.0;
	xmin_b = xmax_b = pxmin_b = pxmax_b = ymin_b = ymax_b = pymin_b = pymax_b = zmin_b = zmax_b = pzmin_b = pzmax_b = Emin_b = Emax_b = gammamin_b = gammamax_b = thetamin_b = thetamax_b = true;
	nbin_x_b = nbin_px_b = nbin_y_b = nbin_z_b = nbin_py_b = nbin_pz_b = nbin_E_b = nbin_theta_b = true;
	for (int i = 0; i < NPARAMETRI; i++) 
	{
		p_b[i] = true;
		p[i] = -1;
	}
	old_fortran_bin = false;
	endian_file = 0;
	endian_machine = is_big_endian();
}

/* costruttore parametrico 1D */
parametri :: parametri(float xm, float xM, float pxm, float pxM, int nx, int npx)
{
//	non implementato
}

float parametri :: dimmi_dimx()
{
	return (xmax - xmin) / static_cast <float> (nbin_x);
}
float parametri :: dimmi_dimy()
{
	return (ymax - ymin) / static_cast <float> (nbin_y);
}
float parametri :: dimmi_dimz()
{
	return (zmax - zmin) / static_cast <float> (nbin_z);
}
float parametri :: dimmi_dimpx()
{
	return (pxmax - pxmin) / static_cast <float> (nbin_px);
}
float parametri :: dimmi_dimpy()
{
	return (pymax - pymin) / static_cast <float> (nbin_py);
}
float parametri :: dimmi_dimpz()
{
	return (pzmax - pzmin) / static_cast <float> (nbin_pz);
}
float parametri :: dimmi_dimgamma()
{
	return (gammamax - gammamin) / static_cast <float> (nbin_gamma);
}
float parametri :: dimmi_dimtheta()
{
	return (thetamax - thetamin) / static_cast <float> (nbin_theta);
}
float parametri :: dimmi_dimE()
{
	return (Emax - Emin) / static_cast <float> (nbin_E);
}
float parametri :: dimmi_dim(int colonna)
{
	if (colonna == 0)		return dimmi_dimx();
	else if (colonna == 1)	return dimmi_dimy();
	else if (colonna == 2)	return dimmi_dimz();
	else if (colonna == 3)	return dimmi_dimpx();
	else if (colonna == 4)	return dimmi_dimpy();
	else if (colonna == 5)	return dimmi_dimpz();
	else if (colonna == 6)	return dimmi_dimgamma();
	else if (colonna == 7)	return dimmi_dimtheta();
	else if (colonna == 8)	return dimmi_dimE();
	else return 1.0;
}


void parametri :: leggi_da_file(const char *nomefile)
{
	bool failed_opening_file = true;
	std::ifstream fileParametri;
	fileParametri.open(nomefile);
	failed_opening_file = fileParametri.fail();
	if (failed_opening_file) 
	{
		std::cout << "Impossibile aprire il file dei minimi/massimi" << std::endl;
		return;
	}

	std::string nomepar, evita, leggi;

	if (p_b[FUNZIONE])
	{
		std::cout << "Che tipo di file e'? 1 campi, 2 particelle: ";
		std::cin >> p[FUNZIONE];
		p_b[FUNZIONE] = false;
	}
	if (p_b[SWAP])
	{
		std::cout << "E' necessario fare lo swap dell'endian? 1 si', 0 no: ";
		std::cin >> p[SWAP];
		p_b[SWAP] = false;
	}
	if (p[FUNZIONE] == 2)
	{
		if (p_b[WEIGHT])
		{
			std::cout << "C'e' la colonna dei weight? 0 per output vecchio (no), 1 per nuovo (si'): ";
			std::cin >> p[WEIGHT];
			p_b[WEIGHT] = false;
		}
		if (p_b[FIND_MINMAX])
		{
			std::cout << "Vuoi cercare massimi e minimi? 0 per no, 1 per si': ";
			std::cin >> p[FIND_MINMAX];
			p_b[FIND_MINMAX] = false;
		}
		if (p[FIND_MINMAX] != 1)
		{
			if (p_b[DO_BINNING])
			{
				std::cout << "Vuoi fare il binnaggio dei dati? 1 si', 0 no: ";
				std::cin >> p[DO_BINNING];
				p_b[DO_BINNING] = false;
			}
			if (p[DO_BINNING] == 1)
			{
				std::cout << "Quanti bin per asse vuoi usare? (nb: Per ora sono tutti uguali): ";
				std::cin >> nbin_x;
				nbin_px = nbin_y = nbin_z = nbin_py = nbin_pz = nbin_E = nbin_theta = nbin_x;
				while(!fileParametri.eof())
				{
					fileParametri >> nomepar >> evita >> leggi;
					if ((nomepar == "xmin" || nomepar == "XMIN") && xmin_b)
					{
						xmin = (float) std::atof(leggi.c_str());
						xmin_b = false;
					}
					else if ((nomepar == "xmax" || nomepar == "XMAX") && xmax_b)
					{
						xmax = (float) std::atof(leggi.c_str());
						xmax_b = false;
					}
					else if ((nomepar == "ymin" || nomepar == "YMIN") && ymin_b)
					{
						ymin = (float) std::atof(leggi.c_str());
						ymin_b = false;
					}
					else if ((nomepar == "ymax" || nomepar == "YMAX") && ymax_b)
					{
						ymax = (float) std::atof(leggi.c_str());
						ymax_b = false;
					}
					else if ((nomepar == "zmin" || nomepar == "ZMIN") && zmin_b)
					{
						zmin = (float) std::atof(leggi.c_str());
						zmin_b = false;
					}
					else if ((nomepar == "zmax" || nomepar == "ZMAX") && zmax_b)
					{
						zmax = (float) std::atof(leggi.c_str());
						zmax_b = false;
					}
					else if ((nomepar == "pxmin" || nomepar == "PXMIN") && pxmin_b)
					{
						pxmin = (float) std::atof(leggi.c_str());
						pxmin_b = false;
					}
					else if ((nomepar == "pxmax" || nomepar == "PXMAX") && pxmax_b)
					{
						pxmax = (float) std::atof(leggi.c_str());
						pxmax_b = false;
					}
					else if ((nomepar == "pymin" || nomepar == "PYMIN") && pymin_b)
					{
						pymin = (float) std::atof(leggi.c_str());
						pymin_b = false;
					}
					else if ((nomepar == "pymax" || nomepar == "PYMAX") && pymax_b)
					{
						pymax = (float) std::atof(leggi.c_str());
						pymax_b = false;
					}
					else if ((nomepar == "pzmin" || nomepar == "PZMIN") && pzmin_b)
					{
						pzmin = (float) std::atof(leggi.c_str());
						pzmin_b = false;
					}
					else if ((nomepar == "pzmax" || nomepar == "PZMAX") && pzmax_b)
					{
						pzmax = (float) std::atof(leggi.c_str());
						pzmax_b = false;
					}
					else if ((nomepar == "gammamin" || nomepar == "GAMMAMIN") && gammamin_b)
					{
						gammamin = (float) std::atof(leggi.c_str());
						gammamin_b = false;
					}
					else if ((nomepar == "gammamax" || nomepar == "GAMMAMAX") && gammamax_b)
					{
						gammamax = (float) std::atof(leggi.c_str());
						gammamax_b = false;
					}
					else if ((nomepar == "thetamin" || nomepar == "THETAMIN") && thetamin_b)
					{
						thetamin = (float) std::atof(leggi.c_str());
						thetamin_b = false;
					}
					else if ((nomepar == "thetamax" || nomepar == "THETAMAX") && thetamax_b)
					{
						thetamax = (float) std::atof(leggi.c_str());
						thetamax_b = false;
					}
					else if ((nomepar == "emin" || nomepar == "EMIN") && Emin_b)
					{
						Emin = (float) std::atof(leggi.c_str());
						Emin_b = false;
					}
					else if ((nomepar == "emax" || nomepar == "EMAX") && Emax_b)
					{
						Emax = (float) std::atof(leggi.c_str());
						Emax_b = false;
					}
					else
					{
						std::cout << "Parametro " << nomepar << " non riconosciuto." << std::endl;
					}
				}
#ifdef ENABLE_DEBUG
				std::cout << "Dal file " << nomefile << " ho letto i seguenti estremi:" << std::endl;
				std::cout << "XMIN = " << xmin << std::endl;
				std::cout << "XMAX = " << xmax << std::endl;
				std::cout << "YMIN = " << ymin << std::endl;
				std::cout << "YMAX = " << ymax << std::endl;
				std::cout << "ZMIN = " << zmin << std::endl;
				std::cout << "ZMAX = " << zmax << std::endl;
				std::cout << "PXMIN = " << pxmin << std::endl;
				std::cout << "PXMAX = " << pxmax << std::endl;
				std::cout << "PYMIN = " << pymin << std::endl;
				std::cout << "PYMAX = " << pymax << std::endl;
				std::cout << "PZMIN = " << pzmin << std::endl;
				std::cout << "PZMAX = " << pzmax << std::endl;
				std::cout << "GAMMAMIN = " << gammamin << std::endl;
				std::cout << "GAMMAMAX = " << gammamax << std::endl;
				std::cout << "THETAMIN = " << thetamin << std::endl;
				std::cout << "THETAMAX = " << thetamax << std::endl;
				std::cout << "EMIN = " << Emin << std::endl;
				std::cout << "EMAX = " << Emax << std::endl;
#endif
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
	fileParametri.close();
}


void parametri :: leggi_interattivo()
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
				if (xmin_b)
				{
					std::cout << "xmin = ";
					std::cin >> xmin;
					xmin_b = false;
				}
				if (xmax_b)
				{
					std::cout << "xmax = ";
					std::cin >> xmax;
					xmax_b = false;
				}
				if (nbin_x_b)
				{
					std::cout << "nbin_x = ";
					std::cin >> nbin_x;
					nbin_x_b = false;
				}
				if (ymin_b)
				{
					std::cout << "ymin = ";
					std::cin >> ymin;
					ymin_b = false;
				}
				if (ymax_b)
				{
					std::cout << "ymax = ";
					std::cin >> ymax;
					ymax_b = false;
				}
				if (nbin_y_b)
				{
					std::cout << "nbin_y = ";
					std::cin >> nbin_y;
					nbin_y_b = false;
				}
				if (zmin_b)
				{
					std::cout << "zmin = ";
					std::cin >> zmin;
					zmin_b = false;
				}
				if (zmax_b)
				{
					std::cout << "zmax = ";
					std::cin >> zmax;
					zmax_b = false;
				}
				if (nbin_z_b)
				{
					std::cout << "nbin_z = ";
					std::cin >> nbin_z;
					nbin_z_b = false;
				}
				if (pxmin_b)
				{
					std::cout << "pxmin = ";
					std::cin >> pxmin;
					pxmin_b = false;
				}
				if (pxmax_b)
				{
					std::cout << "pxmax = ";
					std::cin >> pxmax;
					pxmax_b = false;
				}
				if (nbin_px_b)
				{
					std::cout << "nbin_px = ";
					std::cin >> nbin_px;
					nbin_px_b = false;
				}
				if (pymin_b)
				{
					std::cout << "pymin = ";
					std::cin >> pymin;
					pymin_b = false;
				}
				if (pymax_b)
				{
					std::cout << "pymax = ";
					std::cin >> pymax;
					pymax_b = false;
				}
				if (nbin_py_b)
				{
					std::cout << "nbin_py = ";
					std::cin >> nbin_py;
					nbin_py_b = false;
				}
				if (pzmin_b)
				{
					std::cout << "pzmin = ";
					std::cin >> pzmin;
					pzmin_b = false;
				}
				if (pzmax_b)
				{
					std::cout << "pzmax = ";
					std::cin >> pzmax;
					pzmax_b = false;
				}
				if (nbin_pz_b)
				{
					std::cout << "nbin_pz = ";
					std::cin >> nbin_pz;
					nbin_pz_b = false;
				}
				if (gammamin_b)
				{
					std::cout << "gammamin = ";
					std::cin >> gammamin;
					gammamin_b = false;
				}
				if (gammamax_b)
				{
					std::cout << "gammamax = ";
					std::cin >> gammamax;
					gammamax_b = false;
				}
				if (nbin_gamma_b)
				{
					std::cout << "nbin_gamma = ";
					std::cin >> nbin_gamma;
					nbin_gamma_b = false;
				}
				if (thetamin_b)
				{
					std::cout << "thetamin = ";
					std::cin >> thetamin;
					thetamin_b = false;
				}
				if (thetamax_b)
				{
					std::cout << "thetamax = ";
					std::cin >> thetamax;
					thetamax_b = false;
				}
				if (nbin_theta_b)
				{
					std::cout << "nbin_theta = ";
					std::cin >> nbin_theta;
					nbin_theta_b = false;
				}
				if (Emin_b)
				{
					std::cout << "Emin = ";
					std::cin >> Emin;
					Emin_b = false;
				}
				if (Emax_b)
				{
					std::cout << "Emax = ";
					std::cin >> Emax;
					Emax_b = false;
				}
				if (nbin_E_b)
				{
					std::cout << "nbin_E = ";
					std::cin >> nbin_E;
					nbin_E_b = false;
				}
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

void parametri :: leggi_da_shell(int argc, const char *argv[])
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
			xmin_b = false;
			i++;
		}
		else if (std::string(argv[i]) == "-xmax")
		{
			xmax = (float) atof(argv[i+1]);
			xmax_b = false;
			i++;
		}
		else if (std::string(argv[i]) == "-ymin")
		{
			ymin = (float) atof(argv[i+1]);
			ymin_b = false;
			i++;
		}
		else if (std::string(argv[i]) == "-ymax")
		{
			ymax = (float) atof(argv[i+1]);
			ymax_b = false;
			i++;
		}
		else if (std::string(argv[i]) == "-zmin")
		{
			zmin = (float) atof(argv[i+1]);
			zmin_b = false;
			i++;
		}
		else if (std::string(argv[i]) == "-zmax")
		{
			zmax = (float) atof(argv[i+1]);
			zmax_b = false;
			i++;
		}
		else if (std::string(argv[i]) == "-pxmin")
		{
			pxmin = (float) atof(argv[i+1]);
			pxmin_b = false;
			i++;
		}
		else if (std::string(argv[i]) == "-pxmax")
		{
			pxmax = (float) atof(argv[i+1]);
			pxmax_b = false;
			i++;
		}
		else if (std::string(argv[i]) == "-pymin")
		{
			pymin = (float) atof(argv[i+1]);
			pymin_b = false;
			i++;
		}
		else if (std::string(argv[i]) == "-pymax")
		{
			pymax = (float) atof(argv[i+1]);
			pymax_b = false;
			i++;
		}
		else if (std::string(argv[i]) == "-pzmin")
		{
			pzmin = (float) atof(argv[i+1]);
			pzmin_b = false;
			i++;
		}
		else if (std::string(argv[i]) == "-pzmax")
		{
			pzmax = (float) atof(argv[i+1]);
			pzmax_b = false;
			i++;
		}
		else if (std::string(argv[i]) == "-thetamin")
		{
			thetamin = (float) atof(argv[i+1]);
			thetamin_b = false;
			i++;
		}
		else if (std::string(argv[i]) == "-thetamax")
		{
			thetamax = (float) atof(argv[i+1]);
			thetamax_b = false;
			i++;
		}
		else if (std::string(argv[i]) == "-gammamin")
		{
			gammamin = (float) atof(argv[i+1]);
			gammamin_b = false;
			i++;
		}
		else if (std::string(argv[i]) == "-gammamax")
		{
			gammamax = (float) atof(argv[i+1]);
			gammamax_b = false;
			i++;
		}
		else if (std::string(argv[i]) == "-Emin")
		{
			Emin = (float) atof(argv[i+1]);
			Emin_b = false;
			i++;
		}
		else if (std::string(argv[i]) == "-Emax")
		{
			Emax = (float) atof(argv[i+1]);
			Emax_b = false;
			i++;
		}
		else if (std::string(argv[i]) == "-nbinx")
		{
			nbin_x = atoi(argv[i+1]);
			nbin_x_b = false;
			i++;
		}
		else if (std::string(argv[i]) == "-nbiny")
		{
			nbin_y = atoi(argv[i+1]);
			nbin_y_b = false;
			i++;
		}
		else if (std::string(argv[i]) == "-nbinz")
		{
			nbin_z = atoi(argv[i+1]);
			nbin_z_b = false;
			i++;
		}
		else if (std::string(argv[i]) == "-nbinpx")
		{
			nbin_px = atoi(argv[i+1]);
			nbin_px_b = false;
			i++;
		}
		else if (std::string(argv[i]) == "-nbinpy")
		{
			nbin_py = atoi(argv[i+1]);
			nbin_py_b = false;
			i++;
		}
		else if (std::string(argv[i]) == "-nbinpz")
		{
			nbin_pz = atoi(argv[i+1]);
			nbin_pz_b = false;
			i++;
		}
		else if (std::string(argv[i]) == "-nbintheta")
		{
			nbin_gamma = atoi(argv[i+1]);
			nbin_gamma_b = false;
			i++;
		}
		else if (std::string(argv[i]) == "-nbingamma")
		{
			nbin_theta = atoi(argv[i+1]);
			nbin_theta_b = false;
			i++;
		}
		else if (std::string(argv[i]) == "-nbinE")
		{
			nbin_E = atoi(argv[i+1]);
			nbin_E_b = false;
			i++;
		}
	}
}



bool parametri :: check_parametri()
{
	bool test = true;
	//	( p[FUNZIONE] != 1 && p[FUNZIONE]   != 2						) 
	//	( p[SWAP]     != 0 && p[SWAP]       != 1						) 
	//	( p[FUNZIONE] == 2 && p[WEIGHT]     != 0 && p[WEIGHT]     != 1	)  
	//	( p[FUNZIONE] == 2 && p[OUT_ASCII]  != 0 && p[OUT_ASCII]  != 1	) 
	//	( p[FUNZIONE] == 2 && p[OUT_BINARY] != 0 && p[OUT_BINARY] != 1  ) 
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


void parametri :: organizza_minimi_massimi()
{
	minimi[0] = xmin;
	minimi[1] = ymin;
	minimi[2] = zmin;
	minimi[3] = pxmin;
	minimi[4] = pymin;
	minimi[5] = pzmin;
	minimi[6] = gammamin;
	minimi[7] = thetamin;
	minimi[8] = Emin;

	massimi[0] = xmin;
	massimi[1] = ymin;
	massimi[2] = zmin;
	massimi[3] = pxmin;
	massimi[4] = pymin;
	massimi[5] = pzmin;
	massimi[6] = gammamin;
	massimi[7] = thetamin;
	massimi[8] = Emin;
}

bool parametri :: incompleto()
{
	bool test=true;
	if (!xmin_b && !xmax_b && !pxmin_b && !pxmax_b && !ymin_b && !ymax_b && !pymin_b && !pymax_b && !zmin_b && !zmax_b && !pzmin_b && !pzmax_b && 
		!Emin_b && !Emax_b && !nbin_x_b && !nbin_y_b && !nbin_z_b && !nbin_px_b && !nbin_py_b && !nbin_pz_b && !nbin_E_b)  test = false;
	//	if (!gammamin_b && !gammamax_b && !thetamin_b && !thetamax_b && !nbin_theta_b && !nbin_gamma_b && !test) test = false;
	return test;
}


#endif
