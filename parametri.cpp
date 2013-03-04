#ifndef __BINNING_C
#define __BINNING_C

#include "leggi_binario_ALaDyn_fortran.h"



/* costruttore default - inizializza a zero per ora. Siccome tutto e' inizializzato a zero 
non puo' nemmeno calcolare le le dimensioni dei bin!
Verificare che non ci siano cose intelligenti da poter fare! */
parametri :: parametri()
{
	massa_particella_MeV = 0.;
	nbin_x = nbin_px = nbin_y = nbin_z = nbin_py = nbin_pz = nbin_E = nbin_theta = 100;
	xmin = xmax = pxmin = pxmax = ymin = ymax = pymin = pymax = zmin = zmax = pzmin = pzmax = thetamin = thetamax = Emin = Emax = gammamin = gammamax = 0.0;
	//		dim_x = dim_y = dim_z = dim_px = dim_py = dim_pz = dim_E = dim_theta = 0.0;
	ymin_b = ymax_b = pymin_b = pymax_b = zmin_b = zmax_b = pzmin_b = pzmax_b = gammamin_b = gammamax_b = false;
	xmin_b = xmax_b = pxmin_b = pxmax_b = Emin_b = Emax_b = thetamin_b = thetamax_b = true;
	nbin_y_b = nbin_z_b = nbin_py_b = nbin_pz_b = false;
	nbin_x_b = nbin_px_b = nbin_E_b = nbin_theta_b = true;
	fai_plot_xpx = fai_plot_Espec = fai_plot_Etheta = false;
	for (int i = 0; i < NPARAMETRI; i++) 
	{
		p_b[i] = true;
		p[i] = -1;
	}
	old_fortran_bin = false;
	endian_file = 0;
	endian_machine = is_big_endian();
	file_particelle_P = file_particelle_E = file_particelle_HI = file_particelle_LI = false;
	file_campi_Ex = file_campi_Ey = file_campi_Ez = file_campi_Bx = file_campi_By = file_campi_Bz = false;
	file_densita_elettroni = file_densita_protoni = file_densita_LI = file_densita_HI = false;
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


void parametri :: leggi_endian_ndv(std::ifstream& file_dat)
{
	std::string riga_persa;
	int trascura, ndv, i_end;
	std::getline(file_dat,riga_persa);
	file_dat >> trascura; // 1° parametro
	file_dat >> trascura; // 2° parametro
	file_dat >> trascura; // 3° parametro
	file_dat >> trascura; // 4° parametro
	file_dat >> trascura; // 5° parametro
	file_dat >> trascura; // 6° parametro
	file_dat >> trascura; // 7° parametro
	file_dat >> trascura; // 8° parametro
	file_dat >> trascura; // 9° parametro
	file_dat >> trascura; // 10° parametro
	file_dat >> trascura; // 11° parametro
	file_dat >> trascura; // 12° parametro
	file_dat >> trascura; // 13° parametro
	file_dat >> trascura; // 14° parametro
	file_dat >> trascura; // 15° parametro
	file_dat >> trascura; // 16° parametro
	file_dat >> trascura; // 17° parametro
	file_dat >> ndv;	  // 18° parametro
	file_dat >> trascura; // 19° parametro
	file_dat >> i_end;	  // 20° parametro

	p[NCOLUMNS] = ndv;
	p_b[NCOLUMNS] = false;

	if (ndv == 4 || ndv == 6)
	{
		p[WEIGHT] = 0;
		p_b[WEIGHT] = false;
	}
	else if (ndv == 5 || ndv == 7)
	{
		p[WEIGHT] = 1;
		p_b[WEIGHT] = false;
	}

	endian_file = (i_end-1);

}

void parametri :: chiedi_numero_colonne()
{
	int ncolonne;
	std::cout << "Il file contiene 6 o 7 colonne? ";
	std::cin >> ncolonne;
	if (ncolonne==6) p[WEIGHT] = 0;
	else if (ncolonne==7) p[WEIGHT] = 1;
	else exit(-5);
	p_b[WEIGHT] = false;
}


void parametri :: chiedi_endian_file()
{
	std::cout << "Il file e' little [x86] (0) o big [ppc] (1) endian? ";
	std::cin >> endian_file;
	if (endian_file != 1 && endian_file != 2) exit(-4);
}


void parametri :: check_filename(const char *nomefile)
{
	if (nomefile[0] == 'P')
	{
		if (nomefile[1] == 'r')
		{
			if (nomefile[2] == 'p')
			{
				massa_particella_MeV = (float) MP_MEV;
				file_particelle_P = true;
			}
		}
		else if (nomefile[1] == 'd')
		{
			file_densita_protoni = true;
			sprintf (support_label,"pden");
		}
	}
	else if (nomefile[0] == 'H')
	{
		if (nomefile[1] == 'i')
		{
			if (nomefile[2] == 'd')
			{
				file_densita_HI = true;
				sprintf (support_label,"hiden");
			}
		}
	}
	else if (nomefile[0] == 'L')
	{
		if (nomefile[1] == 'i')
		{
			if (nomefile[2] == 'd')
			{
				file_densita_LI = true;
				sprintf (support_label,"liden");
			}
		}
	}
	else if (nomefile[0] == 'E')
	{
		if (nomefile[1] == 'l')
		{
			if (nomefile[2] == 'p')
			{
				massa_particella_MeV = (float) ME_MEV;
				file_particelle_E = true;
			}
		}
		else if (nomefile[1] == 'x')
		{
			if (nomefile[2] == 'f')
			{
				file_campi_Ex = true;
				sprintf (support_label,"Ex");
			}
		}
		else if (nomefile[1] == 'y')
		{
			if (nomefile[2] == 'f')
			{
				file_campi_Ey = true;
				sprintf (support_label,"Ey");
			}
		}
		else if (nomefile[1] == 'z')
		{
			if (nomefile[2] == 'f')
			{
				file_campi_Ez = true;
				sprintf (support_label,"Ez");
			}
		}
		else if (nomefile[1] == 'd')
		{
			file_campi_Ez = true;
			sprintf (support_label,"eden");
		}
	}
	else if (nomefile[0] == 'B')
	{
		if (nomefile[1] == 'x')
		{
			if (nomefile[2] == 'f')
			{
				file_campi_Bx = true;
				sprintf (support_label,"Bx");
			}
		}
		else if (nomefile[1] == 'y')
		{
			if (nomefile[2] == 'f')
			{
				file_campi_By = true;
				sprintf (support_label,"By");
			}
		}
		else if (nomefile[1] == 'z')
		{
			if (nomefile[2] == 'f')
			{
				file_campi_Bz = true;
				sprintf (support_label,"Bz");
			}
		}
	}
	else
	{
		std::cout << "File non riconosciuto" << std::endl;
		exit(-15);
	}

}

void parametri :: leggi_da_file(const char *nomefile)
{
	bool failed_opening_file = true;
	std::ifstream fileParametri;
	fileParametri.open(nomefile);
	failed_opening_file = fileParametri.fail();
	if (failed_opening_file) 
	{
		std::cout << "Impossibile aprire il file dei parametri per il binnaggio" << std::endl;
		return;
	}

	std::string nomepar, evita, leggi;

	if (file_particelle_P || file_particelle_E || file_particelle_HI || file_particelle_LI)
	{
		if (p_b[FIND_MINMAX])
		{
			std::cout << "Vuoi cercare massimi e minimi? 0 per no, 1 per si': ";
			std::cin >> p[FIND_MINMAX];
			p_b[FIND_MINMAX] = false;
		}
		if (p_b[DO_BINNING])
		{
			std::cout << "Vuoi fare il binnaggio dei dati? 1 si', 0 no: ";
			std::cin >> p[DO_BINNING];
			p_b[DO_BINNING] = false;
		}
		if (p[DO_BINNING] == 1)
		{
			std::cout << "Quanti bin per asse vuoi usare? (consiglio: 100): ";
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
				/*
				else
				{
				std::cout << "Parametro " << nomepar << " non riconosciuto." << std::endl;
				}
				*/
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
		std::cout << "Vuoi l'output completo binario? 1 si', 0 no: ";
		std::cin >> p[OUT_BINARY];
		std::cout << "Vuoi l'output completo ascii? 1 si', 0 no: ";
		std::cin >> p[OUT_ASCII];
		std::cout << "Vuoi l'output dei parametri contenuti nel file? 1 si', 0 no: ";
		std::cin >> p[OUT_PARAMS];
	}
	/*
	else if(file_campi_Ex || file_campi_Ey || file_campi_Ez || file_campi_Bx || file_campi_By || file_campi_Bz)
	{
	std::cout << "Inserisci label per il file: (i.e. Ex) ";
	std::cin >> support_label;
	}
	*/
	fileParametri.close();
}


void parametri :: leggi_interattivo()
{
	int fare_xpx, fare_Espec, fare_Etheta;
	if (file_particelle_P || file_particelle_E || file_particelle_HI || file_particelle_LI)
	{
		std::cout << "Vuoi cercare massimi e minimi? 0 per no, 1 per si': ";
		std::cin >> p[FIND_MINMAX];
		std::cout << "Vuoi fare il binnaggio dei dati? 1 si', 0 no: ";
		std::cin >> p[DO_BINNING];
		if (p[DO_BINNING] == 1)
		{
			std::cout << "Vuoi fare il plot x-px? 1 si', 0 no: ";
			std::cin >> fare_xpx;
			std::cout << "Vuoi fare il plot E-theta? 1 si', 0 no: ";
			std::cin >> fare_Etheta;
			std::cout << "Vuoi fare lo spettro in energia? 1 si', 0 no: ";
			std::cin >> fare_Espec;
			fai_plot_xpx = fare_xpx;
			fai_plot_Espec = fare_Espec;
			fai_plot_Etheta = fare_Etheta;
			if (fare_xpx)
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
			}
			if (fare_Etheta)
			{
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
			}
			if (fare_Etheta || fare_Espec)
			{
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
		std::cout << "Hai scritto Emin = "<< Emin << "  Emax = "<< Emax<< std::endl;
		std::cout << "Vuoi l'output completo binario? 1 si', 0 no: ";
		std::cin >> p[OUT_BINARY];
		std::cout << "Vuoi l'output completo ascii? 1 si', 0 no: ";
		std::cin >> p[OUT_ASCII];
		std::cout << "Vuoi l'output dei parametri contenuti nel file? 1 si', 0 no: ";
		std::cin >> p[OUT_PARAMS];
	}
	/*
	else if(file_campi_Ex || file_campi_Ey || file_campi_Ez || file_campi_Bx || file_campi_By || file_campi_Bz)
	{
	std::cout << "Inserisci label per il file: (i.e. Ex) ";
	std::cin >> support_label;
	}
	*/
}

void parametri :: leggi_da_shell(int argc, const char *argv[])
{
	int ncolumns;
	for (int i = 2; i < argc; i++)	// * We will iterate over argv[] to get the parameters stored inside.
	{								// * Note that we're starting on 1 because we don't need to know the
		// * path of the program, which is stored in argv[0], and the input file,
		// * which is supposed to be given as the first argument and so is in argv[1]
		if (std::string(argv[i]) == "-swap")
		{
			p[SWAP] = 1;
		}
		else if (std::string(argv[i]) == "-ncol")
		{
			ncolumns = atoi(argv[i+1]);
			i++;
			if (ncolumns == 6) p[WEIGHT] = 0;
			else if (ncolumns == 7) p[WEIGHT] = 1;
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
		else if (std::string(argv[i]) == "-nbintheta")
		{
			nbin_gamma = atoi(argv[i+1]);
			nbin_gamma_b = false;
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

	std::cout << "---- organizza_minimi_massimi() -----Hai scritto Emin = "<< Emin << "  Emax = "<< Emax<< std::endl;
		
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
	std::cout << "---- organizza_minimi_massimi() -----Hai scritto Emin = "<< minimi[8] << "  Emax = "<< massimi[8] << std::endl;
	
}

bool parametri :: incompleto()
{
	bool test=true;
	//	if (!xmin_b && !xmax_b && !pxmin_b && !pxmax_b && !ymin_b && !ymax_b && !pymin_b && !pymax_b && !zmin_b && !zmax_b && !pzmin_b && !pzmax_b && 
	//		!Emin_b && !Emax_b && !nbin_x_b && !nbin_y_b && !nbin_z_b && !nbin_px_b && !nbin_py_b && !nbin_pz_b && !nbin_E_b)  test = false;
	//	if (!gammamin_b && !gammamax_b && !thetamin_b && !thetamax_b && !nbin_theta_b && !nbin_gamma_b && !test) test = false;

	test=false;
	return test;
}


#endif
