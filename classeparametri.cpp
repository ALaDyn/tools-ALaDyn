#ifndef __BINNING_C
#define __BINNING_C

#include "leggi_binario_ALaDyn_fortran.h"


class parametri
{
public:
	int nbin_x, nbin_y, nbin_z, nbin_px, nbin_py, nbin_pz, nbin_E, nbin_theta;
	int p[NPARAMETRI];
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

	void leggi_da_file(char *nomefile)
	{
		std::cerr << "Non implementato" << std::endl;
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