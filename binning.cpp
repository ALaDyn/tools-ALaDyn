#ifndef __BINNING_C
#define __BINNING_C

#include<cstdio>
#include<iostream>
#include<cstdlib>
#include<cmath>
#include<cstring>

class parametri_binnaggio
{
public:
	int nbin_x, nbin_y, nbin_z, nbin_px, nbin_py, nbin_pz, nbin_E, nbin_theta;
	float xmin, xmax, pxmin, pxmax, ymin, ymax, pymin, pymax, zmin, zmax, pzmin, pzmax, Emin, Emax, thetamin, thetamax;
//	float dim_x, dim_y, dim_z, dim_px, dim_py, dim_pz, dim_E, dim_theta;

	/* costruttore default - inizializza a zero per ora. Siccome tutto e' inizializzato a zero 
	non puo' nemmeno calcolare le le dimensioni dei bin!
	Verificare che non ci siano cose intelligenti da poter fare! */
	parametri_binnaggio()
	{
		nbin_x = nbin_px = nbin_y = nbin_z = nbin_py = nbin_pz = nbin_E = nbin_theta = 0;
		xmin = xmax = pxmin = pxmax = ymin = ymax = pymin = pymax = zmin = zmax = pzmin = pzmax = thetamin = thetamax = 0.0;
//		dim_x = dim_y = dim_z = dim_px = dim_py = dim_pz = dim_E = dim_theta = 0.0;
	}

	/* costruttore parametrico 1D */
	parametri_binnaggio(float xm, float xM, float pxm, float pxM, int nx, int npx)
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

};


#endif