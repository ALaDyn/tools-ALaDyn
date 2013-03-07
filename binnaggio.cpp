# ifndef __BINNING_CPP
# define __BINNING_CPP
# include "leggi_binario_ALaDyn_fortran.h"

_Binnaggio :: _Binnaggio(float * particelle, int npart, int ndv, Parametri * parametri, float ** data_binned, std::string binx, std::string biny)
{
	int binnare_su_x = 0, binnare_su_y = 0;
	int whichbin_x = 0, whichbin_y = 0;

	if (binx == "x") binnare_su_x = 0;
	else if (binx == "y") binnare_su_x = 1;
	else if (binx == "z") binnare_su_x = 2;
	else if (binx == "px") binnare_su_x = 3;
	else if (binx == "py") binnare_su_x = 4;
	else if (binx == "pz") binnare_su_x = 5;
	else if (binx == "gamma") binnare_su_x = 6;
	else if (binx == "theta") binnare_su_x = 7;
	else if (binx == "E") binnare_su_x = 8;
	else printf("variabile x non riconosciuta\n");

	if (biny == "x") binnare_su_y = 0;
	else if (biny == "y") binnare_su_y = 1;
	else if (biny == "z") binnare_su_y = 2;
	else if (biny == "px") binnare_su_y = 3;
	else if (biny == "py") binnare_su_y = 4;
	else if (biny == "pz") binnare_su_y = 5;
	else if (biny == "gamma") binnare_su_y = 6;
	else if (biny == "theta") binnare_su_y = 7;
	else if (biny == "E") binnare_su_y = 8;
	else printf("variabile y non riconosciuta\n");

	//	float z;
	float x, y, px, py, pz, gamma, theta, E;
	float dato_da_binnare_x = 0., dato_da_binnare_y = 0.;
	for (int i = 0; i < npart; i++)
	{
		if (ndv == 4 || ndv == 5)
		{
			x=*(particelle+i*ndv);
			y=*(particelle+i*ndv+1);
			px=*(particelle+i*ndv+2);
			py=*(particelle+i*ndv+3);
			gamma=(float)(sqrt(1.+px*px+py*py)-1.);				//gamma
			theta=(float)(atan2(py,px)*180./M_PI);				//theta
			E=(float)(gamma*parametri->massa_particella_MeV);	//energia
			if (binnare_su_x == 0) dato_da_binnare_x = x;
			else if (binnare_su_x == 1) dato_da_binnare_x = y;
			else if (binnare_su_x == 3) dato_da_binnare_x = px;
			else if (binnare_su_x == 4) dato_da_binnare_x = py;
			else if (binnare_su_x == 6) dato_da_binnare_x = gamma;
			else if (binnare_su_x == 7) dato_da_binnare_x = theta;
			else if (binnare_su_x == 8) dato_da_binnare_x = E;
			if (binnare_su_y == 0) dato_da_binnare_y = x;
			else if (binnare_su_y == 1) dato_da_binnare_y = y;
			else if (binnare_su_y == 3) dato_da_binnare_y = px;
			else if (binnare_su_y == 4) dato_da_binnare_y = py;
			else if (binnare_su_y == 6) dato_da_binnare_y = gamma;
			else if (binnare_su_y == 7) dato_da_binnare_y = theta;
			else if (binnare_su_y == 8) dato_da_binnare_y = E;
		}
		else if (ndv == 6 || ndv == 7)
		{
			px=*(particelle+i*ndv+3);
			py=*(particelle+i*ndv+4);
			pz=*(particelle+i*ndv+5);
			gamma=(float)(sqrt(1.+px*px+py*py+pz*pz)-1.);			//gamma
			theta=(float)(atan2(sqrt(py*py+pz*pz),px)*180./M_PI);	//theta nb: py e pz sono quelli trasversi in ALaDyn!
			E=(float)(gamma*parametri->massa_particella_MeV);		//energia
			if (binnare_su_x < 6) dato_da_binnare_x = *(particelle+i*ndv+binnare_su_x);
			else if (binnare_su_x == 6) dato_da_binnare_x = gamma;
			else if (binnare_su_x == 7) dato_da_binnare_x = theta;
			else if (binnare_su_x == 8) dato_da_binnare_x = E;
			if (binnare_su_y < 6) dato_da_binnare_y = *(particelle+i*ndv+binnare_su_y);
			else if (binnare_su_y == 6) dato_da_binnare_y = gamma;
			else if (binnare_su_y == 7) dato_da_binnare_y = theta;
			else if (binnare_su_y == 8) dato_da_binnare_y = E;
		}


		if (dato_da_binnare_x < parametri->minimi[binnare_su_x])
		{
			whichbin_x = 0;
		}
		else if (dato_da_binnare_x > parametri->massimi[binnare_su_x])
		{
			whichbin_x = parametri->nbin_x + 2;
		}
		else
		{
			whichbin_x = (int) (((dato_da_binnare_x - parametri->minimi[binnare_su_x]) / parametri->dimmi_dim(binnare_su_x)) +1.0);
		}
		if (dato_da_binnare_y < parametri->minimi[binnare_su_y])
		{
			whichbin_y = 0;
		}
		else if (dato_da_binnare_y > parametri->massimi[binnare_su_y])
		{
			whichbin_y = parametri->nbin_px + 2;
		}
		else
		{
			whichbin_y = (int) (((dato_da_binnare_y - parametri->minimi[binnare_su_y]) / parametri->dimmi_dim(binnare_su_y)) +1.0);
		}
		if (parametri->p[WEIGHT])	data_binned[whichbin_x][whichbin_y] += fabs(*(particelle+i*ndv+6));
		else						data_binned[whichbin_x][whichbin_y] += 1.0;
	}
}

_Binnaggio :: _Binnaggio(float * particelle, int npart, int ndv, Parametri * parametri, float * data_binned, std::string binx)
{
	int binnare_su_x = 0;
	int whichbin_x = 0;

	if (binx == "x") binnare_su_x = 0;
	else if (binx == "y") binnare_su_x = 1;
	else if (binx == "z") binnare_su_x = 2;
	else if (binx == "px") binnare_su_x = 3;
	else if (binx == "py") binnare_su_x = 4;
	else if (binx == "pz") binnare_su_x = 5;
	else if (binx == "gamma") binnare_su_x = 6;
	else if (binx == "theta") binnare_su_x = 7;
	else if (binx == "E") binnare_su_x = 8;
	else printf("variabile x non riconosciuta\n");

	//	float z;
	float x, y, px, py, pz, gamma, theta, E;
	float dato_da_binnare_x = 0.;
	//printf("AIUTO binnare_su_x=%i     min=%g    max=%g   dim=%g\n",binnare_su_x, parametri->minimi[binnare_su_x], parametri->massimi[binnare_su_x], parametri->dimmi_dim(binnare_su_x));
	fflush(stdout);
	for (int i = 0; i < npart; i++)
	{
		if (ndv == 4 || ndv == 5)
		{
			x=*(particelle+i*ndv);
			y=*(particelle+i*ndv+1);
			px=*(particelle+i*ndv+2);
			py=*(particelle+i*ndv+3);
			gamma=(float)(sqrt(1.+px*px+py*py)-1.);				//gamma
			theta=(float)(atan2(py,px)*180./M_PI);				//theta
			E=(float)(gamma*parametri->massa_particella_MeV);	//energia
			if (binnare_su_x == 0) dato_da_binnare_x = x;
			else if (binnare_su_x == 1) dato_da_binnare_x = y;
			else if (binnare_su_x == 3) dato_da_binnare_x = px;
			else if (binnare_su_x == 4) dato_da_binnare_x = py;
			else if (binnare_su_x == 6) dato_da_binnare_x = gamma;
			else if (binnare_su_x == 7) dato_da_binnare_x = theta;
			else if (binnare_su_x == 8) dato_da_binnare_x = E;
		}
		else if (ndv == 6 || ndv == 7)
		{
			px=*(particelle+i*ndv+3);
			py=*(particelle+i*ndv+4);
			pz=*(particelle+i*ndv+5);
			gamma=(float)(sqrt(1.+px*px+py*py+pz*pz)-1.);			//gamma
			theta=(float)(atan2(sqrt(py*py+pz*pz),px)*180./M_PI);	//theta nb: py e pz sono quelli trasversi in ALaDyn!
			E=(float)(gamma*parametri->massa_particella_MeV);		//energia
			if (binnare_su_x < 6) dato_da_binnare_x = *(particelle+i*ndv+binnare_su_x);
			else if (binnare_su_x == 6) dato_da_binnare_x = gamma;
			else if (binnare_su_x == 7) dato_da_binnare_x = theta;
			else if (binnare_su_x == 8) dato_da_binnare_x = E;
		}

		if (dato_da_binnare_x < parametri->minimi[binnare_su_x])
		{
			whichbin_x = 0;
		}
		else if (dato_da_binnare_x > parametri->massimi[binnare_su_x])
		{
			whichbin_x = parametri->nbin_x + 2;
		}
		else
		{
			whichbin_x = (int) (((dato_da_binnare_x - parametri->minimi[binnare_su_x]) / parametri->dimmi_dim(binnare_su_x)) +1.0);
		}
		if (parametri->p[WEIGHT])	data_binned[whichbin_x] += fabs(*(particelle+i*ndv+6));
		else						data_binned[whichbin_x] += 1.0;
	}
}

#endif

