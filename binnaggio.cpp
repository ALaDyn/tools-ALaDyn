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
	else if (binx == "thetaT") binnare_su_x = 9;
	else if (binx == "ty") binnare_su_x = 10;
	else if (binx == "tz") binnare_su_x = 11;
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
	else if (biny == "thetaT") binnare_su_y = 9;
	else if (biny == "ty") binnare_su_y = 10;
	else if (biny == "tz") binnare_su_y = 11;
	else printf("variabile y non riconosciuta\n");

	//	float z;
	float x, y, px, py, pz, gamma, theta, thetaT, E, ty, tz;
	float dato_da_binnare_x = 0., dato_da_binnare_y = 0.;
	for (int i = 0; i < npart; i++)
	{
		if (ndv == 4 || ndv == 5)
		{
			x=*(particelle+i*ndv);
			y=*(particelle+i*ndv+1);
			px=*(particelle+i*ndv+2);
			py=*(particelle+i*ndv+3);
			gamma=(float)(sqrt(1.+px*px+py*py)-1.);				//gamma-1
			theta=(float)(atan2(py,px)*180./M_PI);				//theta sgatto
			thetaT=(float) atan(sqrt((py*py/(px*px))));			//theta turch                                                                        
			//			thetaT=(float) atan(py/px);				//theta turch
			E=(float)(gamma*parametri->massa_particella_MeV);	//energia
			if (px > 0.0) ty = py/px;
			else ty = 1e15;
			tz = 0.0;
			if (binnare_su_x == 0) dato_da_binnare_x = x;
			else if (binnare_su_x == 1) dato_da_binnare_x = y;
			else if (binnare_su_x == 2) std::cout << "Unable to bin on z in 2D" << std::endl;
			else if (binnare_su_x == 3) dato_da_binnare_x = px;
			else if (binnare_su_x == 4) dato_da_binnare_x = py;
			else if (binnare_su_x == 5) std::cout << "Unable to bin on pz in 2D" << std::endl;
			else if (binnare_su_x == 6) dato_da_binnare_x = gamma;
			else if (binnare_su_x == 7) dato_da_binnare_x = theta;
			else if (binnare_su_x == 8) dato_da_binnare_x = E;
			else if (binnare_su_x == 9) dato_da_binnare_x = thetaT;
			else if (binnare_su_x == 10) dato_da_binnare_x = ty;
			else if (binnare_su_x == 11) std::cout << "Unable to bin on tz in 2D" << std::endl;
			if (binnare_su_y == 0) dato_da_binnare_y = x;
			else if (binnare_su_y == 1) dato_da_binnare_y = y;
			else if (binnare_su_y == 2) std::cout << "Unable to bin on z in 2D" << std::endl;
			else if (binnare_su_y == 3) dato_da_binnare_y = px;
			else if (binnare_su_y == 4) dato_da_binnare_y = py;
			else if (binnare_su_y == 5) std::cout << "Unable to bin on pz in 2D" << std::endl;
			else if (binnare_su_y == 6) dato_da_binnare_y = gamma;
			else if (binnare_su_y == 7) dato_da_binnare_y = theta;
			else if (binnare_su_y == 8) dato_da_binnare_y = E;
			else if (binnare_su_y == 9) dato_da_binnare_y = thetaT;
			else if (binnare_su_y == 10) dato_da_binnare_y = ty;
			else if (binnare_su_y == 11) std::cout << "Unable to bin on tz in 2D" << std::endl;
		}
		else if (ndv == 6 || ndv == 7)
		{
			px=*(particelle+i*ndv+3);
			py=*(particelle+i*ndv+4);
			pz=*(particelle+i*ndv+5);
			gamma=(float)(sqrt(1.+px*px+py*py+pz*pz)-1.);					//gamma-1
			theta=(float)(atan2(sqrt(py*py+pz*pz),px)*180./M_PI);			//theta sgatto
			thetaT=(float) atan(sqrt((py*py/(px*px)) + (pz*pz/(px*px))));	//theta turch
			E=(float)(gamma*parametri->massa_particella_MeV);				//energia
			if (px > 0.0)
			{
				ty = py/px;
				tz = pz/px;
			}
			else
			{
				ty = 1e15;
				tz = 1e15;
			}
			if (binnare_su_x < 6) dato_da_binnare_x = *(particelle+i*ndv+binnare_su_x);
			else if (binnare_su_x == 6) dato_da_binnare_x = gamma;
			else if (binnare_su_x == 7) dato_da_binnare_x = theta;
			else if (binnare_su_x == 8) dato_da_binnare_x = E;
			else if (binnare_su_x == 9) dato_da_binnare_x = thetaT;
			else if (binnare_su_x == 10) dato_da_binnare_x = ty;
			else if (binnare_su_x == 11) dato_da_binnare_x = tz;
			if (binnare_su_y < 6) dato_da_binnare_y = *(particelle+i*ndv+binnare_su_y);
			else if (binnare_su_y == 6) dato_da_binnare_y = gamma;
			else if (binnare_su_y == 7) dato_da_binnare_y = theta;
			else if (binnare_su_y == 8) dato_da_binnare_y = E;
			else if (binnare_su_y == 9) dato_da_binnare_y = thetaT;
			else if (binnare_su_y == 10) dato_da_binnare_y = ty;
			else if (binnare_su_y == 11) dato_da_binnare_y = tz;
		}


		if (dato_da_binnare_x < parametri->minimi[binnare_su_x])
		{
			whichbin_x = 0;
		}
		else if (dato_da_binnare_x > parametri->massimi[binnare_su_x])
		{
			whichbin_x = parametri->dimmi_nbin(binnare_su_x) + 2;
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
			whichbin_y = parametri->dimmi_nbin(binnare_su_y) + 2;
		}
		else
		{
			whichbin_y = (int) (((dato_da_binnare_y - parametri->minimi[binnare_su_y]) / parametri->dimmi_dim(binnare_su_y)) +1.0);
		}
		if (parametri->p[WEIGHT] && !parametri->overwrite_weight)	data_binned[whichbin_x][whichbin_y] += fabs(*(particelle+i*ndv+ndv-1));
		else														data_binned[whichbin_x][whichbin_y] += parametri->overwrite_weight_value;
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
	else if (binx == "thetaT") binnare_su_x = 9;
	else if (binx == "ty") binnare_su_x = 10;
	else if (binx == "tz") binnare_su_x = 11;
	else printf("variabile x non riconosciuta\n");

	//	float z;
	float x, y, px, py, pz, gamma, theta, thetaT, E, ty, tz;
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
			thetaT=(float) atan(py/px);							//theta turch
			E=(float)(gamma*parametri->massa_particella_MeV);	//energia
			if (px > 0.0) ty = py/px;
			else ty = 0.0;
			tz = 0.0;
			if (binnare_su_x == 0) dato_da_binnare_x = x;
			else if (binnare_su_x == 1) dato_da_binnare_x = y;
			else if (binnare_su_x == 2) std::cout << "Unable to bin on z in 2D" << std::endl;
			else if (binnare_su_x == 3) dato_da_binnare_x = px;
			else if (binnare_su_x == 4) dato_da_binnare_x = py;
			else if (binnare_su_x == 5) std::cout << "Unable to bin on pz in 2D" << std::endl;
			else if (binnare_su_x == 6) dato_da_binnare_x = gamma;
			else if (binnare_su_x == 7) dato_da_binnare_x = theta;
			else if (binnare_su_x == 8) dato_da_binnare_x = E;
			else if (binnare_su_x == 9) dato_da_binnare_x = thetaT;
			else if (binnare_su_x == 10) dato_da_binnare_x = ty;
			else if (binnare_su_x == 11) std::cout << "Unable to bin on tz in 2D" << std::endl;
		}
		else if (ndv == 6 || ndv == 7)
		{
			px=*(particelle+i*ndv+3);
			py=*(particelle+i*ndv+4);
			pz=*(particelle+i*ndv+5);
			gamma=(float)(sqrt(1.+px*px+py*py+pz*pz)-1.);			//gamma
			theta=(float)(atan2(sqrt(py*py+pz*pz),px)*180./M_PI);	//theta nb: py e pz sono quelli trasversi in ALaDyn!
			thetaT=(float) atan(py/px);							//theta turch
			E=(float)(gamma*parametri->massa_particella_MeV);		//energia
			if (px > 0.0)
			{
				ty = py/px;
				tz = pz/px;
			}
			else
			{
				ty = 0.0;
				tz = 0.0;
			}
			if (binnare_su_x < 6) dato_da_binnare_x = *(particelle+i*ndv+binnare_su_x);
			else if (binnare_su_x == 6) dato_da_binnare_x = gamma;
			else if (binnare_su_x == 7) dato_da_binnare_x = theta;
			else if (binnare_su_x == 8) dato_da_binnare_x = E;
			else if (binnare_su_x == 9) dato_da_binnare_x = thetaT;
			else if (binnare_su_x == 10) dato_da_binnare_x = ty;
			else if (binnare_su_x == 11) dato_da_binnare_x = tz;
		}

		if (dato_da_binnare_x < parametri->minimi[binnare_su_x])
		{
			whichbin_x = 0;
		}
		else if (dato_da_binnare_x > parametri->massimi[binnare_su_x])
		{
			whichbin_x = parametri->dimmi_nbin(binnare_su_x) + 2;
		}
		else
		{
			whichbin_x = (int) (((dato_da_binnare_x - parametri->minimi[binnare_su_x]) / parametri->dimmi_dim(binnare_su_x)) +1.0);
		}
		if (parametri->p[WEIGHT] && !parametri->overwrite_weight)	data_binned[whichbin_x] += fabs(*(particelle+i*ndv+ndv-1));
		else														data_binned[whichbin_x] += parametri->overwrite_weight_value;
	}
}


_Scrittura :: _Scrittura(Parametri * parametri, float ** data_binned, std::string x, std::string y, std::string nomefile_out)
{
	float xmin, xmax, dimx, ymin, ymax, dimy;
	int nbinx, nbiny;
	std::ofstream file_out;
	file_out.open(nomefile_out.c_str(), std::ofstream::out);

	if (x == "x")
	{
		dimx = (parametri->dimmi_dimx());
		xmin = (parametri->xmin)-dimx;
		xmax = (parametri->xmin);
		nbinx = (parametri->nbin_x+3);
	}
	else if (x == "y") 
	{
		dimx = (parametri->dimmi_dimy());
		xmin = (parametri->ymin)-dimx;
		xmax = (parametri->ymin);
		nbinx = (parametri->nbin_y+3);
	}
	else if (x == "z") 
	{
		dimx = (parametri->dimmi_dimz());
		xmin = (parametri->zmin)-dimx;
		xmax = (parametri->zmin);
		nbinx = (parametri->nbin_z+3);
	}
	else if (x == "ty") 
	{
		dimx = (parametri->dimmi_dimty());
		xmin = (parametri->tymin)-dimx;
		xmax = (parametri->tymin);
		nbinx = (parametri->nbin_ty+3);
	}
	else if (x == "tz") 
	{
		dimx = (parametri->dimmi_dimtz());
		xmin = (parametri->tzmin)-dimx;
		xmax = (parametri->tzmin);
		nbinx = (parametri->nbin_tz+3);
	}
	else if (x == "px") 
	{
		dimx = (parametri->dimmi_dimpx());
		xmin = (parametri->pxmin)-dimx;
		xmax = (parametri->pxmin);
		nbinx = (parametri->nbin_px+3);
	}
	else if (x == "py") 
	{
		dimx = (parametri->dimmi_dimpy());
		xmin = (parametri->pymin)-dimx;
		xmax = (parametri->pymin);
		nbinx = (parametri->nbin_py+3);
	}
	else if (x == "pz") 
	{
		dimx = (parametri->dimmi_dimpz());
		xmin = (parametri->pzmin)-dimx;
		xmax = (parametri->pzmin);
		nbinx = (parametri->nbin_pz+3);
	}
	else if (x == "gamma") 
	{
		dimx = (parametri->dimmi_dimgamma());
		xmin = (parametri->gammamin)-dimx;
		xmax = (parametri->gammamin);
		nbinx = (parametri->nbin_gamma+3);
	}
	else if (x == "theta") 
	{
		dimx = (parametri->dimmi_dimtheta());
		xmin = (parametri->thetamin)-dimx;
		xmax = (parametri->thetamin);
		nbinx = (parametri->nbin_theta+3);
	}
	else if (x == "E") 
	{
		dimx = (parametri->dimmi_dimE());
		xmin = (parametri->Emin)-dimx;
		xmax = (parametri->Emin);
		nbinx = (parametri->nbin_E+3);
	}
	else if (x == "thetaT") 
	{
		dimx = (parametri->dimmi_dimthetaT());
		xmin = (parametri->thetaTmin)-dimx;
		xmax = (parametri->thetaTmin);
		nbinx = (parametri->nbin_thetaT+3);
	}
	else printf("variabile x non riconosciuta\n");


	if (y == "x")
	{
		dimy = (parametri->dimmi_dimx());
		ymin = (parametri->xmin)-dimy;
		ymax = (parametri->xmin);
		nbiny = (parametri->nbin_x+3);
	}
	else if (y == "y") 
	{
		dimy = (parametri->dimmi_dimy());
		ymin = (parametri->ymin)-dimy;
		ymax = (parametri->ymin);
		nbiny = (parametri->nbin_y+3);
	}
	else if (y == "z") 
	{
		dimy = (parametri->dimmi_dimz());
		ymin = (parametri->zmin)-dimy;
		ymax = (parametri->zmin);
		nbiny = (parametri->nbin_z+3);
	}
	else if (y == "ty") 
	{
		dimy = (parametri->dimmi_dimty());
		ymin = (parametri->tymin)-dimy;
		ymax = (parametri->tymin);
		nbiny = (parametri->nbin_ty+3);
	}
	else if (y == "tz") 
	{
		dimy = (parametri->dimmi_dimtz());
		ymin = (parametri->tzmin)-dimy;
		ymax = (parametri->tzmin);
		nbiny = (parametri->nbin_tz+3);
	}
	else if (y == "px") 
	{
		dimy = (parametri->dimmi_dimpx());
		ymin = (parametri->pxmin)-dimy;
		ymax = (parametri->pxmin);
		nbiny = (parametri->nbin_px+3);
	}
	else if (y == "py") 
	{
		dimy = (parametri->dimmi_dimpy());
		ymin = (parametri->pymin)-dimy;
		ymax = (parametri->pymin);
		nbiny = (parametri->nbin_py+3);
	}
	else if (y == "pz") 
	{
		dimy = (parametri->dimmi_dimpz());
		ymin = (parametri->pzmin)-dimy;
		ymax = (parametri->pzmin);
		nbiny = (parametri->nbin_pz+3);
	}
	else if (y == "gamma") 
	{
		dimy = (parametri->dimmi_dimgamma());
		ymin = (parametri->gammamin)-dimy;
		ymax = (parametri->gammamin);
		nbiny = (parametri->nbin_gamma+3);
	}
	else if (y == "theta") 
	{
		dimy = (parametri->dimmi_dimtheta());
		ymin = (parametri->thetamin)-dimy;
		ymax = (parametri->thetamin);
		nbiny = (parametri->nbin_theta+3);
	}
	else if (y == "E") 
	{
		dimy = (parametri->dimmi_dimE());
		ymin = (parametri->Emin)-dimy;
		ymax = (parametri->Emin);
		nbiny = (parametri->nbin_E+3);
	}
	else if (y == "thetaT") 
	{
		dimy = (parametri->dimmi_dimthetaT());
		ymin = (parametri->thetaTmin)-dimy;
		ymax = (parametri->thetaTmin);
		nbiny = (parametri->nbin_thetaT+3);
	}
	else printf("variabile y non riconosciuta\n");


	float yminT = ymin;
	float ymaxT = ymax;

	for (int i = 0; i < nbinx; i++)
	{
		for (int j = 0; j < nbiny; j++)
		{
			file_out << std::setprecision(6) << xmin << "\t" << xmax << "\t" << yminT << "\t" << ymaxT << "\t" << data_binned[i][j] << std::endl;
			yminT += dimy;
			ymaxT += dimy;
		}
		xmin += dimx;
		xmax += dimx;
		yminT = ymin;
		ymaxT = ymax;
	}
	file_out.close();
}



_Scrittura :: _Scrittura(Parametri * parametri, float * data_binned, std::string x, std::string nomefile_out)
{
	float xmin, xmax, dimx;
	int nbinx;
	std::ofstream file_out;
	file_out.open(nomefile_out.c_str(), std::ofstream::out);

	if (x == "x")
	{
		dimx = (parametri->dimmi_dimx());
		xmin = (parametri->xmin)-dimx;
		xmax = (parametri->xmin);
		nbinx = (parametri->nbin_x+3);
	}
	else if (x == "y") 
	{
		dimx = (parametri->dimmi_dimy());
		xmin = (parametri->ymin)-dimx;
		xmax = (parametri->ymin);
		nbinx = (parametri->nbin_y+3);
	}
	else if (x == "z") 
	{
		dimx = (parametri->dimmi_dimz());
		xmin = (parametri->zmin)-dimx;
		xmax = (parametri->zmin);
		nbinx = (parametri->nbin_z+3);
	}
	else if (x == "px") 
	{
		dimx = (parametri->dimmi_dimpx());
		xmin = (parametri->pxmin)-dimx;
		xmax = (parametri->pxmin);
		nbinx = (parametri->nbin_px+3);
	}
	else if (x == "py") 
	{
		dimx = (parametri->dimmi_dimpy());
		xmin = (parametri->pymin)-dimx;
		xmax = (parametri->pymin);
		nbinx = (parametri->nbin_py+3);
	}
	else if (x == "pz") 
	{
		dimx = (parametri->dimmi_dimpz());
		xmin = (parametri->pzmin)-dimx;
		xmax = (parametri->pzmin);
		nbinx = (parametri->nbin_pz+3);
	}
	else if (x == "gamma") 
	{
		dimx = (parametri->dimmi_dimgamma());
		xmin = (parametri->gammamin)-dimx;
		xmax = (parametri->gammamin);
		nbinx = (parametri->nbin_gamma+3);
	}
	else if (x == "theta") 
	{
		dimx = (parametri->dimmi_dimtheta());
		xmin = (parametri->thetamin)-dimx;
		xmax = (parametri->thetamin);
		nbinx = (parametri->nbin_theta+3);
	}
	else if (x == "E") 
	{
		dimx = (parametri->dimmi_dimE());
		xmin = (parametri->Emin)-dimx;
		xmax = (parametri->Emin);
		nbinx = (parametri->nbin_E+3);
	}
	else if (x == "thetaT") 
	{
		dimx = (parametri->dimmi_dimthetaT());
		xmin = (parametri->thetaTmin)-dimx;
		xmax = (parametri->thetaTmin);
		nbinx = (parametri->nbin_thetaT+3);
	}
	else if (x == "ty") 
	{
		dimx = (parametri->dimmi_dimty());
		xmin = (parametri->tymin)-dimx;
		xmax = (parametri->tymin);
		nbinx = (parametri->nbin_ty+3);
	}
	else if (x == "tz") 
	{
		dimx = (parametri->dimmi_dimtz());
		xmin = (parametri->tzmin)-dimx;
		xmax = (parametri->tzmin);
		nbinx = (parametri->nbin_tz+3);
	}
	else printf("variabile x non riconosciuta\n");


	for (int i = 0; i < nbinx; i++)
	{
		file_out << std::setprecision(6) << xmin << "\t" << xmax << "\t" << data_binned[i] << std::endl;
		xmin += dimx;
		xmax += dimx;
	}
	file_out.close();
}

#endif

