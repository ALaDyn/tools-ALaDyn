#ifndef __LEGGI_CAMPI_C
#define __LEGGI_CAMPI_C
#include "leggi_binario_ALaDyn_fortran.h"


int leggi_campi(int argc, const char** argv, Parametri * parametri)
{
	std::string basefilename = std::string(argv[1]);
	std::ostringstream nomefile_bin;

	const int out_swap = parametri->p[SWAP];
	const int out_parameters = parametri->p[OUT_PARAMS];
	const int out_vtk = parametri->p[OUT_VTK];
	const int out_vtk_nostretch = parametri->p[OUT_VTK_NOSTRETCH];
	const int out_cutx = parametri->p[OUT_CUTX];
	const int out_cuty = parametri->p[OUT_CUTY];
	const int out_cutz = parametri->p[OUT_CUTZ];
	const int out_lineoutx = parametri->p[OUT_LINEOUT_X];
	const int out_2d = parametri->p[OUT_GRID2D];

	size_t npunti_x = parametri->npunti_x_ricampionati;
	size_t npunti_y = parametri->npunti_y_ricampionati;
	size_t npunti_z = parametri->npunti_z_ricampionati;
	size_t npunti_x_per_cpu = parametri->npx_per_cpu;
	size_t npunti_y_per_cpu = parametri->npy_per_cpu;
	size_t npunti_z_per_cpu = parametri->npz_per_cpu;
	int ncpu_x = parametri->ncpu_x;
	int ncpu_y = parametri->ncpu_y;
	int ncpu_z = parametri->ncpu_z;
	int span = parametri->span;
	float xminimo = parametri->xmin;
	float xmaximo = parametri->xmax;
	float yminimo = parametri->ymin;
	float ymaximo = parametri->ymax;
	float zminimo = parametri->zmin;
	float zmaximo = parametri->zmax;
	float taglio;
	std::vector<float> cutx, cuty, cutz;
	std::vector<size_t> gridIndex_cutx, gridIndex_cuty, gridIndex_cutz;

	int indice_multifile=0;

	int N_param = 0, np_loc, npoint_loc[3] = {0,0,0}, loc_size = 0;
	int segnoy = 0, segnoz = 0, buff = 0;
	int nx = 0, ibx = 0, iby = 0, ibz = 0, model = 0, dmodel = 0, nsp = 0, lpord = 0, deord = 0, npe = 0, fvar = 0;
	int npe_y = 0, npe_z = 0, nxloc = 0, nx1 = 0, ny1 = 0, nyloc = 0, nz1 = 0, nzloc = 0;

	float tnow = 0., w0x = 0., w0y = 0., nrat = 0., a0 = 0., lam0 = 0., B0 = 0., ompe = 0.;
	float xt_in = 0., xt_end = 0., charge = 0., mass = 0., xmin = 0., xmax = 0., ymin = 0., ymax = 0., zmin = 0., zmax = 0.;
	float E0 = 0., dx = 0., dy = 0., dz = 0., xx = 0., yy = 0.;


	int *int_param;
	float *real_param;
	float *field,*buffer;
	double *x_lineout;
	char nomefile_parametri[MAX_LENGTH_FILENAME];
	char nomefile_campi[MAX_LENGTH_FILENAME];

	std::FILE *file_in;
	std::FILE *parameters;
	std::FILE *clean_fields;

	size_t fread_size = 0;


	if (parametri->old_fortran_bin)
	{
		nomefile_bin << basefilename << ".bin";
		file_in=fopen(nomefile_bin.str().c_str(), "rb");
		if(file_in==NULL) std::cout << "Unable to open file!" << std::endl;
		else std::cout << "File opened to read parameters!" << std::endl;

		fread_size = std::fread(&buff,sizeof(int),1,file_in);
		fread_size = std::fread(&N_param,sizeof(int),1,file_in);
		fread_size = std::fread(&buff,sizeof(int),1,file_in);

		if(out_swap) swap_endian_i(&N_param,1);

		int_param = new int[N_param];
		real_param = new float[N_param];
		//	int_param=(int*)malloc(N_param*sizeof(int));
		//	real_param=(float*)malloc(N_param*sizeof(float));
		fread_size = std::fread(&buff,sizeof(int),1,file_in);
		fread_size = std::fread(int_param,sizeof(int),N_param,file_in);
		fread_size = std::fread(&buff,sizeof(int),1,file_in);
		fread_size = std::fread(&buff,sizeof(int),1,file_in);
		fread_size = std::fread(real_param,sizeof(float),N_param,file_in);
		fread_size = std::fread(&buff,sizeof(int),1,file_in);

		if (out_swap) swap_endian_i(int_param,N_param);
		if (out_swap) swap_endian_f(real_param,N_param);

		fclose(file_in);

		npe_y=int_param[0];
		npe_z=int_param[1];
		npe=npe_y*npe_z;
		nx=int_param[2];
		nx1=int_param[3];
		ny1=int_param[4];
		nyloc=int_param[5];
		nz1=int_param[6];
		nzloc=int_param[7];
		ibx=int_param[8];
		iby=int_param[9];
		ibz=int_param[10];
		model=int_param[11];  //modello di laser utilizzato
		dmodel=int_param[12]; //modello di condizioni iniziali
		nsp=int_param[13];    //numero di speci
		np_loc=int_param[14];  //numero di componenti dello spazio dei momenti
		lpord=int_param[15]; //ordine dello schema leapfrog
		deord=int_param[16]; //ordine derivate
		fvar=int_param[17]; 
		tnow=real_param[0];  //tempo dell'output
		xmin=real_param[1];  //estremi della griglia
		xmax=real_param[2];  //estremi della griglia
		ymin=real_param[3];  //estremi della griglia
		ymax=real_param[4];  //estremi della griglia
		zmin=real_param[5];  //estremi della griglia
		zmax=real_param[6];  //estremi della griglia
		w0x=real_param[7];      //waist del laser in x
		w0y=real_param[8];      //waist del laser in y
		nrat=real_param[9];     //n orver n critical
		a0=real_param[10];      // a0 laser
		lam0=real_param[11];    // lambda
		E0=real_param[12];      //conversione da campi numerici a TV/m
		B0=E0+(float)(33.3);
		ompe=real_param[13];    //costante accoppiamento correnti campi
		xt_in=real_param[14];   //inizio plasma
		xt_end=real_param[15];
		charge=real_param[16];  //carica particella su carica elettrone
		mass=real_param[17];    //massa particelle su massa elettrone
	}
	else
	{
		std::cout << "Reading parameters from dat files!" << std::endl;
		npe_y=parametri->intpar[0];
		npe_z=parametri->intpar[1];
		npe=npe_y*npe_z;
		nx=parametri->intpar[2];
		nx1=parametri->intpar[3];
		ny1=parametri->intpar[4];
		nyloc=parametri->intpar[5];
		nz1=parametri->intpar[6];
		nzloc=parametri->intpar[7];
		ibx=parametri->intpar[8];
		iby=parametri->intpar[9];
		ibz=parametri->intpar[10];
		model=parametri->intpar[11];
		dmodel=parametri->intpar[12];
		nsp=parametri->intpar[13];
		np_loc=parametri->intpar[14];
		lpord=parametri->intpar[15];
		deord=parametri->intpar[16];
		fvar=parametri->intpar[17];
		tnow=parametri->realpar[0];
		xmin=parametri->realpar[1];
		xmax=parametri->realpar[2];
		ymin=parametri->realpar[3];
		ymax=parametri->realpar[4];
		zmin=parametri->realpar[5];
		zmax=parametri->realpar[6];
		w0x=parametri->realpar[7];
		w0y=parametri->realpar[8];
		nrat=parametri->realpar[9];
		a0=parametri->realpar[10];
		lam0=parametri->realpar[11];
		E0=parametri->realpar[12];
		B0=E0+(float)(100.0/3.0);
		ompe=parametri->realpar[13];
		xt_in=parametri->realpar[14];
		xt_end=parametri->realpar[15];
		charge=parametri->realpar[16];
		mass=parametri->realpar[17];
	}


	printf("=========INIZIO LETTURE==========\n");
	size_t prodotto=nx1*ny1*nz1;
	printf("nx1*ny1*nz1: %i %i %i = %.10e\n",nx1,ny1,nz1,(double)prodotto);
	fflush(stdout);

	/* Quanto segue non sarebbe necessario se vecchie e nuove variabili fossero unificate*/
	npunti_x = nx1;
	npunti_y = ny1;
	npunti_z = nz1;
	ncpu_x = 1;
	ncpu_y = npe_y;
	ncpu_z = npe_z;
	npunti_x_per_cpu = nx1;
	npunti_y_per_cpu = nyloc;
	npunti_z_per_cpu = nzloc;
	xminimo = xmin;
	xmaximo = xmax;
	yminimo = ymin;
	ymaximo = ymax;
	zminimo = zmin;
	zmaximo = zmax;


	field = new float[npunti_x*npunti_y*npunti_z];
	// field=(float*)malloc(npunti_x*npunti_y*npunti_z*sizeof(float));
	x_lineout=new double[npunti_x];

	if (!parametri->nuovi_dati_su_griglia)
	{
		nomefile_bin.str("");
		nomefile_bin << basefilename << ".bin";
		file_in=fopen(nomefile_bin.str().c_str(), "rb");
		if(file_in==NULL) std::cout << "Unable to open file!" << std::endl;
		else std::cout << "File opened to read data!" << std::endl;

		fread_size = std::fread(&buff,sizeof(int),1,file_in);
		fread_size = std::fread(&N_param,sizeof(int),1,file_in);
		fread_size = std::fread(&buff,sizeof(int),1,file_in);

		if(out_swap) swap_endian_i(&N_param,1);
		int_param = new int[N_param];
		real_param = new float[N_param];

		fread_size = std::fread(&buff,sizeof(int),1,file_in);
		fread_size = std::fread(int_param,sizeof(int),N_param,file_in);
		fread_size = std::fread(&buff,sizeof(int),1,file_in);
		fread_size = std::fread(&buff,sizeof(int),1,file_in);
		fread_size = std::fread(real_param,sizeof(float),N_param,file_in);
		fread_size = std::fread(&buff,sizeof(int),1,file_in);

		for(int ipz = 0; ipz < ncpu_z; ipz++)
		{
			segnoy=0;
			for(int ipy=0; ipy < ncpu_y; ipy++)
			{
				fread_size = std::fread(&buff,sizeof(int),1,file_in);
				fread_size = std::fread(npoint_loc,sizeof(int),3,file_in);
				fread_size = std::fread(&buff,sizeof(int),1,file_in);

				if(out_swap) swap_endian_i(npoint_loc,3);

				loc_size=npoint_loc[0]*npoint_loc[1]*npoint_loc[2];

				nxloc=npoint_loc[0];	// ma non dovrebbero essere gia' noti? possono variare da cpu a cpu? non credo...
				nyloc=npoint_loc[1];	// comunque a questo punto lascio i "vecchi" nxloc, nyloc ed nzloc, in attesa di rimuoverli in futuro con le versioni aggiornate
				nzloc=npoint_loc[2];	// che leggono dal binario se necessario oppure dal dat allegato per i nuovi

#ifdef ENABLE_DEBUG
				printf("processore ipz=%i/%i  ipy=%i/%i     segnoz=%i     segnoy=%i\n",ipz,npe_z, ipy,npe_y,segnoz,segnoy );
#else
				printf("processore ipz=%i/%i  ipy=%i/%i     segnoz=%i     segnoy=%i\r",ipz,npe_z, ipy,npe_y,segnoz,segnoy );
#endif
				fflush(stdout);

				buffer = new float[loc_size];
				//	buffer=(float*)malloc(loc_size*sizeof(float));
				fread_size = std::fread(&buff,sizeof(int),1,file_in);
				fread_size = std::fread(buffer,sizeof(float),loc_size,file_in);
				fread_size = std::fread(&buff,sizeof(int),1,file_in);

				if(out_swap) swap_endian_f(buffer,loc_size);

				for(size_t k=0; k<nzloc; k++)
					for(size_t j=0; j<nyloc; j++)
						for(size_t i=0; i<nxloc; i++)
							field[i+(j+segnoy)*npunti_x+(k+segnoz)*npunti_x*npunti_y]=buffer[i+j*nxloc+k*nxloc*nyloc];
				segnoy += nyloc;
			} 
			segnoz += nzloc;
			delete[] buffer;
			buffer = NULL;
		}

		// leggiamo ora le coordinate dei punti di griglia, presenti solo nelle versioni che possono prevedere griglia stretchata e che ancora non la scrivevano nel .dat
		// se presenti, sovrascrivono quelle lette o precostruite (se non trovate nel file .dat) dalle routine dei parametri


		fread_size = std::fread(&buff,sizeof(int),1,file_in);	// facciamo il test sul buffer Fortran della prima coordinata;
		// se esiste, non e' necessario tornare indietro perche' il buffer fortran che precede i dati non e' di alcun interesse

		if(!std::feof(file_in))
		{
			float *x_coordinates, *y_coordinates, *z_coordinates;
			x_coordinates=new float[npunti_x];
			y_coordinates=new float[npunti_y];
			z_coordinates=new float[npunti_z];
			fread_size = std::fread(x_coordinates,sizeof(float),npunti_x,file_in);
			fread_size = std::fread(&buff,sizeof(int),1,file_in);
			fread_size = std::fread(&buff,sizeof(int),1,file_in);
			fread_size = std::fread(y_coordinates,sizeof(float),npunti_y,file_in);
			fread_size = std::fread(&buff,sizeof(int),1,file_in);
			fread_size = std::fread(&buff,sizeof(int),1,file_in);
			fread_size = std::fread(z_coordinates,sizeof(float),npunti_z,file_in);
			fread_size = std::fread(&buff,sizeof(int),1,file_in);

			if(out_swap)
			{
				swap_endian_f(x_coordinates,npunti_x);
				swap_endian_f(y_coordinates,npunti_y);
				swap_endian_f(z_coordinates,npunti_z);
			}

			parametri->xcoord.resize(npunti_x,0);
			parametri->ycoord.resize(npunti_y,0);
			parametri->zcoord.resize(npunti_z,0);

			for (int i = 0; i < npunti_x; i++)
				parametri->xcoord[i] = x_coordinates[i];
			for (int i = 0; i < npunti_y; i++)
				parametri->ycoord[i] = y_coordinates[i];
			for (int i = 0; i < npunti_z; i++)
				parametri->zcoord[i] = z_coordinates[i];
		}
		else parametri->stretched_grid = false;

		printf("=========FINE LETTURE==========\n");
		fflush(stdout);
	}
	else
	{
		if (!parametri->multifile)
		{
			nomefile_bin.str("");
			nomefile_bin << std::string(argv[1]) << ".bin";
			file_in=fopen(nomefile_bin.str().c_str(), "rb");
			if(file_in==NULL) std::cout << "Unable to open file!" << std::endl;
			else std::cout << "File opened to read data!" << std::endl;
			int header_size = 3;
			int * header = new int[header_size];
			for(int ipz = 0; ipz < ncpu_z; ipz++)
			{
				segnoy=0;
				for(int ipy=0; ipy < ncpu_y; ipy++)
				{
					fread_size = std::fread(header,sizeof(int),header_size,file_in);
					if(out_swap) swap_endian_i(header,header_size);

					loc_size=header[0]*header[1]*header[2];

					nxloc=header[0];
					nyloc=header[1];
					nzloc=header[2];

#ifdef ENABLE_DEBUG
					printf("file %i, processore ipy=%i/%i, reading %i elements\n", indice_multifile, ipy, npe_y, loc_size);
#else
					printf("file %i, processore ipy=%i/%i\r", indice_multifile, ipy, npe_y);
#endif
					fflush(stdout);

					buffer = new float[loc_size];
					fread_size = std::fread(buffer,sizeof(float),loc_size,file_in);

					if(out_swap) swap_endian_f(buffer,loc_size);

					for(size_t k=0; k<nzloc; k++)
						for(size_t j=0; j<nyloc; j++)
							for(size_t i=0; i<nxloc; i++)
								field[i+(j+segnoy)*npunti_x+(k+segnoz)*npunti_x*npunti_y]=buffer[i+j*nxloc+k*nxloc*nyloc];
					segnoy += nyloc;
				}
				segnoz += nzloc;
				delete[] buffer;
				buffer = NULL;
			}

			printf("=========FINE LETTURE==========\n");
			fflush(stdout);
		}
		else
		{
			int header_size = 3;
			int * header = new int[header_size];
			while(1)
			{
				nomefile_bin.str("");
				nomefile_bin << std::string(argv[1]) << "_" << std::setfill('0') << std::setw(3) << indice_multifile <<".bin";
				file_in=fopen(nomefile_bin.str().c_str(), "rb");
				if (file_in == NULL)
				{
					std::cout << "End of files!" << std::endl;
					break;
				}
				else std::cout << "Opened file #" << indice_multifile << " to read data!" << std::endl;

				segnoy=0;
				for(int ipy=0; ipy < ncpu_y; ipy++)
				{
					fread_size = std::fread(header,sizeof(int),header_size,file_in);
					if(out_swap) swap_endian_i(header,header_size);

					loc_size=header[0]*header[1]*header[2];

					nxloc=header[0];
					nyloc=header[1];
					nzloc=header[2];

#ifdef ENABLE_DEBUG
					printf("file %i, processore ipy=%i/%i, reading %i elements\n", indice_multifile, ipy, npe_y, loc_size);
#else
					printf("file %i, processore ipy=%i/%i\r", indice_multifile, ipy, npe_y);
#endif
					fflush(stdout);

					buffer = new float[loc_size];
					fread_size = std::fread(buffer,sizeof(float),loc_size,file_in);

					if(out_swap) swap_endian_f(buffer,loc_size);

					for(size_t k=0; k<nzloc; k++)
						for(size_t j=0; j<nyloc; j++)
							for(size_t i=0; i<nxloc; i++)
								field[i+(j+segnoy)*npunti_x+(k+segnoz)*npunti_x*npunti_y]=buffer[i+j*nxloc+k*nxloc*nyloc];
					segnoy += nyloc;
				}
				segnoz += nzloc;
				delete[] buffer;
				buffer = NULL;
				indice_multifile++;
				fclose(file_in);
			}

			printf("=========FINE LETTURE==========\n");
			fflush(stdout);
		}
	}


	if(npunti_z == 1 && out_2d)
	{
		printf("\nScrittura file gnuplot 2D\n");
		sprintf(nomefile_campi,"%s.txt",argv[1]);
		clean_fields=fopen(nomefile_campi, "w");
		printf("\nWriting the fields file 2D (not vtk)\n");

		//output per gnuplot (x:y:valore) compatibile con programmino passe_par_tout togliendo i #
		fprintf(clean_fields,"#%i\n#%i\n#%i\n",npunti_x, npunti_y, npunti_z); 
		fprintf(clean_fields,"#%f %f\n#%f %f\n",xmin, ymin, xmax, ymax);
		for(int j = 0; j < npunti_y; j++)
		{
			for(int i = 0; i < npunti_x; i++)
			{
				xx=parametri->xcoord[i];//xmin+dx*i;
				yy=parametri->ycoord[j];//ymin+dy*j;
				fprintf(clean_fields,"%.4g %.4g %.4g\n",xx, yy, field[i+j*npunti_x]);
			}
		}
		fclose(clean_fields);
	}


	if(npunti_z == 1 && out_lineoutx)
	{
		int myj;
		printf("\nScrittura lineout 1D\n");
		sprintf(nomefile_campi,"%s_lineout.txt",argv[1]);
		clean_fields=fopen(nomefile_campi, "w");
		printf("\nWriting the lineout file 1D (not vtk)\n");

		for(int i = 0; i < npunti_x; i++) x_lineout[i]=0;

		for(int j = 0; j < npunti_y; j++)
		{
			if(parametri->ycoord[j]>=0)
			{
				myj=j;
				break;
			}
		}

		fprintf(clean_fields,"#");
		for(int j = myj-span; j < (myj+span+1); j++)
		{
			fprintf(clean_fields,"%.4g\t",parametri->ycoord[j]);
			for(int i = 0; i < npunti_x; i++)
			{
				x_lineout[i]+=field[i+j*npunti_x]/(2*span+1.0);
			}
		}
		fprintf(clean_fields,"\n");
		for(int i = 0; i < npunti_x; i++)
		{
			xx=parametri->xcoord[i];//xmin+dx*i;
			fprintf(clean_fields,"%.4g %.4g\n",xx, x_lineout[i]);
		}
		fclose(clean_fields);
	}

	if(npunti_z > 1 && out_lineoutx)
	{
		int myj;
		printf("\nScrittura lineout 1D\n");
		sprintf(nomefile_campi,"%s_lineout.txt",argv[1]);
		clean_fields=fopen(nomefile_campi, "w");
		printf("\nWriting the lineout file 1D (not vtk)\n");

		for(int i = 0; i < npunti_x; i++)
			x_lineout[i]=0;
		for(int j = 0; j < npunti_y; j++)
		{
			if(parametri->ycoord[j]>=0)
			{
				myj=j;
				break;
			}
		}
		fprintf(clean_fields,"#");
		for(size_t k = myj-span; k < (myj+span+1); k++)
		{
			fprintf(clean_fields,"%.4g\t",parametri->zcoord[k]);

			for(size_t j = myj-span; j < (myj+span+1); j++)
			{
				for(size_t i = 0; i < npunti_x; i++)
				{
					x_lineout[i]+=field[i+j*npunti_x+k*npunti_x*npunti_y]/((2*span+1.0)*(2*span+1.0));
				}
			}
		}
		fprintf(clean_fields,"\n");
		for(int i = 0; i < npunti_x; i++)
		{
			xx=parametri->xcoord[i];//xmin+dx*i;
			fprintf(clean_fields,"%.4g %.4g\n",xx, x_lineout[i]);
		}
		fclose(clean_fields);
	}


	if (npunti_z > 1 && out_cutz)
	{
		if (parametri->posizioni_taglio_griglia_z.size() == 0)
		{
			taglio = parametri->zcoord[npunti_z/2];
			cutz.push_back(taglio);
			gridIndex_cutz.push_back(npunti_z/2);
		}
		else
		{
			taglio = zminimo;
			int i = 0;
			for (size_t j = 0; j < parametri->posizioni_taglio_griglia_z.size(); j++)
			{
				while (taglio < parametri->posizioni_taglio_griglia_z.at(j)) taglio = parametri->zcoord[i], i++;
				cutz.push_back(taglio);
				gridIndex_cutz.push_back(i-1);
				i = 0;
			}
		}
		for (size_t n = 0; n < cutz.size(); n++)
		{
			sprintf(nomefile_campi,"%s_cutz_%g.txt",argv[1], cutz[n]);
			printf("\nScrittura file gnuplot taglio z=%g\n",cutz[n]);
			clean_fields=fopen(nomefile_campi, "wb");
			printf("\nWriting the fields file 2D (not vtk)\n");
			//output per gnuplot (x:y:valore) compatibile con programmino passe_par_tout togliendo i #
			fprintf(clean_fields,"# 2D cut at z=%g\n", cutz[n]); 
			fprintf(clean_fields,"# %i\n#%i\n#%i\n",npunti_x, npunti_y, 1); 
			fprintf(clean_fields,"#%f %f\n#%f %f\n",xmin, ymin, xmax, ymax);
			size_t k = gridIndex_cutz[n];
			for(size_t j = 0; j < npunti_y; j++)
			{
				for(size_t i = 0; i < npunti_x; i++)
				{
					xx=parametri->xcoord[i];//xx=xmin+dx*i;
					yy=parametri->ycoord[j];//yy=ymin+dy*j;
					fprintf(clean_fields,"%.4g %.4g %.4g\n",xx, yy, field[i+j*npunti_x+k*npunti_x*npunti_y]);
				}
			}
			fclose(clean_fields);
		}
	}

	if (npunti_z > 1 && out_cuty)
	{
		if (parametri->posizioni_taglio_griglia_y.size() == 0)
		{
			taglio = parametri->ycoord[npunti_y/2];
			cuty.push_back(taglio);
			gridIndex_cuty.push_back(npunti_y/2);
		}
		else
		{
			taglio = yminimo;
			int i = 0;
			for (size_t j = 0; j < parametri->posizioni_taglio_griglia_y.size(); j++)
			{
				while (taglio < parametri->posizioni_taglio_griglia_y.at(j)) taglio = parametri->ycoord[i], i++;
				cuty.push_back(taglio);
				gridIndex_cuty.push_back(i-1);
				i = 0;
			}
		}
		for (size_t n = 0; n < cuty.size(); n++)
		{
			sprintf(nomefile_campi,"%s_cuty_%g.txt",argv[1], cuty[n]);
			printf("\nScrittura file gnuplot taglio y=%g\n",cuty[n]);
			clean_fields=fopen(nomefile_campi, "wb");
			printf("\nWriting the fields file 2D (not vtk)\n");
			//output per gnuplot (x:y:valore) compatibile con programmino passe_par_tout togliendo i #
			fprintf(clean_fields,"# 2D cut at y=%g\n", cuty[n]); 
			fprintf(clean_fields,"# %i\n#%i\n#%i\n",npunti_x, npunti_z, 1); 
			fprintf(clean_fields,"#%f %f\n#%f %f\n",xmin, zmin, xmax, zmax);
			size_t j = gridIndex_cuty[n];
			for(size_t k = 0; k < npunti_z; k++)
			{
				for(size_t i = 0; i < npunti_x; i++)
				{
					xx=parametri->xcoord[i];//xx=xmin+dx*i;
					yy=parametri->ycoord[k];//yy=ymin+dy*j;
					fprintf(clean_fields,"%.4g %.4g %.4g\n",xx, yy, field[i+j*npunti_x+k*npunti_x*npunti_y]);
				}
			}
			fclose(clean_fields);
		}
	}

	if (npunti_z > 1 && out_cutx)
	{
		if (parametri->posizioni_taglio_griglia_x.size() == 0)
		{
			taglio = parametri->xcoord[npunti_x/2];
			cutx.push_back(taglio);
			gridIndex_cutx.push_back(npunti_x/2);
		}
		else
		{
			taglio = xminimo;
			int i = 0;
			for (size_t j = 0; j < parametri->posizioni_taglio_griglia_x.size(); j++)
			{
				while (taglio < parametri->posizioni_taglio_griglia_x.at(j)) taglio = parametri->xcoord[i], i++;
				cutx.push_back(taglio);
				gridIndex_cutx.push_back(i-1);
				i = 0;
			}
		}
		for (size_t n = 0; n < cutx.size(); n++)
		{
			sprintf(nomefile_campi,"%s_cutx_%g.txt",argv[1], cutx[n]);
			printf("\nScrittura file gnuplot taglio x=%g\n",cutx[n]);
			clean_fields=fopen(nomefile_campi, "wb");
			printf("\nWriting the fields file 2D (not vtk)\n");
			//output per gnuplot (x:y:valore) compatibile con programmino passe_par_tout togliendo i #
			fprintf(clean_fields,"# 2D cut at x=%g\n", cutx[n]); 
			fprintf(clean_fields,"# %i\n#%i\n#%i\n",npunti_y, npunti_z, 1); 
			fprintf(clean_fields,"#%f %f\n#%f %f\n",ymin, zmin, ymax, zmax);
			size_t i = gridIndex_cutx[n];
			for(size_t k = 0; k < npunti_z; k++)
			{
				for(size_t j = 0; j < npunti_y; j++)
				{
					xx=parametri->ycoord[j];//xx=xmin+dx*i;
					yy=parametri->zcoord[k];//yy=ymin+dy*j;
					fprintf(clean_fields,"%.4g %.4g %.4g\n",xx, yy, field[i+j*npunti_x+k*npunti_x*npunti_y]);
				}
			}
			fclose(clean_fields);
		}
	}


	if (out_vtk)
	{
		printf("%lu\nScrittura vtk\n\n",(unsigned long) fread_size);

		float *x_coordinates, *y_coordinates, *z_coordinates;
		x_coordinates=new float[npunti_x];
		y_coordinates=new float[npunti_y];
		z_coordinates=new float[npunti_z];

		for (size_t i = 0; i < npunti_x; i++)
			x_coordinates[i] = parametri->xcoord[i];
		for (size_t i = 0; i < npunti_y; i++)
			y_coordinates[i] = parametri->ycoord[i];
		for (size_t i = 0; i < npunti_z; i++)
			z_coordinates[i] = parametri->zcoord[i];

		if(parametri->endian_machine == 0)
		{
			swap_endian_f(field,npunti_x*npunti_y*npunti_z);
			swap_endian_f(x_coordinates,npunti_x);
			swap_endian_f(y_coordinates,npunti_y);
			swap_endian_f(z_coordinates,npunti_z);
		}

		//////// DATASET UNSTRUCTURED_GRID VERSION    ////////

		sprintf(nomefile_campi,"%s.vtk",argv[1]);
		clean_fields=fopen(nomefile_campi, "w");
		printf("\nWriting the fields file\n");
		fprintf(clean_fields,"# vtk DataFile Version 2.0\n");
		fprintf(clean_fields,"titolo mio\n");
		fprintf(clean_fields,"BINARY\n");
		fprintf(clean_fields,"DATASET UNSTRUCTURED_GRID\n");
		fprintf(clean_fields,"POINTS %i float\n",npunti_x*((size_t)npunti_y)*npunti_z);
		float rr[3];
		for(size_t k=0;k<npunti_z;k++)
		{
			rr[2]=z_coordinates[k];
			for(size_t j=0;j<npunti_y;j++)
			{
				rr[1]=y_coordinates[j];
				for(size_t i=0;i<npunti_x;i++)
				{
					rr[0]=x_coordinates[i];
					fwrite((void*)rr,sizeof(float),3,clean_fields);
				}
			}
		}

		fprintf(clean_fields,"POINT_DATA %i\n",npunti_x*npunti_y*npunti_z);
		fprintf(clean_fields,"SCALARS %s float 1\n",parametri->support_label);
		fprintf(clean_fields,"LOOKUP_TABLE default\n");
		fwrite((void*)field,sizeof(float),npunti_x*((size_t)npunti_y)*npunti_z,clean_fields);
		fclose(clean_fields);
	}



	if (out_vtk_nostretch)
	{
		printf("%lu\nScrittura vtk della parte non stretchata della griglia\n\n",(unsigned long) fread_size);

		int inizio_punti_non_stretchati_x, inizio_punti_non_stretchati_y, inizio_punti_non_stretchati_z;
		int fine_punti_non_stretchati_x, fine_punti_non_stretchati_y, fine_punti_non_stretchati_z;
		size_t npunti_non_stretchati_x, npunti_non_stretchati_y, npunti_non_stretchati_z;
		if (parametri->stretched_grid)
		{
			inizio_punti_non_stretchati_x = 0;//(int) (npunti_x/6.0);
			inizio_punti_non_stretchati_y = (int) (npunti_y/6.0);
			inizio_punti_non_stretchati_z = (int) (npunti_z/6.0);
			fine_punti_non_stretchati_x = (int) (npunti_x);
			fine_punti_non_stretchati_y = (int) (npunti_y*5.0/6.0);
			fine_punti_non_stretchati_z = (int) (npunti_z*5.0/6.0);
			npunti_non_stretchati_x = (size_t) (fine_punti_non_stretchati_x - inizio_punti_non_stretchati_x);
			npunti_non_stretchati_y = (size_t) (fine_punti_non_stretchati_y - inizio_punti_non_stretchati_y);
			npunti_non_stretchati_z = (size_t) (fine_punti_non_stretchati_z - inizio_punti_non_stretchati_z);
			
		}
		else
		{
			fine_punti_non_stretchati_x = (int) npunti_x;
			npunti_non_stretchati_x = npunti_x;
			fine_punti_non_stretchati_y = (int) npunti_y;
			npunti_non_stretchati_y = npunti_y;
			fine_punti_non_stretchati_z = (int) npunti_z;
			npunti_non_stretchati_z = npunti_z;
			inizio_punti_non_stretchati_x = inizio_punti_non_stretchati_y = inizio_punti_non_stretchati_z = 0;
		}

		float * field_non_stretchato = new float[npunti_non_stretchati_x*npunti_non_stretchati_y*npunti_non_stretchati_z];
		int a = 0, b = 0, c = 0;

		for (int k = inizio_punti_non_stretchati_z; k < fine_punti_non_stretchati_z; k++)
		{
			c = k - inizio_punti_non_stretchati_z;
			for (int j = inizio_punti_non_stretchati_y; j < fine_punti_non_stretchati_y; j++)
			{
				b = j - inizio_punti_non_stretchati_y;
				for (int i = inizio_punti_non_stretchati_x; i < fine_punti_non_stretchati_x; i++)
				{
					a = i - inizio_punti_non_stretchati_x;
					field_non_stretchato[a+b*npunti_non_stretchati_x+c*npunti_non_stretchati_x*npunti_non_stretchati_y] = field[i+j*npunti_x+k*npunti_x*npunti_y];
				}
			}
		}

		if(parametri->endian_machine == 0) swap_endian_f(field_non_stretchato,npunti_non_stretchati_x*npunti_non_stretchati_y*npunti_non_stretchati_z);



		//////// DATASET STRUCTURED_POINTS VERSION    ////////

		float xmin_non_stretchato = parametri->xcoord[inizio_punti_non_stretchati_x];
		//		float xmax_non_stretchato = parametri->xcoord[fine_punti_non_stretchati_x];
		float ymin_non_stretchato = parametri->ycoord[inizio_punti_non_stretchati_y];
		//		float ymax_non_stretchato = parametri->ycoord[fine_punti_non_stretchati_y];
		float zmin_non_stretchato = parametri->zcoord[inizio_punti_non_stretchati_z];
		//		float zmax_non_stretchato = parametri->zcoord[fine_punti_non_stretchati_z];

		dx=parametri->xcoord[npunti_x/2]-parametri->xcoord[npunti_x/2-1];
		dy=parametri->ycoord[npunti_y/2]-parametri->ycoord[npunti_y/2-1];
		dz=parametri->zcoord[npunti_z/2]-parametri->zcoord[npunti_z/2-1];
		sprintf(nomefile_campi,"%s_nostretch.vtk",argv[1]);
		clean_fields=fopen(nomefile_campi, "wb");
		printf("\nWriting the fields file\n");
		fprintf(clean_fields,"# vtk DataFile Version 2.0\n");
		fprintf(clean_fields,"titolo mio\n");
		fprintf(clean_fields,"BINARY\n");
		fprintf(clean_fields,"DATASET STRUCTURED_POINTS\n");
		fprintf(clean_fields,"DIMENSIONS %i %i %i\n",npunti_non_stretchati_x, npunti_non_stretchati_y, npunti_non_stretchati_z);
		fprintf(clean_fields,"ORIGIN %f %f %f\n",xmin_non_stretchato, ymin_non_stretchato, zmin_non_stretchato);
		fprintf(clean_fields,"SPACING %f %f %f\n",dx, dy, dz);
		fprintf(clean_fields,"POINT_DATA %i\n",npunti_non_stretchati_x*npunti_non_stretchati_y*npunti_non_stretchati_z);
		fprintf(clean_fields,"SCALARS %s float 1\n",parametri->support_label);
		fprintf(clean_fields,"LOOKUP_TABLE default\n");
		fwrite((void*)field_non_stretchato,sizeof(float),npunti_non_stretchati_x*npunti_non_stretchati_y*npunti_non_stretchati_z,clean_fields);




		/******************************************************************************
		//////// DATASET RECTILINEAR_GRID VERSION    ////////

		float *x_coordinates, *y_coordinates, *z_coordinates;
		x_coordinates=new float[npunti_non_stretchati_x];
		y_coordinates=new float[npunti_non_stretchati_y];
		z_coordinates=new float[npunti_non_stretchati_z];
		for (int i = inizio_punti_non_stretchati_x; i < fine_punti_non_stretchati_x; i++) x_coordinates[i] = parametri->xcoord[i];
		for (int i = inizio_punti_non_stretchati_y; i < fine_punti_non_stretchati_y; i++) y_coordinates[i] = parametri->ycoord[i];
		for (int i = inizio_punti_non_stretchati_z; i < fine_punti_non_stretchati_z; i++) z_coordinates[i] = parametri->zcoord[i];
		if(parametri->endian_machine == 0)
		{
		swap_endian_f(x_coordinates,npunti_non_stretchati_x);
		swap_endian_f(y_coordinates,npunti_non_stretchati_y);
		swap_endian_f(z_coordinates,npunti_non_stretchati_z);
		}

		sprintf(nomefile_campi,"%s_out.vtk",argv[1]);
		clean_fields=fopen(nomefile_campi, "wb");
		printf("\nWriting the fields file\n");
		fprintf(clean_fields,"# vtk DataFile Version 2.0\n");
		fprintf(clean_fields,"titolo mio\n");
		fprintf(clean_fields,"BINARY\n");
		fprintf(clean_fields,"DATASET RECTILINEAR_GRID\n");
		fprintf(clean_fields,"DIMENSIONS %i %i %i\n",npunti_non_stretchati_x, npunti_non_stretchati_y, npunti_non_stretchati_z);
		fprintf(clean_fields,"X_COORDINATES %i float\n",npunti_non_stretchati_x);
		fwrite((void*)x_coordinates,sizeof(float),npunti_non_stretchati_x,clean_fields);
		fprintf(clean_fields,"Y_COORDINATES %i float\n",npunti_non_stretchati_y);
		fwrite((void*)y_coordinates,sizeof(float),npunti_non_stretchati_y,clean_fields);
		fprintf(clean_fields,"Z_COORDINATES %i float\n",npunti_non_stretchati_z);
		fwrite((void*)z_coordinates,sizeof(float),npunti_non_stretchati_z,clean_fields);
		fprintf(clean_fields,"POINT_DATA %i\n",npunti_non_stretchati_x*npunti_non_stretchati_y*npunti_non_stretchati_z);
		fprintf(clean_fields,"SCALARS %s float 1\n",parametri->support_label);
		fprintf(clean_fields,"LOOKUP_TABLE default\n");
		fwrite((void*)field_non_stretchato,sizeof(float),npunti_non_stretchati_x*npunti_non_stretchati_y*npunti_non_stretchati_z,clean_fields);
		******************************************************************************/

		fclose(clean_fields);
	}

	if (out_parameters)
	{
		sprintf(nomefile_parametri,"%s.parameters",argv[1]);
		parameters=fopen(nomefile_parametri, "w");
		printf("\nWriting parameters to file\n");
		fprintf(parameters,"interi\n");
		fprintf(parameters,"npe_y=%i\n",npe_y);
		fprintf(parameters,"npe_z=%i\n",npe_z);
		fprintf(parameters,"npe=%i\n",npe);
		fprintf(parameters,"nx=%i\n",nx);
		fprintf(parameters,"nx1=%i\n",nx1);
		fprintf(parameters,"ny1=%i\n",ny1);
		fprintf(parameters,"nyloc=%i\n",nyloc);
		fprintf(parameters,"nz1=%i\n",nz1);
		fprintf(parameters,"nzloc=%i\n",nzloc);
		fprintf(parameters,"ibx=%i\n",ibx);
		fprintf(parameters,"iby=%i\n",iby);
		fprintf(parameters,"ibz=%i\n",ibz);
		fprintf(parameters,"model=%i\n",model);
		fprintf(parameters,"dmodel=%i\n",dmodel);
		fprintf(parameters,"nsp=%i\n",nsp);
		fprintf(parameters,"np_loc=%i\n",np_loc);
		fprintf(parameters,"lpord=%i\n",lpord);
		fprintf(parameters,"deord=%i\n",deord);
		fprintf(parameters,"fvar=%i\n",fvar);
		fprintf(parameters,"========= fine interi\n");
		fprintf(parameters,"\n floating\n");
		fprintf(parameters,"tnow=%f\n",tnow);
		fprintf(parameters,"xmin=%f\n",xmin);
		fprintf(parameters,"xmax=%f\n",xmax);
		fprintf(parameters,"ymin=%f\n",ymin);
		fprintf(parameters,"ymax=%f\n",ymax);
		fprintf(parameters,"zmin=%f\n",zmin);
		fprintf(parameters,"zmax=%f\n",zmax);
		fprintf(parameters,"w0x=%f\n",w0x);
		fprintf(parameters,"w0y=%f\n",w0y);
		fprintf(parameters,"nrat=%f\n",nrat);
		fprintf(parameters,"a0=%f\n",a0);
		fprintf(parameters,"lam0=%f\n",lam0);
		fprintf(parameters,"E0=%f\n",E0);
		fprintf(parameters,"B0=%f\n",B0);
		fprintf(parameters,"ompe=%f\n",ompe);
		fprintf(parameters,"xt_in=%f\n",xt_in);
		fprintf(parameters,"xt_end=%f\n",xt_end);
		fprintf(parameters,"charge=%f\n",charge);
		fprintf(parameters,"mass=%f\n",mass);

		fprintf(parameters,"\n\nGrid along x axis\n");
		for (int i = 0; i < nx1; i++)
		{
			fprintf(parameters,"%.4g  ",parametri->xcoord[i]);
			if (i > 0 && i % 10 == 0) fprintf(parameters,"\n");
		}

		fprintf(parameters,"\n\nGrid along y axis\n");
		for (int i = 0; i < ny1; i++)
		{
			fprintf(parameters,"%.4g  ",parametri->ycoord[i]);
			if (i > 0 && i % 10 == 0) fprintf(parameters,"\n");
		}

		fprintf(parameters,"\n\nGrid along z axis\n");
		for (int i = 0; i < nz1; i++)
		{
			fprintf(parameters,"%.4g  ",parametri->zcoord[i]);
			if (i > 0 && i % 10 == 0) fprintf(parameters,"\n");
		}

		fclose(parameters);
	}

	printf("%lu\nFine\n\n",(unsigned long) fread_size);


	return 0;

}

#endif
