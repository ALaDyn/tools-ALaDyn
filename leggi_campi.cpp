#ifndef __LEGGI_CAMPI_C
#define __LEGGI_CAMPI_C

#include "leggi_binario_ALaDyn_fortran.h"


int leggi_campi(int argc, const char** argv, Parametri * parametri)
{
	std::ostringstream nomefile_bin;
	nomefile_bin << std::string(argv[1]) << ".bin";

	const int out_swap = parametri->p[SWAP];
	const int out_parameters = parametri->p[OUT_PARAMS];
	const int out_vtk = parametri->p[OUT_VTK];
	const int out_cutx = parametri->p[OUT_CUTX];
	const int out_cuty = parametri->p[OUT_CUTY];
	const int out_cutz = parametri->p[OUT_CUTZ];
	const int out_lineoutx = parametri->p[OUT_LINEOUT_X];
	const int out_2d = parametri->p[OUT_GRID2D];

	int npunti_x = parametri->npunti_x_ricampionati;
	int npunti_y = parametri->npunti_y_ricampionati;
	int npunti_z = parametri->npunti_z_ricampionati;
	int npunti_x_per_cpu = parametri->npx_per_cpu;
	int npunti_y_per_cpu = parametri->npy_per_cpu;
	int npunti_z_per_cpu = parametri->npz_per_cpu;
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
	std::vector<int> gridIndex_cutx, gridIndex_cuty, gridIndex_cutz;


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
	float *x_coordinates, *y_coordinates, *z_coordinates;
	double *x_lineout;
	char nomefile_parametri[MAX_LENGTH_FILENAME];
	char nomefile_campi[MAX_LENGTH_FILENAME];

	std::FILE *file_in;
	std::FILE *parameters;
	std::FILE *clean_fields;

	file_in=fopen(nomefile_bin.str().c_str(), "rb");

	size_t fread_size = 0;

	fread_size = std::fread(&buff,sizeof(int),1,file_in);
	fread_size = std::fread(&N_param,sizeof(int),1,file_in);
	fread_size = std::fread(&buff,sizeof(int),1,file_in);

	if(out_swap) swap_endian_i(&N_param,1);

	printf("numero parametri=%i\n",N_param);
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

	npe_y=int_param[0];     //numero processori lungo y
	npe_z=int_param[1];     //numero processori lungo z
	npe=npe_y*npe_z;     //numero processori
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

	if (out_parameters)
	{
		sprintf(nomefile_parametri,"%s.parameters",argv[1]);
		parameters=fopen(nomefile_parametri, "w");
		printf("\nWriting the parameters file\n");
		fprintf(parameters,"interi\n");
		fprintf(parameters,"npe_y=%i\n",npe_y);     //numero processori
		fprintf(parameters,"npe_z=%i\n",npe_z);     //numero processori
		fprintf(parameters,"npe=%i\n",npe);     //numero processori
		fprintf(parameters,"nx=%i\n",nx);
		fprintf(parameters,"nx1=%i\n",nx1);
		fprintf(parameters,"ny1=%i\n",ny1);
		fprintf(parameters,"nyloc=%i\n",nyloc);
		fprintf(parameters,"nz1=%i\n",nz1);
		fprintf(parameters,"nzloc=%i\n",nzloc);
		fprintf(parameters,"ibx=%i\n",ibx);
		fprintf(parameters,"iby=%i\n",iby);
		fprintf(parameters,"ibz=%i\n",ibz);
		fprintf(parameters,"model=%i\n",model);  //modello di laser utilizzato
		fprintf(parameters,"dmodel=%i\n",dmodel); //modello di condizioni iniziali
		fprintf(parameters,"nsp=%i\n",nsp);    //numero di speci
		fprintf(parameters,"np_loc=%i\n",np_loc);  //numero di componenti dello spazio dei momenti
		fprintf(parameters,"lpord=%i\n",lpord); //ordine dello schema leapfrog
		fprintf(parameters,"deord=%i\n",deord); //ordine derivate
		fprintf(parameters,"fvar=%i\n",fvar); 
		fprintf(parameters,"========= fine interi\n");
		fprintf(parameters,"\n floating\n");
		fprintf(parameters,"tnow=%f\n",tnow);  //tempo dell'output
		fprintf(parameters,"xmin=%f\n",xmin);  //estremi della griglia
		fprintf(parameters,"xmax=%f\n",xmax);  //estremi della griglia
		fprintf(parameters,"ymin=%f\n",ymin);  //estremi della griglia
		fprintf(parameters,"ymax=%f\n",ymax);  //estremi della griglia
		fprintf(parameters,"zmin=%f\n",zmin);  //estremi della griglia
		fprintf(parameters,"zmax=%f\n",zmax);  //estremi della griglia
		fprintf(parameters,"w0x=%f\n",w0x);      //waist del laser in x
		fprintf(parameters,"w0y=%f\n",w0y);      //waist del laser in y
		fprintf(parameters,"nrat=%f\n",nrat);     //n orver n critical
		fprintf(parameters,"a0=%f\n",a0);      // a0 laser
		fprintf(parameters,"lam0=%f\n",lam0);    // lambda
		fprintf(parameters,"E0=%f\n",E0);      //conversione da campi numerici a TV/m
		fprintf(parameters,"B0=%f\n",B0);
		fprintf(parameters,"ompe=%f\n",ompe);    //costante accoppiamento correnti campi
		fprintf(parameters,"xt_in=%f\n",xt_in);   //inizio plasma
		fprintf(parameters,"xt_end=%f\n",xt_end);
		fprintf(parameters,"charge=%f\n",charge);  //carica particella su carica elettrone
		fprintf(parameters,"mass=%f\n",mass);    //massa particelle su massa elettrone
		fclose(parameters);
	}


	printf("=========INIZIO LETTURE==========\n");
	printf("nx1*ny1*nz1: %i %i %i = %i\n",nx1,ny1,nz1,nx1*ny1*nz1);
	fflush(stdout);


	if (parametri->old_fortran_bin)
	{
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
	}


	field = new float[npunti_x*npunti_y*npunti_z];
	// field=(float*)malloc(npunti_x*npunti_y*npunti_z*sizeof(float));



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

			printf("processore ipz=%i/%i  ipy=%i/%i     segnoz=%i     segnoy=%i\r",ipz,npe_z, ipy,npe_y,segnoz,segnoy );
			fflush(stdout);

			buffer = new float[loc_size];
			//	buffer=(float*)malloc(loc_size*sizeof(float));
			fread_size = std::fread(&buff,sizeof(int),1,file_in);
			fread_size = std::fread(buffer,sizeof(float),loc_size,file_in);
			fread_size = std::fread(&buff,sizeof(int),1,file_in);

			if(out_swap) swap_endian_f(buffer,loc_size);

			for(int k=0; k<nzloc; k++)
				for(int j=0; j<nyloc; j++)
					for(int i=0; i<nxloc; i++)
						field[i+(j+segnoy)*nx1+(k+segnoz)*nx1*ny1]=buffer[i+j*nxloc+k*nxloc*nyloc];
			segnoy += nyloc;
		} 
		segnoz += nzloc;
	}



	dx=(xmax-xmin)/(npunti_x-1);
	dy=(ymax-ymin)/(npunti_y-1);
	if(npunti_z == 1)
		dz = (zmax-zmin);
	else
		dz = (zmax-zmin) / (npunti_z-1);

	x_coordinates=new float[npunti_x];
	y_coordinates=new float[npunti_y];
	z_coordinates=new float[npunti_z];
	x_lineout=new double[npunti_x];



	// leggiamo ora le coordinate dei punti di griglia, presenti solo nelle nuove versioni che possono prevedere griglia stretchata
	// se non ci sono, le posizioni vengono create, sotto, nell'"else", come struttura regolare a partire dagli indici di griglia e dai dx, dy, dz

	fread_size = std::fread(&buff,sizeof(int),1,file_in);	// facciamo il test sul buffer Fortran della prima coordinata, cosi' se c'e' non dobbiamo nemmeno tornare indietro
	// in futuro, con un output pulito, andra' modificata, facendo il test di lettura della prima coordinata ed eventualmente
	// tornando indietro se c'e' per poi entrare correttamente nella lettura seguente

	if(!std::feof(file_in))
	{
		fread_size = std::fread(x_coordinates,sizeof(float),nx1,file_in);
		fread_size = std::fread(&buff,sizeof(int),1,file_in);
		fread_size = std::fread(&buff,sizeof(int),1,file_in);
		fread_size = std::fread(y_coordinates,sizeof(float),ny1,file_in);
		fread_size = std::fread(&buff,sizeof(int),1,file_in);
		fread_size = std::fread(&buff,sizeof(int),1,file_in);
		fread_size = std::fread(z_coordinates,sizeof(float),nz1,file_in);
		fread_size = std::fread(&buff,sizeof(int),1,file_in);

		if(out_swap){
			swap_endian_f(x_coordinates,nx1);
			swap_endian_f(y_coordinates,ny1);
			swap_endian_f(z_coordinates,nz1);
		}

	}
	else
	{
		for(int i = 0; i < nx1; i++)
			x_coordinates[i] = xmin + dx*i;
		for(int i = 0; i < ny1; i++)
			y_coordinates[i] = ymin + dy*i;
		for(int i = 0; i < nz1; i++)
			z_coordinates[i] = zmin + dz*i;
	}

	printf("=========FINE LETTURE==========\n");
	fflush(stdout);



	if(npunti_z == 1 && out_2d)
	{
		printf("\nScrittura file gnuplot 2D\n");
		sprintf(nomefile_campi,"%s.txt",argv[1]);
		clean_fields=fopen(nomefile_campi, "w");
		printf("\nWriting the fields file 2D (not vtk)\n");

		//output per gnuplot (x:y:valore) compatibile con programmino passe_par_tout togliendo i #
		fprintf(clean_fields,"# %i\n#%i\n#%i\n",nx1, ny1, 1); 
		fprintf(clean_fields,"#%f %f\n#%f %f\n",xmin, ymin, xmax, ymax);
		for(int j = 0; j < ny1; j++)
		{
			for(int i = 0; i < nx1; i++)
			{
				xx=x_coordinates[i];//xmin+dx*i;
				yy=y_coordinates[j];//ymin+dy*j;
				fprintf(clean_fields,"%.4g %.4g %.4g\n",xx, yy, field[i+j*nx1]);
			}
		}
		fclose(clean_fields);
	}


	//avg of the 10 points aorund the focal axis
	if(npunti_z == 1 && out_lineoutx)
	{
		int myj;
		printf("\nScrittura lineout 1D\n");
		sprintf(nomefile_campi,"lineout_%s.txt",argv[1]);
		clean_fields=fopen(nomefile_campi, "w");
		printf("\nWriting the lineout file 1D (not vtk)\n");

		for(int i = 0; i < nx1; i++) x_lineout[i]=0;

		for(int j = 0; j < ny1; j++)
		{
			if(y_coordinates[j]>=0)
			{
				myj=j;
				break;
			}
		}

		fprintf(clean_fields,"#");
		for(int j = myj-span; j < (myj+span+1); j++)
		{
			fprintf(clean_fields,"%.4g\t",y_coordinates[j]);
			for(int i = 0; i < nx1; i++)
			{
				x_lineout[i]+=field[i+j*nx1]/(2*span+1.0);
			}
		}
		fprintf(clean_fields,"\n");
		for(int i = 0; i < nx1; i++)
		{
			xx=x_coordinates[i];//xmin+dx*i;
			fprintf(clean_fields,"%.4g %.4g\n",xx, x_lineout[i]);
		}
		fclose(clean_fields);
	}

	if(npunti_z > 1 && out_lineoutx)
	{
		int myj;
		printf("\nScrittura lineout 1D\n");
		sprintf(nomefile_campi,"lineout_%s.txt",argv[1]);
		clean_fields=fopen(nomefile_campi, "w");
		printf("\nWriting the lineout file 1D (not vtk)\n");

		for(int i = 0; i < nx1; i++)
			x_lineout[i]=0;
		for(int j = 0; j < ny1; j++)
		{
			if(y_coordinates[j]>=0)
			{
				myj=j;
				break;
			}
		}
		fprintf(clean_fields,"#");
		for(int k = myj-span; k < (myj+span+1); k++)
		{
			fprintf(clean_fields,"%.4g\t",z_coordinates[k]);

			for(int j = myj-span; j < (myj+span+1); j++)
			{
				for(int i = 0; i < nx1; i++)
				{
					x_lineout[i]+=field[i+j*nx1+k*nx1*ny1]/((2*span+1.0)*(2*span+1.0));
				}
			}
		}
		fprintf(clean_fields,"\n");
		for(int i = 0; i < nx1; i++)
		{
			xx=x_coordinates[i];//xmin+dx*i;
			fprintf(clean_fields,"%.4g %.4g\n",xx, x_lineout[i]);
		}
		fclose(clean_fields);
	}


	if (npunti_z > 1 && out_cutz)
	{
		if (parametri->posizioni_taglio_griglia_z.size() == 0)
		{
			taglio = z_coordinates[nz1/2];
			cutz.push_back(taglio);
			gridIndex_cutz.push_back(nz1/2);
		}
		else
		{
			taglio = 0.0;
			int i = 0;
			for (size_t j = 0; j < parametri->posizioni_taglio_griglia_z.size(); j++)
			{
				while (taglio < parametri->posizioni_taglio_griglia_z.at(j)) taglio = z_coordinates[i], i++;
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
			fprintf(clean_fields,"# %i\n#%i\n#%i\n",nx1, ny1, 1); 
			fprintf(clean_fields,"#%f %f\n#%f %f\n",xmin, ymin, xmax, ymax);
			int k = gridIndex_cutz[n];
			for(int j = 0; j < ny1; j++)
			{
				for(int i = 0; i < nx1; i++)
				{
					xx=x_coordinates[i];//xx=xmin+dx*i;
					yy=y_coordinates[j];//yy=ymin+dy*j;
					fprintf(clean_fields,"%.4g %.4g %.4g\n",xx, yy, field[i+j*nx1+k*nx1*ny1]);
				}
			}
			fclose(clean_fields);
		}
	}

	if (npunti_z > 1 && out_cuty)
	{
		if (parametri->posizioni_taglio_griglia_y.size() == 0)
		{
			taglio = y_coordinates[ny1/2];
			cuty.push_back(taglio);
			gridIndex_cuty.push_back(ny1/2);
		}
		else
		{
			taglio = 0.0;
			int i = 0;
			for (size_t j = 0; j < parametri->posizioni_taglio_griglia_y.size(); j++)
			{
				while (taglio < parametri->posizioni_taglio_griglia_y.at(j)) taglio = y_coordinates[i], i++;
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
			fprintf(clean_fields,"# %i\n#%i\n#%i\n",nx1, nz1, 1); 
			fprintf(clean_fields,"#%f %f\n#%f %f\n",xmin, zmin, xmax, zmax);
			int j = gridIndex_cuty[n];
			for(int k = 0; k < nz1; k++)
			{
				for(int i = 0; i < nx1; i++)
				{
					xx=x_coordinates[i];//xx=xmin+dx*i;
					yy=z_coordinates[k];//yy=ymin+dy*j;
					fprintf(clean_fields,"%.4g %.4g %.4g\n",xx, yy, field[i+j*nx1+k*nx1*ny1]);
				}
			}
			fclose(clean_fields);
		}
	}

	if (npunti_z > 1 && out_cutx)
	{
		if (parametri->posizioni_taglio_griglia_x.size() == 0)
		{
			taglio = x_coordinates[nx1/2];
			cutx.push_back(taglio);
			gridIndex_cutx.push_back(nx1/2);
		}
		else
		{
			taglio = 0.0;
			int i = 0;
			for (size_t j = 0; j < parametri->posizioni_taglio_griglia_x.size(); j++)
			{
				while (taglio < parametri->posizioni_taglio_griglia_x.at(j)) taglio = x_coordinates[i], i++;
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
			fprintf(clean_fields,"# %i\n#%i\n#%i\n",ny1, nz1, 1); 
			fprintf(clean_fields,"#%f %f\n#%f %f\n",ymin, zmin, ymax, zmax);
			int i = gridIndex_cutx[n];
			for(int k = 0; k < nz1; k++)
			{
				for(int j = 0; j < ny1; j++)
				{
					xx=y_coordinates[j];//xx=xmin+dx*i;
					yy=z_coordinates[k];//yy=ymin+dy*j;
					fprintf(clean_fields,"%.4g %.4g %.4g\n",xx, yy, field[i+j*nx1+k*nx1*ny1]);
				}
			}
			fclose(clean_fields);
		}
	}


	if (out_vtk)
	{
		printf("%lu\nScrittura vtk\n\n",(unsigned long) fread_size);

		if(parametri->endian_machine == 0)
		{
			swap_endian_f(field,nx1*ny1*nz1);
			swap_endian_f(x_coordinates,nx1);
			swap_endian_f(y_coordinates,ny1);
			swap_endian_f(z_coordinates,nz1);
		}

		//////// DATASET STRUCTURED_POINTS VERSION    ////////

		sprintf(nomefile_campi,"%s_out.vtk",argv[1]);
		clean_fields=fopen(nomefile_campi, "w");
		printf("\nWriting the fields file\n");
		fprintf(clean_fields,"# vtk DataFile Version 2.0\n");
		fprintf(clean_fields,"titolo mio\n");
		fprintf(clean_fields,"BINARY\n");
		//fprintf(clean_fields,"DATASET STRUCTURED_POINTS\n");
		fprintf(clean_fields,"DATASET UNSTRUCTURED_GRID\n");
		fprintf(clean_fields,"POINTS %i float\n",nx1*ny1*nz1);
		float rr[3];
		for(int k=0;k<nz1;k++)
		{
			rr[2]=z_coordinates[k];
			for(int j=0;j<ny1;j++)
			{
				rr[1]=y_coordinates[j];
				for(int i=0;i<nx1;i++)
				{
					rr[0]=x_coordinates[i];
					fwrite((void*)rr,sizeof(float),3,clean_fields);
				}
			}
		}

		fprintf(clean_fields,"POINT_DATA %i\n",nx1*ny1*nz1);
		fprintf(clean_fields,"SCALARS %s float 1\n",parametri->support_label);
		fprintf(clean_fields,"LOOKUP_TABLE default\n");
		fwrite((void*)field,sizeof(float),nx1*ny1*nz1,clean_fields);
		fclose(clean_fields);


		/*
		sprintf(nomefile_campi,"%s_out.vtk",argv[1]);
		clean_fields=fopen(nomefile_campi, "wb");
		printf("\nWriting the fields file\n");
		fprintf(clean_fields,"# vtk DataFile Version 2.0\n");
		fprintf(clean_fields,"titolo mio\n");
		fprintf(clean_fields,"BINARY\n");
		//fprintf(clean_fields,"DATASET STRUCTURED_POINTS\n");
		fprintf(clean_fields,"DATASET RECTILINEAR_GRID\n");
		fprintf(clean_fields,"DIMENSIONS %i %i %i\n",nx1, ny1, nz1);
		//fprintf(clean_fields,"ORIGIN %f %f %f\n",xmin, ymin, zmin);
		//fprintf(clean_fields,"SPACING %f %f %f\n",dx, dy, dz);
		fprintf(clean_fields,"X_COORDINATES %i float\n",nx1);
		fwrite((void*)x_coordinates,sizeof(float),nx1,clean_fields);
		fprintf(clean_fields,"Y_COORDINATES %i float\n",ny1);
		fwrite((void*)y_coordinates,sizeof(float),ny1,clean_fields);
		fprintf(clean_fields,"Z_COORDINATES %i float\n",nz1);
		fwrite((void*)z_coordinates,sizeof(float),nz1,clean_fields);

		fprintf(clean_fields,"POINT_DATA %i\n",nx1*ny1*nz1);
		fprintf(clean_fields,"SCALARS %s float 1\n",parametri->support_label);
		fprintf(clean_fields,"LOOKUP_TABLE default\n");
		fwrite((void*)field,sizeof(float),nx1*ny1*nz1,clean_fields);
		fclose(clean_fields);
		*/
		fclose(file_in);
	}


	printf("%lu\nFine\n\n",(unsigned long) fread_size);


	return 0;

}

#endif
