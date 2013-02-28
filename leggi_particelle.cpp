#ifndef __LEGGI_PARTICELLE_C
#define __LEGGI_PARTICELLE_C

#include "leggi_binario_ALaDyn_fortran.h"


int leggi_particelle(int argc, const char ** argv, parametri * parametri)
{
	std::ostringstream nomefile_bin, nomefile_dat, nomefile_xpx, nomefile_Etheta, nomefile_Espec, nomefile_Estremi;
	nomefile_bin << std::string(argv[1]) << ".bin";
	nomefile_dat << std::string(argv[1]) << ".dat";
	std::string riga_persa;
	char* trascura;
	trascura = new char[MAX_LENGTH_FILENAME];
	FILE *file_in;
	if ( (file_in=fopen(nomefile_bin.str().c_str(), "r")) == NULL )
	{
		printf ( "file non-existant!\n" );
	}
	else
	{
		printf ( "file exists!\n" );
	}
	FILE *binary_all_out;
	FILE *ascii_all_out;
	FILE *parameters;
	std::ofstream xpx_out, Etheta_out, Espec_out, Estremi_out;
	std::ifstream file_dat;
	if (!(parametri->old_fortran_bin)) 
	{
		file_dat.open(nomefile_dat.str().c_str());
		if (file_dat.fail()) printf ( "file .dat non-existant!\n" );
	}


	int out_swap = parametri->p[SWAP];
	int out_binary = parametri->p[OUT_BINARY];
	int out_ascii = parametri->p[OUT_ASCII];
	int out_parameters = parametri->p[OUT_PARAMS];
	int fai_binning = parametri->p[DO_BINNING];
	int cerca_minmax = parametri->p[FIND_MINMAX];

	int ipc,N_param, *int_param,npart_loc;
	int buff, pID;

	int nelab=3;	//3 valori per ora: gamma, theta ed energia
	float x,y,z,px,py,pz;
	float *estremi_min, *estremi_max;

	short buffshort[2];
	float *particelle, *real_param;
	float gamma, theta, E;
	char nomefile_binary[MAX_LENGTH_FILENAME];
	char nomefile_ascii[MAX_LENGTH_FILENAME];
	char nomefile_parametri[MAX_LENGTH_FILENAME];

	size_t fread_size;

	int npe,nx,nz,ibx,iby,ibz,model,dmodel,nsp,ndim,lpord,deord,nptot, ny_loc, np_loc,ndv, i_end;
	//	int ny;
	float tnow,xmin,xmax,ymin,ymax,zmin,zmax,w0x,w0y,nrat,a0,lam0,E0,ompe,xt_in,xt_end,charge,mass, np_over_nm;
	float rx, ry, rz, ux, uy, uz, wgh;
	int tipo = 0;
	if (argv[1][0] == 'E') tipo = 3;
	else if (argv[1][0] == 'P') tipo = 1;
	else printf("Tipo non riconosciuto!\n");

	if (parametri->old_fortran_bin) 
	{
		fread_size = std::fread((void*) &buff,sizeof(int),1,file_in);
		fread_size = std::fread((void*) &N_param,sizeof(int),1,file_in);
		fread_size = std::fread((void*) &buff,sizeof(int),1,file_in);
		if (out_swap) swap_endian_i(&N_param,1);
		printf("numero parametri %i\n",N_param);
		fflush(stdout);
		int_param=(int*)malloc(N_param*sizeof(int));
		real_param=(float*)malloc(N_param*sizeof(float));
		fread_size = std::fread(&buff,sizeof(int),1,file_in);
		fread_size = std::fread(int_param,sizeof(int),N_param,file_in);
		fread_size = std::fread(&buff,sizeof(int),1,file_in);
		fread_size = std::fread(&buff,sizeof(int),1,file_in);
		fread_size = std::fread(real_param,sizeof(float),N_param,file_in);
		fread_size = std::fread(&buff,sizeof(int),1,file_in);
		if (out_swap) swap_endian_i(int_param,N_param);
		if (out_swap) swap_endian_f(real_param,N_param);
		npe=int_param[0];     //numero processori
		nx=int_param[1];
		ny_loc=int_param[2];
		nz=int_param[3];
		ibx=int_param[4];
		iby=int_param[5];
		ibz=int_param[6];
		model=int_param[7];  //modello di laser utilizzato
		dmodel=int_param[8]; //modello di condizioni iniziali
		nsp=int_param[9];    //numero di speci
		ndim=int_param[10];  //numero di componenti dello spazio dei momenti
		np_loc=int_param[11];  //
		lpord=int_param[12]; //ordine dello schema leapfrog
		deord=int_param[13]; //ordine derivate
		//current type
		pID=int_param[15];  //tipo delle particelle
		nptot=int_param[16]; //numero di particelle da leggere nella specie
		ndv=int_param[17]; //numero di particelle da leggere nella specie
		tnow=real_param[0];  /*tempo dell'output*/
		xmin=real_param[1];  //estremi della griglia
		xmax=real_param[2];  //estremi della griglia
		ymin=real_param[3];  //estremi della griglia
		ymax=real_param[4];  //estremi della griglia
		zmin=real_param[5];  //estremi della griglia
		zmax=real_param[6];  //estremi della griglia
		w0x=real_param[7];      //waist del laser in x
		w0y=real_param[8];      //waist del laser in y
		nrat=real_param[9];     //n over n critical
		a0=real_param[10];      // a0 laser
		lam0=real_param[11];    // lambda
		E0=real_param[12];      //conversione da campi numerici a TV/m
		ompe=real_param[13];    //costante accoppiamento correnti campi
		np_over_nm=real_param[14];		//numerical2physical particles 14 
		xt_in=real_param[15];			//inizio plasma
		xt_end=real_param[16];
		charge=real_param[17];			//carica particella su carica elettrone
		mass=real_param[18];			//massa particelle su massa elettrone
	}
	else
	{
		std::getline(file_dat,riga_persa);
		file_dat >> npe;
		file_dat >> nx;
		file_dat >> ny_loc;
		file_dat >> nz;
		file_dat >> ibx;
		file_dat >> iby;
		file_dat >> ibz;
		file_dat >> model;
		file_dat >> dmodel;
		file_dat >> nsp;
		file_dat >> ndim;
		file_dat >> np_loc;
		file_dat >> lpord;
		file_dat >> deord;
		file_dat >> trascura;
		file_dat >> pID;
		file_dat >> nptot;
		file_dat >> ndv;
		file_dat >> i_end;
		file_dat >> trascura;
		std::getline(file_dat,riga_persa); // nb: servono due getline per eliminare la riga di testo nel .dat e non capisco il perche'...
		std::getline(file_dat,riga_persa);
		file_dat >> tnow;
		file_dat >> xmin;
		file_dat >> xmax;
		file_dat >> ymin;
		file_dat >> ymax;
		file_dat >> zmin;
		file_dat >> zmax;
		file_dat >> w0x;
		file_dat >> w0y;
		file_dat >> nrat;
		file_dat >> a0;
		file_dat >> lam0;
		file_dat >> E0;
		file_dat >> ompe;
		file_dat >> np_over_nm;
		file_dat >> xt_in;
		file_dat >> xt_end;
		file_dat >> charge;
		file_dat >> mass;
		file_dat >> trascura;
	}



	estremi_min = new float[2*ndim+nelab];
	estremi_max = new float[2*ndim+nelab];
	for (int i = 0; i < (2*ndim+nelab); i++)
	{
		estremi_min[i] = (float) NUMERO_MASSIMO;
		estremi_max[i] = (float) -NUMERO_MASSIMO;
	}

	//	ny=ny_loc*npe;
	printf("nptot=%i\n",nptot); 
	printf("\ninizio processori \n");
	fflush(stdout);
	float **xpx = new float* [parametri->nbin_x+3];
	for (int i = 0; i < parametri->nbin_x+3; i++)
	{
		xpx[i] = new float [parametri->nbin_px+3];
		for (int j = 0; j < parametri->nbin_px+3; j++) xpx[i][j] = 0.0;
	}
	float **Etheta = new float* [parametri->nbin_E+3];
	for (int i = 0; i < parametri->nbin_E+3; i++)
	{
		Etheta[i] = new float [parametri->nbin_theta+3];
		for (int j = 0; j < parametri->nbin_theta+3; j++) Etheta[i][j] = 0.0;
	}
	float *Espec = new float [parametri->nbin_E+3];

	sprintf(nomefile_ascii,"%s.ascii",argv[1]);
	sprintf(nomefile_binary,"%s.clean",argv[1]);


#ifdef ENABLE_DEBUG
	if (fai_binning)
	{
		std::cout << "XMIN = " << parametri->xmin << std::endl;
		std::cout << "XMAX = " << parametri->xmax << std::endl;
		std::cout << "YMIN = " << parametri->ymin << std::endl;
		std::cout << "YMAX = " << parametri->ymax << std::endl;
		std::cout << "ZMIN = " << parametri->zmin << std::endl;
		std::cout << "ZMAX = " << parametri->zmax << std::endl;
		std::cout << "PXMIN = " << parametri->pxmin << std::endl;
		std::cout << "PXMAX = " << parametri->pxmax << std::endl;
		std::cout << "PYMIN = " << parametri->pymin << std::endl;
		std::cout << "PYMAX = " << parametri->pymax << std::endl;
		std::cout << "PZMIN = " << parametri->pzmin << std::endl;
		std::cout << "PZMAX = " << parametri->pzmax << std::endl;
		std::cout << "GAMMAMIN = " << parametri->gammamin << std::endl;
		std::cout << "GAMMAMAX = " << parametri->gammamax << std::endl;
		std::cout << "THETAMIN = " << parametri->thetamin << std::endl;
		std::cout << "THETAMAX = " << parametri->thetamax << std::endl;
		std::cout << "EMIN = " << parametri->Emin << std::endl;
		std::cout << "EMAX = " << parametri->Emax << std::endl;
	}
#endif


	for(ipc=0;ipc<npe;ipc++)
	{

		if (parametri->old_fortran_bin)
		{
			fread_size = std::fread(&buff,sizeof(int),1,file_in); 
			fread_size = std::fread(&npart_loc,sizeof(int),1,file_in);
			fread_size = std::fread(&buff,sizeof(int),1,file_in);
			if (out_swap) swap_endian_i(&npart_loc,1);
		}
		else fread_size = std::fread(&npart_loc,sizeof(int),1,file_in);

		particelle=(float*)malloc(npart_loc*(2*ndim+parametri->p[WEIGHT])*sizeof(float));
		printf("proc number %i\tnpart=%i\n",ipc,npart_loc);
		unsigned int val[] = {(unsigned int)npart_loc, (unsigned int)(2*ndim+parametri->p[WEIGHT])};
		if(npart_loc>0)
		{
			printf("\tentro\t");
			fflush(stdout);

			if (parametri->old_fortran_bin)
			{
				fread_size = std::fread(buffshort,sizeof(short),2,file_in);
				// buffshort[0] = ???
				// buffshort[1] = ???
				if (out_swap) swap_endian_s(buffshort,2);
				fread_size = std::fread(particelle,sizeof(float),npart_loc*ndv,file_in);
				fread_size = std::fread(&buff,sizeof(int),1,file_in);
				if (out_swap) swap_endian_f(particelle,npart_loc*ndv);
			}
			else fread_size = std::fread(particelle,sizeof(float),npart_loc*ndv,file_in);

			printf("lunghezza=%i    %hu\t%hu\n",npart_loc*(2*ndim+parametri->p[WEIGHT]),buffshort[0],buffshort[1]);
			_Filtro(parametri, particelle,val,_Filtro::costruisci_filtro(argc, argv));

			if (cerca_minmax)
			{
				for (int i = 0; i < npart_loc; i++)
				{
					if (ndv == 6 || ndv == 7)
					{
						x=*(particelle+i*ndv);
						y=*(particelle+i*ndv+1);
						z=*(particelle+i*ndv+2);
						px=*(particelle+i*ndv+3);
						py=*(particelle+i*ndv+4);
						pz=*(particelle+i*ndv+5);
						gamma=(float)(sqrt(1.+px*px+py*py+pz*pz)-1.);			//gamma
						theta=(float)(atan2(sqrt(py*py+pz*pz),px)*180./M_PI);	//theta nb: py e pz sono quelli trasversi in ALaDyn!
						E=(float)(gamma*parametri->massa_particella_MeV);		//energia
						if (x < estremi_min[0]) estremi_min[0] = x;
						if (x > estremi_max[0]) estremi_max[0] = x;
						if (y < estremi_min[1]) estremi_min[1] = y;
						if (y > estremi_max[1]) estremi_max[1] = y;
						if (z < estremi_min[2]) estremi_min[2] = z;
						if (z > estremi_max[2]) estremi_max[2] = z;
						if (px < estremi_min[3]) estremi_min[3] = px;
						if (px > estremi_max[3]) estremi_max[3] = px;
						if (py < estremi_min[4]) estremi_min[4] = py;
						if (py > estremi_max[4]) estremi_max[4] = py;
						if (pz < estremi_min[5]) estremi_min[5] = pz;
						if (pz > estremi_max[5]) estremi_max[5] = pz;
						if (gamma < estremi_min[6]) estremi_min[6] = gamma;
						if (gamma > estremi_max[6]) estremi_max[6] = gamma;
						if (theta < estremi_min[7]) estremi_min[7] = theta;
						if (theta > estremi_max[7]) estremi_max[7] = theta;
						if (E < estremi_min[8]) estremi_min[8] = E;
						if (E > estremi_max[8]) estremi_max[8] = E;
					}
					else if (ndv == 4 || ndv == 5)
					{
						x=*(particelle+i*ndv);
						y=*(particelle+i*ndv+1);
						px=*(particelle+i*ndv+2);
						py=*(particelle+i*ndv+3);
						gamma=(float)(sqrt(1.+px*px+py*py)-1.);				//gamma
						theta=(float)(atan2(py,px)*180./M_PI);				//theta
						E=(float)(gamma*parametri->massa_particella_MeV);	//energia
						if (x < estremi_min[0]) estremi_min[0] = x;
						if (x > estremi_max[0]) estremi_max[0] = x;
						if (y < estremi_min[1]) estremi_min[1] = y;
						if (y > estremi_max[1]) estremi_max[1] = y;
						if (px < estremi_min[3]) estremi_min[3] = px;
						if (px > estremi_max[3]) estremi_max[3] = px;
						if (py < estremi_min[4]) estremi_min[4] = py;
						if (py > estremi_max[4]) estremi_max[4] = py;
						if (gamma < estremi_min[6]) estremi_min[6] = gamma;
						if (gamma > estremi_max[6]) estremi_max[6] = gamma;
						if (theta < estremi_min[7]) estremi_min[7] = theta;
						if (theta > estremi_max[7]) estremi_max[7] = theta;
						if (E < estremi_min[8]) estremi_min[8] = E;
						if (E > estremi_max[8]) estremi_max[8] = E;
					}
				}
				_Binnaggio(particelle,npart_loc,ndv,parametri,xpx,"x","px");
				_Binnaggio(particelle,npart_loc,ndv,parametri,Espec,"E");
				_Binnaggio(particelle,npart_loc,ndv,parametri,Etheta,"E","theta");
			}

                        int contatori[] = {0,0,0};
			if (out_binary)
			{
				binary_all_out=fopen(nomefile_binary, "ab");
                                contatori[0] +=
                                fprintf(binary_all_out, "# vtk DataFile Version 2.0\n");
                                contatori[0] +=
                                fprintf(binary_all_out, "titolo nostro\n");
                                contatori[0] +=
                                fprintf(binary_all_out, "BINARY\n");
                                contatori[0] +=
                                fprintf(binary_all_out, "DATASET UNSTRUCTURED_GRID\n");
                                contatori[0] +=
                                fprintf(binary_all_out, "POINTS %i float\n", nptot);
				printf("\nWriting the C binary file\n");
                                float zero = 0.0f;
                                switch(ndv)
                                  {
                                   case 4:
                                   case 5:
                                    for(int i=0; i < nptot; i += ndv)
                                     fwrite((void*)(particelle+i),sizeof(float),2,binary_all_out),
                                     fwrite((void*)&zero, sizeof(float), 1, binary_all_out);
                                   break;
                                   case 6:
                                   case 7:
                                    for(int i=0; i < nptot; i += ndv)
                                     fwrite((void*)(particelle+i),sizeof(float),3,binary_all_out);
                                   break;
                                   default:
                                    printf("Sinigardi è pessimista %i\n", ndv);
                                   }
                                
                                contatori[1] +=
                                fprintf(binary_all_out, "DATASET UNSTRUCTURED_GRID\n");
                                contatori[1] +=
                                fprintf(binary_all_out, "POINTS %i float\n", nptot);
                                switch(ndv)
                                  {
                                   case 4:
                                   case 5:
                                    for(int i=2; i < nptot; i += ndv)
                                     fwrite((void*)(particelle+i),sizeof(float),2,binary_all_out),
                                     fwrite((void*)&zero, sizeof(float), 1, binary_all_out);
                                   break;
                                   case 6:
                                   case 7:
                                    for(int i=3; i < nptot; i += ndv)
                                     fwrite((void*)(particelle+i),sizeof(float),3,binary_all_out);
                                   break;
                                   default:
                                    printf("Sinigardi è pessimista %i\n", ndv);
                                   }
                                
				if(ndv == 5 || ndv == 7)
                                  {
                                   contatori[2] +=
	                           fprintf(binary_all_out,"DATASET STRUCTURED_POINTS\n");
                                   contatori[2] +=
	                           fprintf(binary_all_out,"DIMENSIONS %i %i %i\n",nptot, 1, 1);
                                   contatori[2] +=
	                           fprintf(binary_all_out,"ORIGIN 0 0 0\n");
                                   contatori[2] +=
	                           fprintf(binary_all_out,"SPACING 1 1 1\n");
                                   contatori[2] +=
	                           fprintf(binary_all_out,"POINT_DATA %i\n",nptot);
                                   contatori[2] +=
	                           fprintf(binary_all_out,"SCALARS w float 1\n");
                                   contatori[2] +=
	                           fprintf(binary_all_out,"LOOKUP_TABLE default\n");
                                   for(int i=ndv-1; i < nptot; i += ndv)
                                     fwrite((void*)(particelle+i),sizeof(float),1,binary_all_out);
                                  }

				fflush(binary_all_out);
				fclose(binary_all_out);
			}




			if(out_ascii)
			{
				ascii_all_out=fopen(nomefile_ascii, "a");

				if (ndv == 6 || ndv == 7)
				{
					for(int i=0;i<nptot;i++)
					{
						rx=particelle[i*ndv+1]*((float)1.e-4);
						ry=particelle[i*ndv+2]*((float)1.e-4);
						rz=particelle[i*ndv+0]*((float)1.e-4);
						ux=particelle[i*ndv+4];
						uy=particelle[i*ndv+5];
						uz=particelle[i*ndv+3];
						if (parametri->p[WEIGHT])
						{
							wgh=particelle[i*ndv+6];
							fprintf(ascii_all_out,"%e %e %e %e %e %e %d %e 0 %d\n",rx, ry, rz, ux, uy, uz, tipo, wgh, i+1);
						}
						else
						{
							fprintf(ascii_all_out,"%e %e %e %e %e %e %d 1 0 %d\n",rx, ry, rz, ux, uy, uz, tipo, i+1);
						}
					}
				}
				else if (ndv == 4 || ndv == 5)
				{
					for(int i=0;i<nptot;i++)
					{
						rx=particelle[i*ndv+1]*((float)1.e-4);
						rz=particelle[i*ndv+0]*((float)1.e-4);
						ux=particelle[i*ndv+3];
						uz=particelle[i*ndv+2];
						if (parametri->p[WEIGHT])
						{
							wgh=particelle[i*ndv+4];
							fprintf(ascii_all_out,"%e %e %e %e %d %e 0 %d\n",rx, rz, ux, uz, tipo, wgh, i+1);
						}
						else
						{
							fprintf(ascii_all_out,"%e %e %e %e %d 1 0 %d\n",rx, rz, ux, uz, tipo, i+1);
						}
					}
				}
				fflush(ascii_all_out);
				fclose(ascii_all_out);
			}

			free(particelle);
		}
	}


	if (out_parameters)
	{
		sprintf(nomefile_parametri,"%s.parameters",argv[1]);
		parameters=fopen(nomefile_parametri, "w");
		printf("\nWriting the parameters file\n");
		fprintf(parameters,"interi\n");
		fprintf(parameters,"npe=%i\n",npe);     //numero processori
		fprintf(parameters,"nx=%i\n",nx);
		fprintf(parameters,"ny=%i\n",ny_loc);
		fprintf(parameters,"nz=%i\n",nz);
		fprintf(parameters,"ibx=%i\n",ibx);
		fprintf(parameters,"iby=%i\n",iby);
		fprintf(parameters,"ibz=%i\n",ibz);
		fprintf(parameters,"model=%i\n",model);  //modello di laser utilizzato
		fprintf(parameters,"dmodel=%i\n",dmodel); //modello di condizioni iniziali
		fprintf(parameters,"nsp=%i\n",nsp);    //numero di speci
		fprintf(parameters,"ndim=%i\n",ndim);   
		fprintf(parameters,"np_loc=%i\n",np_loc);  
		fprintf(parameters,"lpord=%i\n",lpord); //ordine dello schema leapfrog
		fprintf(parameters,"deord=%i\n",deord); //ordine derivate
		fprintf(parameters,"pID=%i\n",pID); 
		fprintf(parameters,"nptot=%i\n",nptot); 
		fprintf(parameters,"ndv=%i\n",ndv); 
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
		fprintf(parameters,"nrat=%f\n",nrat);     //n over n critical
		fprintf(parameters,"a0=%f\n",a0);      // a0 laser
		fprintf(parameters,"lam0=%f\n",lam0);    // lambda
		fprintf(parameters,"E0=%f\n",E0);      //conversione da campi numerici a TV/m
		fprintf(parameters,"ompe=%f\n",ompe);    //costante accoppiamento correnti campi
		fprintf(parameters,"np_over_nm=%f\n",np_over_nm);   //numerical2physical particles 14 
		fprintf(parameters,"xt_in=%f\n",xt_in);
		fprintf(parameters,"xt_end=%f\n",xt_end);
		fprintf(parameters,"charge=%f\n",charge);  //carica particella su carica elettrone
		fprintf(parameters,"mass=%f\n",mass);    //massa particelle su massa elettrone
		//		if(WEIGHT) fprintf(parameters,"weight=%f\n",particelle[6]);    //massa particelle su massa elettrone
		fclose(parameters);
	}

	if (cerca_minmax)
	{
		nomefile_Estremi << argv[1] << ".extremes";
		Estremi_out.open(nomefile_Estremi.str().c_str());
		if (Estremi_out.fail()) printf("unable to create .extremes file");
		Estremi_out << "XMIN = " << estremi_min[0] << std::endl;
		Estremi_out << "XMAX = " << estremi_max[0] << std::endl;
		Estremi_out << "YMIN = " << estremi_min[1] << std::endl;
		Estremi_out << "YMAX = " << estremi_max[1] << std::endl;
		Estremi_out << "ZMIN = " << estremi_min[2] << std::endl;
		Estremi_out << "ZMAX = " << estremi_max[2] << std::endl;
		Estremi_out << "PXMIN = " << estremi_min[3] << std::endl;
		Estremi_out << "PXMAX = " << estremi_max[3] << std::endl;
		Estremi_out << "PYMIN = " << estremi_min[4] << std::endl;
		Estremi_out << "PYMAX = " << estremi_max[4] << std::endl;
		Estremi_out << "PZMIN = " << estremi_min[5] << std::endl;
		Estremi_out << "PZMAX = " << estremi_max[5] << std::endl;
		Estremi_out << "GAMMAMIN = " << estremi_min[6] << std::endl;
		Estremi_out << "GAMMAMAX = " << estremi_max[6] << std::endl;
		Estremi_out << "THETAMIN = " << estremi_min[7] << std::endl;
		Estremi_out << "THETAMAX = " << estremi_max[7] << std::endl;
		Estremi_out << "EMIN = " << estremi_min[8] << std::endl;
		Estremi_out << "EMAX = " << estremi_max[8] << std::endl;
		Estremi_out.close();
	}


	if (fai_binning)
	{
		nomefile_xpx << argv[1] << "_xpx";
		nomefile_Etheta << argv[1] << "_Etheta";
		nomefile_Espec << argv[1] << "_Espec";
		xpx_out.open(nomefile_xpx.str().c_str());
		Etheta_out.open(nomefile_Etheta.str().c_str());
		Espec_out.open(nomefile_Espec.str().c_str());

		float min1, min2, max1, max2;

		min1=parametri->Emin-parametri->dimmi_dimE();
		max1=parametri->Emin;
		for (int i = 0; i < parametri->nbin_E+3; i++)
		{
			Espec_out << std::setprecision(7) << min1 << "\t" << max1 << "\t" << Espec[i] << std::endl;

			min1 += parametri->dimmi_dimE();
			max1 += parametri->dimmi_dimE();
		}


		min1=parametri->Emin-parametri->dimmi_dimE();
		max1=parametri->Emin;
		min2=parametri->thetamin-parametri->dimmi_dimtheta();
		max2=parametri->thetamin;

		for (int i = 0; i < parametri->nbin_E+3; i++)
		{
			for (int j = 0; j < parametri->nbin_theta+3; j++)
			{
				Etheta_out << std::setprecision(7) << min1 << "\t" << max1 << "\t" << min2 << "\t" << max2 << "\t" << Etheta[i][j] << std::endl;
				min2 += parametri->dimmi_dimtheta();
				max2 += parametri->dimmi_dimtheta();
			}
			min1 += parametri->dimmi_dimE();
			max1 += parametri->dimmi_dimE();
			min2 = parametri->thetamin-parametri->dimmi_dimtheta();
			max2 = parametri->thetamin;

		}


		min1=parametri->xmin-parametri->dimmi_dimx();
		max1=parametri->xmin;
		min2=parametri->pxmin-parametri->dimmi_dimpx();
		max2=parametri->pxmin;

		for (int i = 0; i < parametri->nbin_x+3; i++)
		{
			for (int j = 0; j < parametri->nbin_px+3; j++)
			{
				xpx_out << std::setprecision(7) << min1 << "\t" << max1 << "\t" << min2 << "\t" << max2 << "\t" << xpx[i][j] << std::endl;
				min2 += parametri->dimmi_dimpx();
				max2 += parametri->dimmi_dimpx();
			}
			min1 += parametri->dimmi_dimx();
			max1 += parametri->dimmi_dimx();
			min2 = parametri->pxmin-parametri->dimmi_dimpx();
			max2 = parametri->pxmin;

		}

		xpx_out.close();
		Etheta_out.close();
		Espec_out.close();

	}

	printf("%lu\nFine\n\n",(unsigned long) fread_size);

	if (!(parametri->old_fortran_bin)) file_dat.close();
	fclose(file_in);
	return 0;
}


#endif
