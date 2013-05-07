#ifndef __LEGGI_CAMPI_C
#define __LEGGI_CAMPI_C

#include "leggi_binario_ALaDyn_fortran.h"


int leggi_campi(int argc, const char** argv, Parametri * parametri)
{
	std::ostringstream nomefile_bin;
	nomefile_bin << std::string(argv[1]) << ".bin";

	int out_swap = parametri->p[SWAP];
	//	int out_parameters = parametri.p[OUT_PARAMS];

	//	int kk, segnox=0;
	int i, ipy, ipz,  j, k, N_param, *int_param, np_loc, npoint_loc[3], loc_size;
	int segnoy=0,segnoz=0, buff;
	float *field,*buffer,  *real_param;
	float *x_coordinates, *y_coordinates, *z_coordinates;
	char nomefile_parametri[MAX_LENGTH_FILENAME];
	char nomefile_campi[MAX_LENGTH_FILENAME];

	FILE *file_in;
	FILE *parameters;
	FILE *clean_fields;

	int nx,ibx,iby,ibz,model,dmodel,nsp,lpord,deord,npe,fvar;
	int npe_y, npe_z;
	int nxloc, nx1, ny1, nyloc, nz1, nzloc;
	float tnow,w0x,w0y,nrat,a0,lam0,B0,ompe,xt_in,xt_end,charge,mass;
	float xmin,xmax,ymin,ymax,zmin,zmax,E0;
	float dx, dy, dz, xx, yy;
	//	float zz;
	file_in=fopen(nomefile_bin.str().c_str(), "r");

	size_t fread_size = 0;

	fread_size = std::fread(&buff,sizeof(int),1,file_in);
	fread_size = std::fread(&N_param,sizeof(int),1,file_in);
	fread_size = std::fread(&buff,sizeof(int),1,file_in);

	if(out_swap) swap_endian_i(&N_param,1);

	printf("numero parametri=%i\n",N_param);
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

	npe_y=int_param[0];     //numero processori
	npe_z=int_param[1];     //numero processori
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

	printf("=========INIZIO LETTURE==========\n");
	printf("nx1*ny1*nz1: %i %i %i = %i\n",nx1,ny1,nz1,nx1*ny1*nz1);
	fflush(stdout);

	field=(float*)malloc(nx1*ny1*nz1*sizeof(float));
	//	segnox=0;
	segnoy=segnoz=0;
	for(ipz=0;ipz<npe_z;ipz++)
	{
		segnoy=0;
		for(ipy=0;ipy<npe_y;ipy++)
		{
			fread_size = std::fread(&buff,sizeof(int),1,file_in);
			fread_size = std::fread(npoint_loc,sizeof(int),3,file_in);
			fread_size = std::fread(&buff,sizeof(int),1,file_in);

			if(out_swap) swap_endian_i(npoint_loc,3);

			loc_size=npoint_loc[0]*npoint_loc[1]*npoint_loc[2];
			nxloc=npoint_loc[0];
			nyloc=npoint_loc[1];
			nzloc=npoint_loc[2];

			printf("processore ipz=%i/%i  ipy=%i/%i     segnoz=%i     segnoy=%i\r",ipz,npe_z, ipy,npe_y,segnoz,segnoy );
			fflush(stdout);

			buffer=(float*)malloc(loc_size*sizeof(float));
			fread_size = std::fread(&buff,sizeof(int),1,file_in);
			fread_size = std::fread(buffer,sizeof(float),loc_size,file_in);
			fread_size = std::fread(&buff,sizeof(int),1,file_in);

			if(out_swap) swap_endian_f(buffer,loc_size);

			for(k=0;k<nzloc;k++)
				for(j=0;j<nyloc;j++)
					for(i=0;i<nxloc;i++)
						field[i+(j+segnoy)*nx1+(k+segnoz)*nx1*ny1]=buffer[i+j*nxloc+k*nxloc*nyloc];
			segnoy+=nyloc;
		} 
		segnoz+=nzloc;
	}

	dx=(xmax-xmin)/(nx1-1);
	dy=(ymax-ymin)/(ny1-1);
	if(nz1<=1)
		dz=(zmax-zmin);
	else
		dz=(zmax-zmin)/(nz1-1);

	x_coordinates=new float[nx1];
	y_coordinates=new float[ny1];
	z_coordinates=new float[nz1];

	fread_size = std::fread(&buff,sizeof(int),1,file_in);

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
		for(i=0;i<nx1;i++)
			x_coordinates[i]=xmin+dx*i;
		for(i=0;i<ny1;i++)
			y_coordinates[i]=ymin+dy*i;
		for(i=0;i<nz1;i++)
			z_coordinates[i]=zmin+dz*i;
	}

	printf("=========FINE LETTURE==========\n");
	fflush(stdout);



	if(nz1<=1)
	{
		sprintf(nomefile_campi,"%s_out.2D",argv[1]);
		clean_fields=fopen(nomefile_campi, "w");
		printf("\nWriting the fields file 2D (not vtk)\n");

		//output per gnuplot (x:y:valore) compatibile con programmino passe_par_tout togliendo i #
		fprintf(clean_fields,"# %i\n#%i\n#%i\n",nx1, ny1, 1); 
		fprintf(clean_fields,"#%f %f\n#%f %f\n",xmin, ymin, xmax, ymax);
		for(k=0;k<nz1;k++)
			for(j=0;j<ny1;j++)
				for(i=0;i<nx1;i++)
				{
					xx=x_coordinates[i];
					//xmin+dx*i;
					yy=y_coordinates[j];//ymin+dy*j;
					fprintf(clean_fields,"%.4g %.4g %.4g\n",xx, yy, field[i+j*nx1+k*nx1*ny1]);
				}
				fclose(clean_fields);
	}
	else
	{
		sprintf(nomefile_campi,"%s_out.2D",argv[1]);
		clean_fields=fopen(nomefile_campi, "w");
		printf("\nWriting the fields file 2D (not vtk)\n");

		//output per gnuplot (x:y:valore) compatibile con programmino passe_par_tout togliendo i #
		k=nz1/2;
		fprintf(clean_fields,"# 2D cut at z=%g\n",k*dz+zmin); 
		fprintf(clean_fields,"# %i\n#%i\n#%i\n",nx1, ny1, 1); 
		fprintf(clean_fields,"#%f %f\n#%f %f\n",xmin, ymin, xmax, ymax);
		for(j=0;j<ny1;j++)
			for(i=0;i<nx1;i++)
			{
				xx=x_coordinates[i];
				yy=y_coordinates[j];
				//xx=xmin+dx*i;
				//yy=ymin+dy*j;
				fprintf(clean_fields,"%.4g %.4g %.4g\n",xx, yy, field[i+j*nx1+k*nx1*ny1]);
			}
			fclose(clean_fields);
	}
	printf("%lu\nFine\n\n",(unsigned long) fread_size);


	if(parametri->endian_machine == 0)
	{
		swap_endian_f(field,nx1*ny1*nz1);
		swap_endian_f(x_coordinates,nx1);
		swap_endian_f(y_coordinates,ny1);
		swap_endian_f(z_coordinates,nz1);
	}

	//////// DATASET STRUCTURED_POINTS VERSION    ////////
	/*
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
	for(k=0;k<nz1;k++)
	{
	rr[2]=z_coordinates[k];
	for(j=0;j<ny1;j++)
	{
	rr[1]=y_coordinates[j];
	for(i=0;i<nx1;i++)
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
	*/
	if(int new_vtk=0)
	  {
	    printf("new vtk ...\n");
	    fflush(stdout);
	    sprintf(nomefile_campi,"%s_out.vtk",argv[1]);
	    clean_fields=fopen(nomefile_campi, "w");
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

	fclose(file_in);
	}
	else
	{
	  printf("OLD vtk ...\n");
	    fflush(stdout);
	    sprintf(nomefile_campi,"%s_out.vtk",argv[1]);
	  clean_fields=fopen(nomefile_campi, "w");
	  printf("\nWriting the fields file\n");
	  fprintf(clean_fields,"# vtk DataFile Version 2.0\n");
	  fprintf(clean_fields,"titolo mio\n");
	  fprintf(clean_fields,"BINARY\n");
	  fprintf(clean_fields,"DATASET STRUCTURED_POINTS\n");
	  
	  fprintf(clean_fields,"DIMENSIONS %i %i %i\n",nx1, ny1, nz1);
	  fprintf(clean_fields,"ORIGIN %f %f %f\n",xmin, ymin, zmin);
	  fprintf(clean_fields,"SPACING %f %f %f\n",dx, dy, dz);
	  	  
	  fprintf(clean_fields,"POINT_DATA %i\n",nx1*ny1*nz1);
	  fprintf(clean_fields,"SCALARS %s float 1\n",parametri->support_label);
	  fprintf(clean_fields,"LOOKUP_TABLE default\n");
	  fwrite((void*)field,sizeof(float),nx1*ny1*nz1,clean_fields);
	  fclose(clean_fields);
	  
	fclose(file_in);
	}
	return 0;

}

#endif
