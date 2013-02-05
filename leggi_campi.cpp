#ifndef __LEGGI_CAMPI_C
#define __LEGGI_CAMPI_C

#include "leggi_campi.h"

int leggi_campi(char* fileIN, parametri binning)
{
	int out_swap = binning.p[SWAP];

	int i, ipy, ipz,  j, k,np_loc, N_param, *int_param,npoint_loc[3], loc_size, kk;
	int segnox=0,segnoy=0,segnoz=0, buff;
	float *field,*buffer,  *real_param;
	char nomefile_parametri[MAX_LENGTH_FILENAME];
	char nomefile_campi[MAX_LENGTH_FILENAME];

	FILE *file_in;
	FILE *parameters;
	FILE *clean_fields;

	int npe,nx,ibx,iby,ibz,model,dmodel,nsp,lpord,deord,npe_y, npe_z;
	int nxloc, nx1, ny1, nyloc, nz1, nzloc, fvar;
	float tnow,xmin,xmax,ymin,ymax,zmin,zmax,w0x,w0y,nrat,a0,lam0,E0,B0,ompe,xt_in,xt_end,charge,mass;
	float dx, dy, dz;
	file_in=fopen(fileIN, "r");


	fread(&buff,sizeof(int),1,file_in);
	fread(&N_param,sizeof(int),1,file_in);
	fread(&buff,sizeof(int),1,file_in);

	if(out_swap) swap_endian_i(&N_param,1);

	printf("numero parametri=%i\n",N_param);
	int_param=(int*)malloc(N_param*sizeof(int));
	real_param=(float*)malloc(N_param*sizeof(float));
	fread(&buff,sizeof(int),1,file_in);
	fread(int_param,sizeof(int),N_param,file_in);
	fread(&buff,sizeof(int),1,file_in);
	fread(&buff,sizeof(int),1,file_in);
	fread(real_param,sizeof(float),N_param,file_in);
	fread(&buff,sizeof(int),1,file_in);

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

	sprintf(nomefile_parametri,"%s.parameters",fileIN);
	parameters=fopen(nomefile_parametri, "w");
	printf("\nWriting the parameters file\n");
	fprintf(parameters,"interi\n");
	fprintf(parameters,"npe_y=%i\n",int_param[0]);     //numero processori
	fprintf(parameters,"npe_z=%i\n",int_param[1]);     //numero processori
	fprintf(parameters,"npe=%i\n",npe_y*npe_z);     //numero processori
	fprintf(parameters,"nx=%i\n",int_param[2]);
	fprintf(parameters,"nx1=%i\n",int_param[3]);
	fprintf(parameters,"ny1=%i\n",int_param[4]);
	fprintf(parameters,"nyloc=%i\n",int_param[5]);
	fprintf(parameters,"nz1=%i\n",int_param[6]);
	fprintf(parameters,"nzloc=%i\n",int_param[7]);
	fprintf(parameters,"ibx=%i\n",int_param[8]);
	fprintf(parameters,"iby=%i\n",int_param[9]);
	fprintf(parameters,"ibz=%i\n",int_param[10]);
	fprintf(parameters,"model=%i\n",int_param[11]);  //modello di laser utilizzato
	fprintf(parameters,"dmodel=%i\n",int_param[12]); //modello di condizioni iniziali
	fprintf(parameters,"nsp=%i\n",int_param[13]);    //numero di speci
	fprintf(parameters,"np_loc=%i\n",int_param[14]);  //numero di componenti dello spazio dei momenti
	fprintf(parameters,"lpord=%i\n",int_param[15]); //ordine dello schema leapfrog
	fprintf(parameters,"deord=%i\n",int_param[16]); //ordine derivate
	fprintf(parameters,"fvar=%i\n",int_param[17]); 
	fprintf(parameters,"========= fine interi\n");
	fprintf(parameters,"\n floating\n");
	fprintf(parameters,"tnow=%f\n",real_param[0]);  //tempo dell'output
	fprintf(parameters,"xmin=%f\n",real_param[1]);  //estremi della griglia
	fprintf(parameters,"xmax=%f\n",real_param[2]);  //estremi della griglia
	fprintf(parameters,"ymin=%f\n",real_param[3]);  //estremi della griglia
	fprintf(parameters,"ymax=%f\n",real_param[4]);  //estremi della griglia
	fprintf(parameters,"zmin=%f\n",real_param[5]);  //estremi della griglia
	fprintf(parameters,"zmax=%f\n",real_param[6]);  //estremi della griglia
	fprintf(parameters,"w0x=%f\n",real_param[7]);      //waist del laser in x
	fprintf(parameters,"w0y=%f\n",real_param[8]);      //waist del laser in y
	fprintf(parameters,"nrat=%f\n",real_param[9]);     //n orver n critical
	fprintf(parameters,"a0=%f\n",real_param[10]);      // a0 laser
	fprintf(parameters,"lam0=%f\n",real_param[11]);    // lambda
	fprintf(parameters,"E0=%f\n",real_param[12]);      //conversione da campi numerici a TV/m
	fprintf(parameters,"ompe=%f\n",real_param[13]);    //costante accoppiamento correnti campi
	fprintf(parameters,"xt_in=%f\n",real_param[14]);   //inizio plasma
	fprintf(parameters,"xt_end=%f\n",real_param[15]);
	fprintf(parameters,"charge=%f\n",real_param[16]);  //carica particella su carica elettrone
	fprintf(parameters,"mass=%f\n",real_param[17]);    //massa particelle su massa elettrone
	fclose(parameters);

	printf("=========INIZIO LETTURE==========\n");
	printf("nx1*ny1*nz1: %i %i %i = %i\n",nx1,ny1,nz1,nx1*ny1*nz1);
	fflush(stdout);

	field=(float*)malloc(nx1*ny1*nz1*sizeof(float));
	segnox=segnoy=segnoz=0;
	for(ipz=0;ipz<npe_z;ipz++)
	{
		segnoy=0;
		for(ipy=0;ipy<npe_y;ipy++)
		{
			fread(&buff,sizeof(int),1,file_in);
			fread(npoint_loc,sizeof(int),3,file_in);
			fread(&buff,sizeof(int),1,file_in);

			if(out_swap) swap_endian_i(npoint_loc,3);

			loc_size=npoint_loc[0]*npoint_loc[1]*npoint_loc[2];
			nxloc=npoint_loc[0];
			nyloc=npoint_loc[1];
			nzloc=npoint_loc[2];

			printf("processore ipz=%i/%i  ipy=%i/%i     segnoz=%i     segnoy=%i\r",ipz,npe_z, ipy,npe_y,segnoz,segnoy );
			fflush(stdout);

			buffer=(float*)malloc(loc_size*sizeof(float));
			fread(&buff,sizeof(int),1,file_in);
			fread(buffer,sizeof(float),loc_size,file_in);
			fread(&buff,sizeof(int),1,file_in);

			if(out_swap) swap_endian_f(buffer,loc_size);

			kk=0;
			for(k=0;k<nzloc;k++)
				for(j=0;j<nyloc;j++)
					for(i=0;i<nxloc;i++)
						field[i+(j+segnoy)*nx1+(k+segnoz)*nx1*ny1]=buffer[i+j*nxloc+k*nxloc*nyloc];
			segnoy+=nyloc;
		} 
		segnoz+=nzloc;
	}

	if(out_swap)
	{
		swap_endian_f(field,nx1*ny1*nz1);
	}

	printf("=========FINE LETTURE==========\n");
	fflush(stdout);

	sprintf(nomefile_campi,"%s_out.vtk",fileIN);
	clean_fields=fopen(nomefile_campi, "w");
	printf("\nWriting the fields file\n");
	fprintf(clean_fields,"# vtk DataFile Version 2.0\n");
	fprintf(clean_fields,"titolo mio\n");
	fprintf(clean_fields,"BINARY\n");
	fprintf(clean_fields,"DATASET STRUCTURED_POINTS\n");
	fprintf(clean_fields,"DIMENSIONS %i %i %i\n",nx1, ny1, nz1);
	fprintf(clean_fields,"ORIGIN %f %f %f\n",xmin, ymin, zmin);
	dx=(xmax-xmin)/(nx1-1);
	dy=(ymax-ymin)/(ny1-1);
	dz=(zmax-zmin)/(nz1-1);
	fprintf(clean_fields,"SPACING %f %f %f\n",dx, dy, dz);
	fprintf(clean_fields,"POINT_DATA %i\n",nx1*ny1*nz1);
	fprintf(clean_fields,"SCALARS %s float 1\n",binning.support_label);
	fprintf(clean_fields,"LOOKUP_TABLE default\n");
	fwrite((void*)field,sizeof(float),nx1*ny1*nz1,clean_fields);
	fclose(clean_fields);

	return 0;
}

#endif