
#include "leggi_particelle.h"


int leggi_particelle(char* fileIN, int WEIGHT, int FLAG_ENDIAN, int out_swap, int out_binary, int out_ascii)
{
	int i,ipc,N_param, *int_param,npart_loc;
	int segnalibro=0, buff, pID;
	short buffshort[2];
	float *particelle, *real_param;
	char nomefile_binary[MAX_LENGTH_FILENAME];
	char nomefile_ascii[MAX_LENGTH_FILENAME];
//	char nomefile_parametri[MAX_LENGTH_FILENAME];

	FILE *file_in;
	file_in=fopen(fileIN, "rb");

	FILE *binary_all_out;
	FILE *ascii_all_out;
//	FILE *parameters;

	int npe,nx,ny,nz,ibx,iby,ibz,model,dmodel,nsp,ndim,lpord,deord,nptot, ny_loc, np_loc,ndv;
	float tnow,xmin,xmax,ymin,ymax,zmin,zmax,w0x,w0y,nrat,a0,lam0,E0,ompe,xt_in,xt_end,charge,mass, np_over_nm;
	float rx, ry, rz,ux,uy,uz;

	fread(&buff,sizeof(int),1,file_in);
	fread(&N_param,sizeof(int),1,file_in);
	fread(&buff,sizeof(int),1,file_in);
	if (out_swap) swap_endian_i(&N_param,1);
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
	nrat=real_param[9];     //n orver n critical
	a0=real_param[10];      // a0 laser
	lam0=real_param[11];    // lambda
	E0=real_param[12];      //conversione da campi numerici a TV/m
	ompe=real_param[13];    //costante accoppiamento correnti campi
	np_over_nm=real_param[14];   //numerical2physical particles 14 
	xt_in=real_param[15];   //inizio plasma
	xt_end=real_param[16];
	charge=real_param[17];  //carica particella su carica elettrone
	mass=real_param[18];    //massa particelle su massa elettrone
	ny=ny_loc*npe;
	printf("\ninizio processori \n");
	fflush(stdout);
	particelle=(float*)malloc(nptot*(2*ndim+WEIGHT)*sizeof(float));
	segnalibro=0;

	for(ipc=0;ipc<npe;ipc++)
	{
		fread(&buff,sizeof(int),1,file_in); 
		fread(&npart_loc,sizeof(int),1,file_in);
		fread(&buff,sizeof(int),1,file_in);
		if (out_swap) swap_endian_i(&npart_loc,1);
		printf("proc number %i\tnpart=%i     segnalibro=%i\n",ipc,npart_loc,segnalibro);
		if(npart_loc>0)
		{
			printf("\tentro\t");
			fflush(stdout);
			fread(buffshort,sizeof(short),2,file_in);
			if (out_swap) swap_endian_s(buffshort,2);
			printf("lunghezza=%i    %hu\t%hu\n",npart_loc*2*ndim,buffshort[0],buffshort[1]);
			fread(particelle+segnalibro,sizeof(float),npart_loc*(2*ndim+WEIGHT),file_in);
			fread(&buff,sizeof(int),1,file_in);
			if (out_swap) swap_endian_f(particelle+segnalibro,npart_loc*(2*ndim+WEIGHT));
			segnalibro+=npart_loc*(2*ndim+WEIGHT);
		}
	}    

	printf("=========FINE LETTURE==========\n");
	printf("n_proc=%i,   n_dim=%i,  n_ptot=%i\n",npe,ndim,nptot);
	fflush(stdout);

	/*
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
	*/

	if (out_binary)
	{
		sprintf(nomefile_binary,"%s_7_out",fileIN);
		binary_all_out=fopen(nomefile_binary, "wb");
		printf("\nWriting the clean binary file\n");
		fwrite((void*)particelle,sizeof(float),nptot*(2*ndim+WEIGHT),binary_all_out);
		fflush(binary_all_out);
		fclose(binary_all_out);
	}




	if(out_ascii)
	{
		sprintf(nomefile_ascii,"%s.ascii",fileIN);
		ascii_all_out=fopen(nomefile_ascii, "w");
		printf("\nWriting the ascii file\n");

		for(i=0;i<nptot;i++)
		{
			rx=particelle[i*(6+WEIGHT)+0];
			ry=particelle[i*(6+WEIGHT)+1];
			rz=particelle[i*(6+WEIGHT)+2];
			ux=particelle[i*(6+WEIGHT)+3];
			uy=particelle[i*(6+WEIGHT)+4];
			uz=particelle[i*(6+WEIGHT)+5];

			fprintf(ascii_all_out,"%e %e %e %e %e %e\n",rx, ry, rz, ux, uy, uz);
		}
		fflush(ascii_all_out);
		fclose(ascii_all_out);
	}

	return 0;
}
