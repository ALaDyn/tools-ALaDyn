
#include "leggi_binario_ALaDyn_fortran.h"
#include "leggi_particelle.h"



int leggi_particelle(char* fileIN, int WEIGHT, int FLAG_ENDIAN, int out_swap, int out_file)
{
//	int N_particles, j, k;
	int i,ipc,N_param, *int_param,npart_loc;
	int segnalibro=0, buff, pID;
	short buffshort[2];
	float *particelle, *real_param;
//	char nome[100],nome4[100];
	char nome1[200],nome2[100],nome3[100],nome5[100],nome6[100];
//	FILE *f, *f3, *f4;
	FILE *file_in,*f1,*f2,*f5,*f6;
//  int nppc,rkord;
	int npe,nx,ny,nz,ibx,iby,ibz,model,dmodel,nsp,ndim,lpord,deord,nptot, ny_loc, np_loc,ndv;
//	float num2phys;
	float tnow,xmin,xmax,ymin,ymax,zmin,zmax,w0x,w0y,nrat,a0,lam0,E0,ompe,xt_in,xt_end,charge,mass, np_over_nm;
	float rx, ry, rz,ux,uy,uz,gamma;
//	double theta;

	file_in=fopen(fileIN, "r");

/*
	sscanf(args[2],"%i",&FLAG_ENDIAN);
	printf("FLAG_ENDIAN=%i\n",FLAG_ENDIAN);
	for (i=3; i<narg; i++)
		if (!strncmp(args[i], "-swap",5))
			out_swap=1;
	for (i=3; i<narg; i++)
		if (!strncmp(args[i], "-outfiles",9))
			out_file=1;
*/


	if(FLAG_ENDIAN==1)
	{
		fread(&buff,sizeof(int),1,file_in);
		fread(&N_param,sizeof(int),1,file_in);
		fread(&buff,sizeof(int),1,file_in);
		printf("numero parametri=%i\n",N_param);
		int_param=(int*)malloc(N_param*sizeof(int));
		real_param=(float*)malloc(N_param*sizeof(float));
		fread(&buff,sizeof(int),1,file_in);
		fread(int_param,sizeof(int),N_param,file_in);
		fread(&buff,sizeof(int),1,file_in);
		fread(&buff,sizeof(int),1,file_in);
		fread(real_param,sizeof(float),N_param,file_in);
		fread(&buff,sizeof(int),1,file_in);
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
			printf("proc number %i\tnpart=%i     segnalibro=%i\n",ipc,npart_loc,segnalibro);
			if(npart_loc>0)
			{
				printf("\tentro\t");
				fflush(stdout);
//				fread(&buff,sizeof(int),1,file_in);
				fread(buffshort,sizeof(short),2,file_in);
//				swap_endian_s(buffshort,2);
				printf("lunghezza=%i    %hu\t%hu\n",npart_loc*2*ndim,buffshort[0],buffshort[1]);
				fread(particelle+segnalibro,sizeof(float),npart_loc*(2*ndim+WEIGHT),file_in);
				fread(&buff,sizeof(int),1,file_in);
				segnalibro+=npart_loc*(2*ndim+WEIGHT);
			}
		}
		printf("=========FINE LETTURE==========\n");
		fflush(stdout);
		sprintf(nome1,"%s_7_out",fileIN);
		if(out_swap)
		{
			swap_endian_f(particelle,nptot*(2*ndim+WEIGHT));
			sprintf(nome1,"%s_7_out_swap",fileIN);
		}
		f1=fopen(nome1, "w");
		fwrite((void*)particelle,sizeof(float),nptot*(2*ndim+WEIGHT),f1);
		fflush(f1);
    }


	else if(FLAG_ENDIAN==2)
	{
		fread(&buff,sizeof(int),1,file_in);
		fread(&N_param,sizeof(int),1,file_in);
		fread(&buff,sizeof(int),1,file_in);
		swap_endian_i(&N_param,1);
		printf("numero parametri=%i\n",N_param);
		int_param=(int*)malloc(N_param*sizeof(int));
		real_param=(float*)malloc(N_param*sizeof(float));
		fread(&buff,sizeof(int),1,file_in);
		fread(int_param,sizeof(int),N_param,file_in);
		fread(&buff,sizeof(int),1,file_in);
		fread(&buff,sizeof(int),1,file_in);
		fread(real_param,sizeof(float),N_param,file_in);
		fread(&buff,sizeof(int),1,file_in);
		swap_endian_i(int_param,N_param);
		swap_endian_f(real_param,N_param);
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
			swap_endian_i(&npart_loc,1);
			printf("proc number %i\tnpart=%i     segnalibro=%i\n",ipc,npart_loc,segnalibro);
			if(npart_loc>0)
			{
				printf("\tentro\t");
				fflush(stdout);
//				fread(&buff,sizeof(int),1,file_in);
				fread(buffshort,sizeof(short),2,file_in);
				swap_endian_s(buffshort,2);
				printf("lunghezza=%i    %hu\t%hu\n",npart_loc*2*ndim,buffshort[0],buffshort[1]);
				fread(particelle+segnalibro,sizeof(float),npart_loc*(2*ndim+WEIGHT),file_in);
				fread(&buff,sizeof(int),1,file_in);
				swap_endian_f(particelle+segnalibro,npart_loc*(2*ndim+WEIGHT));
				segnalibro+=npart_loc*(2*ndim+WEIGHT);
			}
		}    

		printf("=========FINE LETTURE==========\n");
		fflush(stdout);
		sprintf(nome1,"%s_7_out",fileIN);
		f1=fopen(nome1, "w");
		fwrite((void*)particelle,sizeof(float),nptot*(2*ndim+WEIGHT),f1);
		fflush(f1);
	}
	printf("n_proc=%i,   n_dim=%i,  n_ptot=%i\n",npe,ndim,nptot);
	fflush(stdout);
	if(out_file)
	{
//		sprintf(nome,"%s_XY",fileIN);
//		f=fopen(nome, "w");
//		sprintf(nome1,"%s_XPX",fileIN);
//		f1=fopen(nome1, "w");
		sprintf(nome2,"%s_PXPYPZ",fileIN);
		f2=fopen(nome2, "w");
		sprintf(nome3,"%s_thetaENERGY",fileIN);
//		f3=fopen(nome3, "w");
//		sprintf(nome4,"%s_XYENERGY",fileIN);
//		f4=fopen(nome4, "w");
		sprintf(nome5,"%s_ALL_7_bin",fileIN);
		f5=fopen(nome5, "w");
		sprintf(nome6,"%s_ALL_6",fileIN);
		f6=fopen(nome6, "w");

		for(i=0;i<nptot;i++)
		{
			rx=particelle[i*(6+WEIGHT)+0];
			ry=particelle[i*(6+WEIGHT)+1];
			rz=particelle[i*(6+WEIGHT)+2];
			ux=particelle[i*(6+WEIGHT)+3];
			uy=particelle[i*(6+WEIGHT)+4];
			uz=particelle[i*(6+WEIGHT)+5];
			gamma=(float) sqrt(1+ux*ux+uy*uy+uz*uz)-1;
//			fprintf(f,"%e %e\n",rx[i],ry[i]);
//			fprintf(f1,"%e %e\n",rx[i],ux[i]);
			fprintf(f2,"%e %e %e\n",ux,uy,uz);
//			fprintf(f3,"%e %e\n",theta,(gamma[i]-1)*938);//*mass*0.511);
//			fprintf(f4,"%e %e %e\n",ry[i],rx[i],(gamma[i]-1)*938);//*mass*0.511);
//			tagli in energia commentati perche' li faccio poi con altri tools
//			if(gamma*938>2)
//			{
				fprintf(f6,"%e %e %e %e %e %e\n",rx, ry, rz, ux, uy, uz);//mass*0.511);
				fwrite(particelle+i*(6+WEIGHT),sizeof(float),7,f5);
//			}
		}

		fclose(f2);
		fclose(f6);
	}
	return 0;
}
