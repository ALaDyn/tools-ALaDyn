#ifndef __LEGGI_PARTICELLE_C
#define __LEGGI_PARTICELLE_C

#include "leggi_particelle.h"



int leggi_particelle(char* fileIN, int WEIGHT, int out_swap, int out_binary, int out_ascii, parametri_binnaggio binning)
{
	bool out_parameters = true;

	bool fai_slices = true;
	float whichbinx = 0., whichbinpx = 0.;
	int whichbinx_int = (int) whichbinx;
	int whichbinpx_int = (int) whichbinpx;
	float whichbinE = 0., whichbintheta = 0.;
	int whichbinE_int = (int) whichbinE;
	int whichbintheta_int = (int) whichbintheta;

	int ipc,N_param, *int_param,npart_loc;
	int buff, pID;

	int nelab=3;	//3 valori per ora: gamma, theta ed energia
	int conta_nptot = 0;
	float x,y,z,px,py,pz;

	short buffshort[2];
	float *particelle, *real_param;
	//	float *particelle_elaborate;
	float gamma, theta, E;
	char nomefile_binary[MAX_LENGTH_FILENAME];
	char nomefile_ascii[MAX_LENGTH_FILENAME];
	char nomefile_parametri[MAX_LENGTH_FILENAME];

	FILE *file_in;
	file_in=fopen(fileIN, "r");
	FILE *binary_all_out;
	FILE *ascii_all_out;
	FILE *parameters;

	int npe,nx,ny,nz,ibx,iby,ibz,model,dmodel,nsp,ndim,lpord,deord,nptot, ny_loc, np_loc,ndv;
	float tnow,xmin,xmax,ymin,ymax,zmin,zmax,w0x,w0y,nrat,a0,lam0,E0,ompe,xt_in,xt_end,charge,mass, np_over_nm;
	float rx, ry, rz, ux, uy, uz, wgh;
	int tipo = 0;
	if (fileIN[0] == 'E') tipo = 3;
	else if (fileIN[0] == 'P') tipo = 1;
	else printf("Tipo non riconosciuto!\n");

	fread((void*) &buff,sizeof(int),1,file_in);
	fread((void*) &N_param,sizeof(int),1,file_in);
	fread((void*) &buff,sizeof(int),1,file_in);
	if (out_swap) swap_endian_i(&N_param,1);
	printf("numero parametri=%i\n",N_param);
	fflush(stdout);
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
	nrat=real_param[9];     //n over n critical
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
	printf("nptot=%i\n",int_param[16]); 
	printf("\ninizio processori \n");
	fflush(stdout);
	float **xpx = new float* [binning.nbin_x+3];
	for (int i = 0; i < binning.nbin_x+3; i++)
	{
		xpx[i] = new float [binning.nbin_px+3];
		for (int j = 0; j < binning.nbin_px+3; j++) xpx[i][j] = 0.0;
	}
	float **Etheta = new float* [binning.nbin_E+3];
	for (int i = 0; i < binning.nbin_E+3; i++)
	{
		Etheta[i] = new float [binning.nbin_theta+3];
		for (int j = 0; j < binning.nbin_theta+3; j++) Etheta[i][j] = 0.0;
	}

	for(ipc=0;ipc<npe;ipc++)
	{
		fread(&buff,sizeof(int),1,file_in); 
		fread(&npart_loc,sizeof(int),1,file_in);
		fread(&buff,sizeof(int),1,file_in);
		if (out_swap) swap_endian_i(&npart_loc,1);
		particelle=(float*)malloc(npart_loc*(2*ndim+WEIGHT)*sizeof(float));
		//		particelle_elaborate=(float*)malloc(npart_loc*nelab*sizeof(float));
		printf("proc number %i\tnpart=%i\n",ipc,npart_loc);
		if(npart_loc>0)
		{
			printf("\tentro\t");
			fflush(stdout);
			fread(buffshort,sizeof(short),2,file_in);
			if (out_swap) swap_endian_s(buffshort,2);
			printf("lunghezza=%i    %hu\t%hu\n",npart_loc*2*ndim,buffshort[0],buffshort[1]);
			fread(particelle,sizeof(float),npart_loc*(2*ndim+WEIGHT),file_in);
			fread(&buff,sizeof(int),1,file_in);
			if (out_swap) swap_endian_f(particelle,npart_loc*(2*ndim+WEIGHT));
		}


		if (fai_slices && ndim == 3)
		{
			for (int i = 0; i < npart_loc; i++)
			{
				x=*(particelle+i*(2*ndim+WEIGHT));
				y=*(particelle+i*(2*ndim+WEIGHT)+1);
				z=*(particelle+i*(2*ndim+WEIGHT)+2);
				px=*(particelle+i*(2*ndim+WEIGHT)+3);
				py=*(particelle+i*(2*ndim+WEIGHT)+4);
				pz=*(particelle+i*(2*ndim+WEIGHT)+5);
				//				particelle_elaborate[i]=(float)(sqrt(1.+px*px+py*py+pz*pz)-1.);				//gamma
				//				particelle_elaborate[i+1]=(float)(atan2(sqrt(py*py+pz*pz),px)*180./M_PI);	//theta nb: py e pz sono quelli trasversi in ALaDyn!
				//				particelle_elaborate[i+2]=(float)(particelle_elaborate[i]*P_MASS);			//energia
				gamma=(float)(sqrt(1.+px*px+py*py+pz*pz)-1.);				//gamma
				theta=(float)(atan2(sqrt(py*py+pz*pz),px)*180./M_PI);	//theta nb: py e pz sono quelli trasversi in ALaDyn!
				E=(float)(gamma*P_MASS);			//energia

				// xpx
				if (x < binning.xmin)
				{
					//					whichbinx = 0.0;
					whichbinx_int = 0;
				}
				else if (x > binning.xmax)
				{
					//				  whichbinx = (float) (binning.nbin_x + 2);
					whichbinx_int = binning.nbin_x + 2;
				}
				else
				{
					whichbinx = (x - binning.xmin) / binning.dimmi_dimx();
					whichbinx_int = (int)(whichbinx+1.0);
				}
				if (px < binning.pxmin)
				{
					//					whichbinpx = 0.0;
					whichbinpx_int = 0;
				}
				else if (px > binning.pxmax)
				{
					//				  whichbinpx = (float) (binning.nbin_px + 2);
					whichbinpx_int = binning.nbin_px + 2;
				}
				else
				{
					whichbinpx = (px - binning.pxmin) / binning.dimmi_dimpx();
					whichbinpx_int = (int)(whichbinpx+1.0);
				}
				if (WEIGHT) xpx[whichbinx_int][whichbinpx_int] += *(particelle+i*(2*ndim+WEIGHT)+6);
				else		xpx[whichbinx_int][whichbinpx_int] += 1.0;

				// Etheta
				if (E < binning.Emin)
				{
					//					whichbinE = 0.0;
					whichbinE_int = 0;
				}
				else if (E > binning.Emax)
				{
					//				  whichbinE = (float) (binning.nbin_E + 2);
					whichbinE_int = binning.nbin_E + 2;
				}
				else
				{
					whichbinE = (E - binning.Emin) / binning.dimmi_dimE();
					whichbinE_int = (int)(whichbinE+1.0);
				}
				if (theta < binning.thetamin)
				{
					//					whichbintheta = 0.0;
					whichbintheta_int = 0;
				}
				else if (theta > binning.thetamax)
				{
					//				  whichbintheta = (float) (binning.nbin_theta + 2);
					whichbintheta_int = binning.nbin_theta + 2;
				}
				else
				{
					whichbintheta = (theta - binning.thetamin) / binning.dimmi_dimtheta();
					whichbintheta_int = (int)(whichbintheta+1.0);
				}
				if (WEIGHT) Etheta[whichbinE_int][whichbintheta_int] += *(particelle+i*(2*ndim+WEIGHT)+6);
				else		Etheta[whichbinE_int][whichbintheta_int] += 1.0;

			}
		}

		if (out_binary)
		{
			sprintf(nomefile_binary,"%s_7_out",fileIN);
			binary_all_out=fopen(nomefile_binary, "ab");
			printf("\nWriting the C binary file\n");
			fwrite((void*)particelle,sizeof(float),nptot*(2*ndim+WEIGHT),binary_all_out);
			fflush(binary_all_out);
			fclose(binary_all_out);
		}




		if(out_ascii)
		{
			sprintf(nomefile_ascii,"%s.ascii",fileIN);
			ascii_all_out=fopen(nomefile_ascii, "a");

			for(int i=0;i<nptot;i++)
			{
				rx=particelle[i*(6+WEIGHT)+1]*((float)1.e-4);
				ry=particelle[i*(6+WEIGHT)+2]*((float)1.e-4);
				rz=particelle[i*(6+WEIGHT)+0]*((float)1.e-4);
				ux=particelle[i*(6+WEIGHT)+4];
				uy=particelle[i*(6+WEIGHT)+5];
				uz=particelle[i*(6+WEIGHT)+3];
				if (WEIGHT)
				{
					wgh=particelle[i*(6+WEIGHT)+6];
					fprintf(ascii_all_out,"%e %e %e %e %e %e %d %e 0 %d\n",rx, ry, rz, ux, uy, uz, tipo, wgh, i+1);
				}
				else
				{
					fprintf(ascii_all_out,"%e %e %e %e %e %e %d 1 0 %d\n",rx, ry, rz, ux, uy, uz, tipo, i+1);
				}
			}
			fflush(ascii_all_out);
			fclose(ascii_all_out);
		}

		free(particelle);
		//		free(particelle_elaborate);
	}


	if (out_parameters)
	{
		sprintf(nomefile_parametri,"%s.parameters",fileIN);
		parameters=fopen(nomefile_parametri, "w");
		printf("\nWriting the parameters file\n");
		fprintf(parameters,"interi\n");
		fprintf(parameters,"npe=%i\n",int_param[0]);     //numero processori
		fprintf(parameters,"nx=%i\n",int_param[1]);
		fprintf(parameters,"ny=%i\n",int_param[2]);
		fprintf(parameters,"nz=%i\n",int_param[3]);
		fprintf(parameters,"ibx=%i\n",int_param[4]);
		fprintf(parameters,"iby=%i\n",int_param[5]);
		fprintf(parameters,"ibz=%i\n",int_param[6]);
		fprintf(parameters,"model=%i\n",int_param[7]);  //modello di laser utilizzato
		fprintf(parameters,"dmodel=%i\n",int_param[8]); //modello di condizioni iniziali
		fprintf(parameters,"nsp=%i\n",int_param[9]);    //numero di speci
		fprintf(parameters,"ndim=%i\n",int_param[10]);   
		fprintf(parameters,"np_loc=%i\n",int_param[11]);  
		fprintf(parameters,"lpord=%i\n",int_param[12]); //ordine dello schema leapfrog
		fprintf(parameters,"deord=%i\n",int_param[13]); //ordine derivate
		fprintf(parameters,"pID=%i\n",int_param[15]); 
		fprintf(parameters,"nptot=%i\n",int_param[16]); 
		fprintf(parameters,"ndv=%i\n",int_param[17]); 
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
		fprintf(parameters,"nrat=%f\n",real_param[9]);     //n over n critical
		fprintf(parameters,"a0=%f\n",real_param[10]);      // a0 laser
		fprintf(parameters,"lam0=%f\n",real_param[11]);    // lambda
		fprintf(parameters,"E0=%f\n",real_param[12]);      //conversione da campi numerici a TV/m
		fprintf(parameters,"ompe=%f\n",real_param[13]);    //costante accoppiamento correnti campi
		fprintf(parameters,"np_over_nm=%f\n",real_param[14]);   //numerical2physical particles 14 
		fprintf(parameters,"xt_in=%f\n",real_param[15]);
		fprintf(parameters,"xt_end=%f\n",real_param[16]);
		fprintf(parameters,"charge=%f\n",real_param[17]);  //carica particella su carica elettrone
		fprintf(parameters,"mass=%f\n",real_param[18]);    //massa particelle su massa elettrone
		//	if(WEIGHT) fprintf(parameters,"weight=%f\n",particelle[6]);    //massa particelle su massa elettrone
		fclose(parameters);
	}


	return 0;
}


#endif