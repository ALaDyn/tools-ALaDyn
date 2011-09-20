
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))
#define TRUE 1
#define FALSE 0

void swap_endian_s(short* in_s,int n)
{
  int i;
  union {short smio; char arr[2];}x;
  char buff;

  for(i=0;i<n;i++)
    {
      x.smio=in_s[i];
      
      buff=x.arr[0];
      x.arr[0]=x.arr[1];
      x.arr[1]=buff;
      
      in_s[i]=x.smio;
    }
}
void swap_endian_i(int* in_i,int n)
{
  int i;
  union {int imio; float fmio; char arr[4];}x;
  char buff;

  for(i=0;i<n;i++)
    {
      x.imio=in_i[i];
      
      buff=x.arr[0];
      x.arr[0]=x.arr[3];
      x.arr[3]=buff;
      buff=x.arr[1];
      x.arr[1]=x.arr[2];
      x.arr[2]=buff;
      
      in_i[i]=x.imio;
    }
}
void swap_endian_f(float* in_f, int n)
{
  int i;
  union {int imio; float fmio; char arr[4];}x;
  char buff;
  for(i=0;i<n;i++)
    {
      x.fmio=in_f[i];
      
      buff=x.arr[0];
      x.arr[0]=x.arr[3];
      x.arr[3]=buff;
      buff=x.arr[1];
      x.arr[1]=x.arr[2];
      x.arr[2]=buff;
      
      in_f[i]=x.fmio;
    }
}

//typedef  union {int imio, float fmio, char arr[4]} piripu;


int main(int narg, char **args)
{
  int i,ipx, ipy, ipz,  j, k,np_loc, N_particles, N_param, *int_param,npoint_loc[3], loc_size, kk;
  int out_swap=0, segnox=0,segnoy=0,segnoz=0,FLAG_ENDIAN=0, buff, pID;
  short buffshort[2];
  float *field,*buffer,  *real_param;
  char nome[100];
  FILE *f,*file_in,*f1,*f2,*f3,*f4,*f5,*f6;
  int npe,nx,ny,nz,ibx,iby,ibz,model,dmodel,nsp,ndim,nppc,rkord,lpord,deord,nptot, npe_y, npe_z;
  int nxloc, nx1, ny1, nyloc, nz1, nzloc, fvar;
  float tnow,xmin,xmax,ymin,ymax,zmin,zmax,w0x,w0y,nrat,a0,lam0,E0,B0,ompe,xt_in,xt_end,charge,mass, num2phys;
  float *rx, *ry, *rz,*ux,*uy,*uz,*gamma, dx, dy, dz;
  double theta;
  file_in=fopen(args[1], "r");
  
  for (i=2; i<narg; i++)
    if (!strncmp(args[i], "-endian",7))
      FLAG_ENDIAN=1;
  printf("FLAG_ENDIAN=%i\n",FLAG_ENDIAN);
  
  for (i=2; i<narg; i++)
    if (!strncmp(args[i], "-swap",5))
      out_swap=1;
  
  fread(&buff,sizeof(int),1,file_in);
  fread(&N_param,sizeof(int),1,file_in);
  fread(&buff,sizeof(int),1,file_in);
  
  if(FLAG_ENDIAN)
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
  
  if(FLAG_ENDIAN)
    swap_endian_i(int_param,N_param),
      swap_endian_f(real_param,N_param);
            
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
  B0=E0+33.3;
  ompe=real_param[13];    //costante accoppiamento correnti campi
  xt_in=real_param[14];   //inizio plasma
  xt_end=real_param[15];
  charge=real_param[16];  //carica particella su carica elettrone
  mass=real_param[17];    //massa particelle su massa elettrone
  

  sprintf(nome,"%s_parametri",args[1]);
  f1=fopen(nome, "w");
  fprintf(f1,"interi\n");
  fprintf(f1,"npe_y=%i\n",int_param[0]);     //numero processori
  fprintf(f1,"npe_z=%i\n",int_param[1]);     //numero processori
  fprintf(f1,"npe=%i\n",npe_y*npe_z);     //numero processori
  fprintf(f1,"nx=%i\n",int_param[2]);
  fprintf(f1,"nx1=%i\n",int_param[3]);
  fprintf(f1,"ny1=%i\n",int_param[4]);
  fprintf(f1,"nyloc=%i\n",int_param[5]);
  fprintf(f1,"nz1=%i\n",int_param[6]);
  fprintf(f1,"nzloc=%i\n",int_param[7]);
  fprintf(f1,"ibx=%i\n",int_param[8]);
  fprintf(f1,"iby=%i\n",int_param[9]);
  fprintf(f1,"ibz=%i\n",int_param[10]);
  fprintf(f1,"model=%i\n",int_param[11]);  //modello di laser utilizzato
  fprintf(f1,"dmodel=%i\n",int_param[12]); //modello di condizioni iniziali
  fprintf(f1,"nsp=%i\n",int_param[13]);    //numero di speci
  fprintf(f1,"np_loc=%i\n",int_param[14]);  //numero di componenti dello spazio dei momenti
  fprintf(f1,"lpord=%i\n",int_param[15]); //ordine dello schema leapfrog
  fprintf(f1,"deord=%i\n",int_param[16]); //ordine derivate
  fprintf(f1,"fvar=%i\n",int_param[17]); 
  fprintf(f1,"========= fine interi\n");
  fprintf(f1,"\n floating\n");
  fprintf(f1,"tnow=%f\n",real_param[0]);  //tempo dell'output
  fprintf(f1,"xmin=%f\n",real_param[1]);  //estremi della griglia
  fprintf(f1,"xmax=%f\n",real_param[2]);  //estremi della griglia
  fprintf(f1,"ymin=%f\n",real_param[3]);  //estremi della griglia
  fprintf(f1,"ymax=%f\n",real_param[4]);  //estremi della griglia
  fprintf(f1,"zmin=%f\n",real_param[5]);  //estremi della griglia
  fprintf(f1,"zmax=%f\n",real_param[6]);  //estremi della griglia
  fprintf(f1,"w0x=%f\n",real_param[7]);      //waist del laser in x
  fprintf(f1,"w0y=%f\n",real_param[8]);      //waist del laser in y
  fprintf(f1,"nrat=%f\n",real_param[9]);     //n orver n critical
  fprintf(f1,"a0=%f\n",real_param[10]);      // a0 laser
  fprintf(f1,"lam0=%f\n",real_param[11]);    // lambda
  fprintf(f1,"E0=%f\n",real_param[12]);      //conversione da campi numerici a TV/m
  fprintf(f1,"ompe=%f\n",real_param[13]);    //costante accoppiamento correnti campi
  fprintf(f1,"xt_in=%f\n",real_param[14]);   //inizio plasma
  fprintf(f1,"xt_end=%f\n",real_param[15]);
  fprintf(f1,"charge=%f\n",real_param[16]);  //carica particella su carica elettrone
  fprintf(f1,"mass=%f\n",real_param[17]);    //massa particelle su massa elettrone
  

  fclose(f1);

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
	  
	  if(FLAG_ENDIAN)
	    swap_endian_i(npoint_loc,3);
	  loc_size=npoint_loc[0]*npoint_loc[1]*npoint_loc[2];
	  nxloc=npoint_loc[0];
	  nyloc=npoint_loc[1];
	  nzloc=npoint_loc[2];
	  
	  printf("processore ipz=%i/%i  ipy=%i/%i     segnoz=%i     segnoy=%i\r",ipz,npe_z, ipy,npe_y,segnoz,segnoy );
	  //printf("\t\t nxloc=%i  nyloc=%i  nzloc=%i\n",nxloc,nyloc,nzloc);
	  fflush(stdout);
 	  free((void*)buffer);
	  buffer=(float*)malloc(loc_size*sizeof(float));
	  
	  fread(&buff,sizeof(int),1,file_in);
	  fread(buffer,sizeof(float),loc_size,file_in);
	  fread(&buff,sizeof(int),1,file_in);
	  if(FLAG_ENDIAN)
	    swap_endian_f(buffer,loc_size);
	  kk=0;
	  for(k=0;k<nzloc;k++)
	    for(j=0;j<nyloc;j++)
	      for(i=0;i<nxloc;i++)
		field[i+(j+segnoy)*nx1+(k+segnoz)*nx1*ny1]=buffer[i+j*nxloc+k*nxloc*nyloc];
	  
	  segnoy+=nyloc;
	} 
      segnoz+=nzloc;
    }
        
  printf("=========FINE LETTURE==========\n");
  fflush(stdout);
 
  sprintf(nome,"%s_out",args[1]);
  f1=fopen(nome, "w");
  if(out_swap)
    {
      swap_endian_f(field,nx1*ny1*nz1);
      sprintf(nome,"%s_out_swap",args[1]);
    }
  fprintf(f1,"# vtk DataFile Version 2.0\n");
  fprintf(f1,"titolo mio\n");
  fprintf(f1,"BINARY\n");
  fprintf(f1,"DATASET STRUCTURED_POINTS\n");
  fprintf(f1,"DIMENSIONS %i %i %i\n",nx1, ny1, nz1);
  fprintf(f1,"ORIGIN %f %f %f\n",xmin, ymin, zmin);
  dx=(xmax-xmin)/(nx1-1);
  dy=(ymax-ymin)/(ny1-1);
  dz=(zmax-zmin)/(nz1-1);
  fprintf(f1,"SPACING %f %f %f\n",dx, dy, dz);
  fprintf(f1,"POINT_DATA %i\n",nx1*ny1*nz1);
  fprintf(f1,"SCALARS Ex float 1\n");
  fprintf(f1,"LOOKUP_TABLE default\n");

  fwrite((void*)field,sizeof(float),nx1*ny1*nz1,f1);
  fclose(f1);



  sprintf(nome,"%s_out_meta.vtk",args[1]);
  f1=fopen(nome, "w");
    fprintf(f1,"# vtk DataFile Version 2.0\n");
  fprintf(f1,"titolo mio\n");
  fprintf(f1,"BINARY\n");
  fprintf(f1,"DATASET STRUCTURED_POINTS\n");
  fprintf(f1,"DIMENSIONS %i %i %i\n",nx1, ny1, nz1/2);
  fprintf(f1,"ORIGIN %f %f %f\n",xmin, ymin, zmin);
  dx=(xmax-xmin)/(nx1-1);
  dy=(ymax-ymin)/(ny1-1);
  dz=(zmax-zmin)/(nz1-1);
  fprintf(f1,"SPACING %f %f %f\n",dx, dy, dz);
  fprintf(f1,"POINT_DATA %i\n",nx1*ny1*(nz1/2));
  fprintf(f1,"SCALARS Ex float 1\n");
  fprintf(f1,"LOOKUP_TABLE default\n");
  fwrite((void*)field,sizeof(float),nx1*ny1*(nz1/2),f1);
  fclose(f1);

  sprintf(nome,"%s_cut_x",args[1]);
  f1=fopen(nome, "w");
  dx=(xmax-xmin)/(nx1-1);
  dy=(ymax-ymin)/(ny1-1);
  dz=(zmax-zmin)/(nz1-1);
  for(i=0;i<nxloc;i++)
    {
      fprintf(f1,"%f  %f\n",i*dx+xmin,E0*field[i+nx1*(ny1/2)+nx1*ny1*(nz1/2)]);
    }
  
  fclose(f1);
  
}
