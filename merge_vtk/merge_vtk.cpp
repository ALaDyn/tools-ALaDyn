
#include<cstdio>
#include<iostream>
#include<cstdlib>
#include<cmath>
#include<cstring>
#include<string>
#include<fstream>
#include<sstream>
#include<iomanip>
#include<cstdarg>


// #define ENABLE_DEBUG


#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))
#define TRUE 1
#define FALSE 0

int main(int argc, char *argv[])
{
  FILE *inputfile[3], *outputfile;
  int nx1[3], ny1[3], nz1[3], npoints[3];
  float dx[3], dy[3], dz[3], xmin[3], ymin[3], zmin[3];
  //	float x,y,z;
  float E[3];
  char nomefile[100], str[3][600], astr[3][600];

  printf("VectorVTK v1\n");
  fflush(stdout);

  for (int f = 0; f < 3; f++)
  {
    inputfile[f] = fopen(argv[1 + f], "r");
  }
  sprintf(nomefile, "%s_new.vtk", argv[1]);
  outputfile = fopen(nomefile, "w");
  printf("hello! 2\n");
  fflush(stdout);



  for (int f = 0; f < 3; f++)
  {
    printf("file_in[%i]\n", f);
    fscanf(inputfile[f], "%[^\n]\n", str[f]);  //"# vtk DataFile Version 2.0\n"
    printf("%s\n", str[f]);
    fscanf(inputfile[f], "%[^\n]\n", astr[f]);  //titolo mio\n
    printf("%s\n", astr[f]);
    fscanf(inputfile[f], "%[^\n]\n", str[f]);  //BINARY 
    printf("%s\n", str[f]);
    fscanf(inputfile[f], "%[^\n]\n", str[f]);  //DATASET STRUCTURED_POINTS\n
    printf("%s\n", str[f]);
    fscanf(inputfile[f], "%s %i %i %i\n", str[f], &nx1[f], &ny1[f], &nz1[f]);   //DIMENSIONS %i %i %i\n
    fscanf(inputfile[f], "%s %f %f %f\n", str[f], &xmin[f], &ymin[f], &zmin[f]);  //ORIGIN %f %f %f\n
    fscanf(inputfile[f], "%s %f %f %f\n", str[f], &dx[f], &dy[f], &dz[f]);       //SPACING %f %f %f\n
    fscanf(inputfile[f], "%s %i\n", str[f], &npoints[f]);                        //POINT_DATA %i\n
    fscanf(inputfile[f], "%[^\n]\n", str[f]);                                    //SCALARS %s float 1
    fscanf(inputfile[f], "%[^\n]\n", str[f]);                                    //LOOKUP_TABLE default
  }
  printf("hello! 3\n");
  fflush(stdout);

  fprintf(outputfile, "# vtk DataFile Version 2.0\n");
  fprintf(outputfile, "titolo mio\n");
  fprintf(outputfile, "BINARY\n");
  fprintf(outputfile, "DATASET STRUCTURED_POINTS\n");
  fprintf(outputfile, "DIMENSIONS %i %i %i\n", nx1[0], ny1[0], nz1[0]);
  fprintf(outputfile, "ORIGIN %f %f %f\n", xmin[0], ymin[0], zmin[0]);
  fprintf(outputfile, "SPACING %f %f %f\n", dx[0], dy[0], dz[0]);
  fprintf(outputfile, "POINT_DATA %i\n", npoints[0]);
  fprintf(outputfile, "VECTORS E float\n");
  fprintf(outputfile, "LOOKUP_TABLE default\n");
  printf("hello! 4\n");
  fflush(stdout);

  printf("# vtk DataFile Version 2.0\n");
  printf("titolo mio\n");
  printf("BINARY\n");
  printf("DATASET STRUCTURED_POINTS\n");
  printf("DIMENSIONS %i %i %i\n", nx1[0], ny1[0], nz1[0]);
  printf("ORIGIN %f %f %f\n", xmin[0], ymin[0], zmin[0]);
  printf("SPACING %f %f %f\n", dx[0], dy[0], dz[0]);
  printf("POINT_DATA %i\n", npoints[0]);
  printf("VECTORS E float\n");
  printf("LOOKUP_TABLE default\n");



  for (int n = 0; n < npoints[0]; n++)
  {
    for (int f = 0; f < 3; f++)
      fread(&E[f], sizeof(float), 1, inputfile[f]);
    fwrite((void*)E, sizeof(float), 3, outputfile);
  }
  fclose(outputfile);
  for (int f = 0; f < 3; f++) fclose(inputfile[f]);
}


