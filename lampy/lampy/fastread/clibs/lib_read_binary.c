#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <errno.h>
#include "lib_read_binary.h"


void read_binary(float *field, float *x, float*y, float *z, char *file_pointer){ 

	FILE *binary;
	int a, b, c, totlocalgrid, nx, ny, nz, i, j, k, ii;
	int npx, npy, npz, nproc_y, nproc_z, nproc_x, nparam;
	int offset, offsetx, offsety, offsetz;
	int *integer_params;
	float *real_params, ***field_temp;
	errno = 0;
	
	binary = fopen(file_pointer, "rb");
	
	if(binary == NULL){
		printf("Can't read %s with error %d \n", file_pointer, errno);
		exit(0);
	}
	
	fseek(binary, 4, 1);
	fread(&nparam, sizeof(int), 1, binary);
	fseek(binary, 8, SEEK_CUR);
	
	integer_params = (int *) malloc(nparam*sizeof(int));
	real_params = (float *) malloc(nparam*sizeof(float));
	
	fread(integer_params, sizeof(int), nparam,binary);
	fseek(binary, 8, SEEK_CUR);
	fread(real_params, sizeof(float), nparam, binary);
	fseek(binary, 4, SEEK_CUR);
	
	nx = integer_params[3];
	ny = integer_params[4];
	nz = integer_params[6];
	nproc_y = integer_params[0];
	nproc_z = integer_params[1];
	nproc_x = integer_params[2];
	
	field_temp = (float ***) malloc(nz*sizeof(float **));
	for(k = 0; k < nz; k++){
		field_temp[k] = (float **) malloc(ny*sizeof(float *));
			for(j = 0; j < ny; j++){
				field_temp[k][j] = (float *) malloc(nx*sizeof(float ));
			}
	
	}
	
	
	offset = 0;
	offsetx=0;
	for(c = 0; c < nproc_x; c++){
	
		offsetz = 0;
		for(a = 0; a < nproc_z; a++){
		
			offsety = 0;
			for(b = 0; b < nproc_y; b++){
			
				fseek(binary, 4, SEEK_CUR);
				fread(&npx, sizeof(int), 1, binary);
				fread(&npy, sizeof(int), 1, binary);
				fread(&npz, sizeof(int), 1, binary);
				totlocalgrid = npx*npy*npz;
				fseek(binary, 8, SEEK_CUR);
				fread(&field[offset], sizeof(float), totlocalgrid, binary);
				ii = 0;
	
				for(k = 0; k < npz; k++){
					for(j = 0; j < npy; j++){
						for(i = 0; i < npx; i++){
							field_temp[k + offsetz][j + offsety][i + offsetx] = field[ii+offset];
							ii += 1;
						}
					}
				}		
				offset+=totlocalgrid;
				fseek(binary,4,SEEK_CUR);
				offsety+=npy;
			}
			offsetz+=npz;
		}
		offsetx+=npx;
	}
	
	ii = 0;
	for(k = 0; k < nz; k++){
		for(j = 0; j < ny; j++){
			for(i = 0; i < nx; i++){
				field[ii++] = field_temp[k][j][i];
			}
		}
	}
	
	fseek(binary, 4, SEEK_CUR);
	fread(x, sizeof(float), nx, binary);
	fseek(binary, 8, SEEK_CUR);
	fread(y, sizeof(float), ny, binary);
	fseek(binary, 8, SEEK_CUR);
	fread(z, sizeof(float), nz, binary);
	
	fclose(binary);

}
