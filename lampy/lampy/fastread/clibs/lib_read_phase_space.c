#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include "lib_read_phase_space.h"


void read_phasespace(float *ps_component, int dimensionality, char *file_pointer){ 

    FILE *binary;
    int npart, a;
    int offset, jump;

    binary=fopen(file_pointer,"rb");

    if(binary==NULL){
    	printf("Can't read\n");
    	exit(0);
    }
    jump=2*dimensionality+1;
    offset=0;
    while(fread(&npart,sizeof(int),1,binary)==1){

        for(a=0;a<npart;a++){
        fread(&ps_component[offset],sizeof(float),jump,binary);
        fseek(binary,4,SEEK_CUR);
        offset+=jump;
        }
        
    }

}

int count_particles(int dimensionality, char *file_pointer){

    FILE *binary;
    int total_particles_number, npart, nelements;
    int jump;

    binary=fopen(file_pointer,"rb");

    if(binary==NULL){
    	printf("Can't read\n");
    	exit(0);
    }
    jump=2*dimensionality+2;
    total_particles_number=0;
    while(fread(&npart,sizeof(int),1,binary)==1){
        nelements=jump*npart;

        fseek(binary,sizeof(float)*nelements,SEEK_CUR);


        total_particles_number+=npart;


    }

    return total_particles_number;

}