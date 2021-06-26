#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "local_pluto.h"

#define DATAWIDTH  4
#define DATAHEIGHT 2000

float array[DATAHEIGHT][DATAWIDTH];


void readICfile(char* filename){

    FILE* fp;
    int height, width, ii, jj;
	
	print("Reading IC from: %s\n",filename);
    if((fp = fopen(filename, "r")) == NULL)
        exit(1);

    for(jj=0; jj<DATAHEIGHT; jj++)
        for(ii=0; ii<DATAWIDTH; ii++)
            if(fscanf(fp, "%e", &array[jj][ii])!=1) exit(1);
                
    fclose(fp);
/*
    for(jj=0; jj<DATAHEIGHT; jj++){
        for(ii=0; ii<DATAWIDTH; ii++)
            printf ("%e    ", array[jj][ii]);
        printf("\n");
    }
*/
	//return array;
}

double* initdpv(double r){
	
	int height, width, ii, jj;
	static double  dpv[3];
	
    for(jj=1; jj<DATAHEIGHT; jj++){
        if ((r-array[jj][0])<0 ){ 
          if (abs(r-array[jj-1][0])>(r-array[jj][0])){ dpv[0]=array[jj][1];dpv[1]=array[jj][2];dpv[2]=array[jj][3]; return dpv;}
          else { dpv[0]=array[jj-1][1];dpv[1]=array[jj-1][2];dpv[2]=array[jj-1][3]; return dpv;}
        }
    }
	dpv[0]=array[DATAHEIGHT-1][1];dpv[1]=array[DATAHEIGHT-1][2];dpv[2]=array[DATAHEIGHT-1][3];
    return dpv;
}
