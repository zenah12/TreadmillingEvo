#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include "Global.hpp"
#include "Functions.hpp"
#include "InitiationFunctions.hpp"
#include "WriteFiles.hpp"


// stores the filament length over time from time tMin to the end of the simulation
// needs to be adjusted so only the length every dt time is stored 
void  StoreFilamentDynamics(double mat[numt][2],  double tMin, double dt){
	char name[250];
    int j=0;
	FILE *pfile=NULL;
    
    
    sprintf(name, "Filament_%.d.txt",1);
	
	pfile = fopen(name, "w+");
    
    // find i_min
    int i_min = 0;
    while(j<numt && i_min==0){
        if(mat[j][0]>tMin) i_min = j;
        j++;
    }
    
    double timePassed = 0.0;
    int count = 0;
    fprintf(pfile, "%.10f   %.10f\n", mat[i_min][0], mat[i_min][1]);
    for(j=i_min; j<numt; j++){
        timePassed = timePassed + mat[j][0];
        // store entry every dt time itnerval
        if(timePassed>dt){
        //printf("here\n");
            fprintf(pfile, "%.10f   %.10f\n", mat[j][0], mat[j][1]);
            timePassed = 0.0;
            count ++;
        }
        if(count>1000000) break;
    }
    
	fclose(pfile);
}//WriteDelta


void  StoreParametersEvo(double mat[numt][paramsN],  int maxi, int index){
    char name[500];
    int i=0;
    int j=0;
    FILE *pfile=NULL;
    

    sprintf(name, "Evo_targetMean%.2f_win%.2f_mu%.2f_semI%.2f_maxEvo%.d_rMean%.2f_samplingN%.d_%.d.txt", targetMean, win, mutationStep, semIndex, numt, rMean, samplingSteps, index);
    
    pfile = fopen(name, "w+");

    for(i=2; i<(maxi+1); i++){
        for(j=0; j<paramsN; j++){
            fprintf(pfile, "%.10f   ", mat[i][j]);
        }
        fprintf(pfile, "\n");
    }
    
    fclose(pfile);
}//WriteDelta


void  StoreMfilaments(int mat[M][maxN], int index){
    char name[250];
    int i=0;
    int j=0;
    FILE *pfile=NULL;
    
    sprintf(name, "FilamentsM_targetMean%.2f_mutStep%.2f_semI%.2f_maxEvo%.d_mukonTp%.d_mukonTm%.d_mukoffTp%.d_mukonDp%.d_mukonDm%.d_mukoffDp%.d_muwde%.d_muwre%.d_rMean%.2f_%.d.txt", targetMean, mutationStep,semIndex, numt, mu_konTp, mu_konTm, mu_koffTp, mu_konDp, mu_konDm, mu_koffDp, mu_w_de, mu_w_re, rMean, index);
    
    pfile = fopen(name, "w+");
    
    for(i=0; i<M; i++){
        for(j=0; j<maxN; j++){
            fprintf(pfile, "%.d   ", mat[i][j]+2);
        }
        fprintf(pfile, "\n");
    }
    
    fclose(pfile);
}//StoreMfilaments


