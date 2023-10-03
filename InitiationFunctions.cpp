#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include "Global.hpp"
#include "Functions.hpp"
#include "InitiationFunctions.hpp"


///////////////////////////////////
/////////Initiating matirces//////
/////////////////////////////////

// initiate spatial profiels
void Init_numt(double mat[numt]){
    int i=0;
    for (i=0; i<numt; i++) mat[i] = 0.0;
}//Init_Concentration

void Initiate_fil(int mat[maxN]){
    int i=0;
    
    for(i=0; i<maxN; i++) mat[i] = -1; // -1 means no block
    //initiate a filament of size two
    mat[0] = 1;
    mat[maxN-1] = 1;
}

void Initiate_maxN(int mat[maxN]){
    int i=0;
    
    for(i=0; i<maxN; i++) mat[i] = -1; // -1 means no block
}

void Initiate_2x2(double mat[2][2]){
    int i=0;
    int j=0;
    
    for(i=0; i<2; i++)
        for(j=0; j<2; j++) mat[i][j] = 10000.0;
}

void Initiate_ParsEvo(double mat[numt][paramsN]){
    int i=0;
    int j=0;
    
    for(i=0; i<numt; i++)
        for(j=0; j<paramsN; j++) mat[i][j] = 0.0;
}

void Initiate_MmaxN(int mat[numt/saveEnd][maxN]){
    int i=0;
    int j=0;
    
    for(i=0; i<numt/saveEnd; i++)
        for(j=0; j<maxN; j++)
            mat[i][j] = 0;
}

void Initiate_debug(double mat[numt][6]){
    int i=0;
    int j=0;
    
    for(i=0; i<numt; i++)
        for(j=0; j<6; j++)
            mat[i][j] = 0;
}

void InitialConditionsRandomize(){
    // asign random values for wde, wre in (min, max)
    w_de = pow(10.0, unifRand()*(log10(kmax)-log10(kmin)) + log10(kmin));
    w_re = pow(10.0, unifRand()*(log10(kmax)-log10(kmin)) + log10(kmin));
    // ks for D
    kon_D_p = pow(10.0, unifRand()*(log10(kmax)-log10(kmin)) + log10(kmin));
    koff_D_p = pow(10.0, unifRand()*(log10(kmax)-log10(kmin)) + log10(kmin));
    kon_D_m = pow(10.0, unifRand()*(log10(kmax)-log10(kmin)) + log10(kmin));
    koff_D_m = pow(10.0, unifRand()*(log10(kmax)-log10(kmin)) + log10(kmin));
    // ks for T
    kon_T_p = pow(10.0, unifRand()*(log10(kmax)-log10(kmin)) + log10(kmin));
    koff_T_p = pow(10.0, unifRand()*(log10(kmax)-log10(kmin)) + log10(kmin));
    kon_T_m = pow(10.0, unifRand()*(log10(kmax)-log10(kmin)) + log10(kmin));
    koff_T_m = pow(10.0, unifRand()*(log10(kmax)-log10(kmin)) + log10(kmin));
}
