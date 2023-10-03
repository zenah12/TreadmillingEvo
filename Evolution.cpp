#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include <ctime>
#include <algorithm>
#include <vector>
#include "Global.hpp"
#include "Functions.hpp"
#include "InitiationFunctions.hpp"
#include "Evolution.hpp"

using namespace std;


/////////////////////////////////////////////////////////////////////
//Here all the functions used in the object Evolution are defined //
////////////////////////////////////////////////////////////////////

/// LOCAL FUNCTIONS ///
double dEnergy(double meansSEs[2][2]){
    double Enow = 0.0;
    double Ebefore = 0.0;
        
    //Ebefore = sqrt((meansSEs[0][0]/meansSEs[0][1] - targetRatio)*(meansSEs[0][0]/meansSEs[0][1] - targetRatio));
    //Enow = sqrt((meansSEs[1][0]/meansSEs[1][1] - targetRatio)*(meansSEs[1][0]/meansSEs[1][1] - targetRatio));
    Ebefore = sqrt((meansSEs[0][0] - targetMean)*(meansSEs[0][0] - targetMean));//+ semIndex*meansSEs[0][0]/meansSEs[0][1];
    Enow = sqrt((meansSEs[1][0] - targetMean)*(meansSEs[1][0] - targetMean)); //+ semIndex*meansSEs[1][0]/meansSEs[1][1];
    
    //printf("se0: %.5f    se1: %.5f     energy before: %.5f     energy now: %.5f\n\n", meansSEs[0][1], meansSEs[1][1], Ebefore, Enow);
    
    return(Enow/Ebefore);
}

/// GLOBAL FUNCTIONS ///

// decide whether to keep current parameter set or revert to the one prior to mutation according
// to the minimization of the energy function
void Evolution::selection(double meansSEs[2][2], double parsEvo[numt][paramsN], double ratio){
    double dE=0.0;
    
    dE = dEnergy(meansSEs);
    
//    printf("Selection step. evoTimeNow:%.d  mean1:%.4f  mean2:%.4f  dE:%.4f ratio:%.4f \n", evoTimeNow, meansSEs[0][0], meansSEs[1][0], dE, ratio);

    // if energy is lower for mew parameters keep current parameters and update evolution matrix
    //if(dE < 0.98 && meansSEs[1][0]<maxN*0.1 && ratio<1.0+rMean && ratio>1.0-rMean){// additional contraint added to avoid infinite filament
    
    if(dE < 0.98 && ratio<1+rMean && ratio>1-rMean){// additional contraint added to avoid infinite filament

        // update meansSEs
        meansSEs[0][0] = meansSEs[1][0]; //selected becomes current
        meansSEs[0][1] = meansSEs[1][1]; //selected becomes current
        // update parameters
        updateParamsEvo(parsEvo, meansSEs[0][0], meansSEs[0][1]);

        
        printf("selected! w_re': %.7f   w_de:%.10f  kon_D_p:%.10f   kon_D_m:%.10f   koff_D_p:%.10f  koff_D_m:%.10f  mean:%.5f   s.d.:%.5f\n", w_re, w_de, kon_D_p, kon_D_m, koff_D_p, koff_D_m, meansSEs[0][0], meansSEs[0][1]);
        
        printf("selected! w_re': %.7f   w_de:%.10f  kon_T_p:%.10f   kon_T_m:%.10f   koff_T_p:%.10f  koff_T_m:%.10f  mean:%.10f  s.d.:%.5f\n\n", w_re, w_de, kon_T_p, kon_T_m, koff_T_p, koff_T_m, meansSEs[0][0], meansSEs[0][1]);
    }
    // else revert to previous
    else{
        kon_T_p = parsEvo[evoTimeNow][0];
        kon_D_p = parsEvo[evoTimeNow][1];
        kon_T_m = parsEvo[evoTimeNow][2];
        kon_D_m = parsEvo[evoTimeNow][3];
        koff_T_p = parsEvo[evoTimeNow][4];
        koff_D_p = parsEvo[evoTimeNow][5];
        koff_T_m = parsEvo[evoTimeNow][6];
        koff_D_m = parsEvo[evoTimeNow][7];
        w_de = parsEvo[evoTimeNow][8];
        w_re = parsEvo[evoTimeNow][9];
    }
    
}

void Evolution::mutation(){
    double factor = 0.0;
    double konkoffTnew = 0.0;
    double konkoffDnew = 0.0;
    
    // mutate w rates //
//    if(mu_w_re == 1){
//        factor = pow(10.0, (mutationStep - 2.0*mutationStep*unifRand()));
//        if(factor*w_re<kmax && factor*w_re>kmin) w_re = factor*w_re;
//    }
//    if(mu_w_de == 1){
//        factor = pow(10.0, (mutationStep - 2.0*mutationStep*unifRand()));
//        if(factor*w_de<kmax && factor*w_de>kmin) w_de = factor*w_de;
//    }
    ////////////////////
    // mutate k rates //
    ///////////////////
    // START WITH T RATES//
    // ratio change //
    if(mu_konkoffT == 1){
        // mutate ratio konkoffT
        factor = pow(10.0, (mutationStep - 2.0*mutationStep*unifRand()));
        if ((kon_T_p/koff_T_p)*factor>kmin/kmax && (kon_T_p/koff_T_p)*factor<kmax/kmin) konkoffTnew = (kon_T_p/koff_T_p)*factor;
        else konkoffTnew = (kon_T_p/koff_T_p);
        
        // Plus side //
        // randomly choose whether to mutate konTp or koffTp
        if(unifRand()<0.5){ //konTp
            factor = pow(10.0, (mutationStep - 2.0*mutationStep*unifRand()));
            while(factor*kon_T_p>kmax || factor*kon_T_p<kmin || factor*kon_T_p/konkoffTnew>kmax || factor*kon_T_p/konkoffTnew<kmin){
                factor = pow(10.0, (mutationStep - 2.0*mutationStep*unifRand()));
//                printf("here 1, factor: %.4f\n", factor);
            }//while
            kon_T_p = factor*kon_T_p;
            koff_T_p = kon_T_p/konkoffTnew;
        }//if
        else{//koffTp
            factor = pow(10.0, (mutationStep - 2.0*mutationStep*unifRand()));
            while(factor*koff_T_p>kmax || factor*koff_T_p<kmin || factor*koff_T_p*konkoffTnew>kmax || factor*koff_T_p*konkoffTnew<kmin){
                factor = pow(10.0, (mutationStep - 2.0*mutationStep*unifRand()));
//                printf("here 2, factor: %.4f\n", factor);
            }// while
            koff_T_p = factor*koff_T_p;
            kon_T_p = koff_T_p*konkoffTnew;
        }//else
        // Minus side //
        // randomly choose whether to mutate konTm or koffTm
        if(unifRand()<0.5){ //konTm
            factor = pow(10.0, (mutationStep - 2.0*mutationStep*unifRand()));
            while(factor*kon_T_m>kmax || factor*kon_T_m<kmin || factor*kon_T_m/konkoffTnew>kmax || factor*kon_T_m/konkoffTnew<kmin){
                factor = pow(10.0, (mutationStep - 2.0*mutationStep*unifRand()));
//                printf("here 3, factor: %.4f\n", factor);
            }//while
            kon_T_m = factor*kon_T_m;
            koff_T_m = kon_T_m/konkoffTnew;
        }//if
        else{//koffTm
            factor = pow(10.0, (mutationStep - 2.0*mutationStep*unifRand()));
            while(factor*koff_T_m>kmax || factor*koff_T_m<kmin || factor*koff_T_m*konkoffTnew>kmax || factor*koff_T_m*konkoffTnew<kmin){
                factor = pow(10.0, (mutationStep - 2.0*mutationStep*unifRand()));
//                printf("here 4, factor: %.4f\n", factor);
            }//while
            koff_T_m = factor*koff_T_m;
            kon_T_m = koff_T_m*konkoffTnew;
        }//else
    }

    // REPEAT FOR D RATES//
    // ratio change //
    if(mu_konkoffD == 1){
        factor = pow(10.0, (mutationStep - 2.0*mutationStep*unifRand()));
        if ((kon_D_p/koff_D_p)*factor>kmin/kmax && (kon_D_p/koff_D_p)*factor<kmax/kmin)  konkoffDnew = (kon_D_p/koff_D_p)*factor;
        else konkoffDnew = kon_D_p/koff_D_p;
    
        // Plus side //
        // randomly choose whether to mutate konDp or koffDp
        if(unifRand()<0.5){ //konDp
            factor = pow(10.0, (mutationStep - 2.0*mutationStep*unifRand()));
            while(factor*kon_D_p>kmax || factor*kon_D_p<kmin || factor*kon_D_p/konkoffDnew>kmax || factor*kon_D_p/konkoffDnew<kmin){
                factor = pow(10.0, (mutationStep - 2.0*mutationStep*unifRand()));
//                printf("here 5, factor: %.4f     konkoffDnew: %.8f \n", factor, konkoffDnew);
            }//while
            kon_D_p = factor*kon_D_p;
            koff_D_p = kon_D_p/konkoffDnew;
        }//if
        else{//koffDp
            factor = pow(10.0, (mutationStep - 2.0*mutationStep*unifRand()));
            while(factor*koff_D_p>kmax || factor*koff_D_p<kmin || factor*koff_D_p*konkoffDnew>kmax || factor*koff_D_p*konkoffDnew<kmin){
                factor = pow(10.0, (mutationStep - 2.0*mutationStep*unifRand()));
//                printf("here 6, factor: %.4f     konkoffDnew: %.8f \n", factor, konkoffDnew);
            }// while
            koff_D_p = factor*koff_D_p;
            kon_D_p = koff_D_p*konkoffDnew;
        }//else
        // Minus side //
        // randomly choose whether to mutate konDm or koffDm
        if(unifRand()<0.5){ //konTm
            factor = pow(10.0, (mutationStep - 2.0*mutationStep*unifRand()));
            while(factor*kon_D_m>kmax || factor*kon_D_m<kmin || factor*kon_D_m/konkoffDnew>kmax || factor*kon_D_m/konkoffDnew<kmin){
                factor = pow(10.0, (mutationStep - 2.0*mutationStep*unifRand()));
//               printf("here 7, factor: %.4f     konkoffDnew: %.8f \n", factor, konkoffDnew);
            }//while
            kon_D_m = factor*kon_D_m;
            koff_D_m = kon_D_m/konkoffDnew;
        }//if
        else{//koffTm
            factor = pow(10.0, (mutationStep - 2.0*mutationStep*unifRand()));
            while(factor*koff_D_m>kmax || factor*koff_D_m<kmin || factor*koff_D_m*konkoffDnew>kmax || factor*koff_D_m*konkoffDnew<kmin){
                factor = pow(10.0, (mutationStep - 2.0*mutationStep*unifRand()));
//               printf("here 8, factor: %.4f     konkoffDnew: %.8f \n", factor, konkoffDnew);
            }//while
            koff_D_m = factor*koff_D_m;
            kon_D_m = koff_D_m*konkoffDnew;
        }//else
    }
}


void Evolution::updateMeanSE(double mean_now, double mean2_now, double meanSE[2][2], int index){

    // assign mean and SE
    switch (index) {
        case 0:
            meanSE[0][0] = mean_now;
            meanSE[0][1] = sqrt(mean2_now - mean_now*mean_now);
            break;
        default:
            // if after the first step, current mean and se saved in entry [1]
            meanSE[1][0] = mean_now;
            meanSE[1][1] = sqrt(mean2_now - mean_now*mean_now);
            break;
    }
    
   // printf("mean_0:%.10f     mean_1:%.10f\n\n", meanSE[0][0],meanSE[1][0]);
}//updateMeanSE

// This function is called when the current parameter set is selected for and it updates the evolution step
// and the accepted parameters matrix
void Evolution::updateParamsEvo(double parsEvo[numt][paramsN], double meanNow, double seNow){
    
    evoTimeNow = evoTimeNow + 1;
    printf("evoTimeNow: %.d\n", evoTimeNow);
    
    parsEvo[evoTimeNow][0] = kon_T_p;
    parsEvo[evoTimeNow][1] = kon_D_p;
    parsEvo[evoTimeNow][2] = kon_T_m;
    parsEvo[evoTimeNow][3] = kon_D_m;
    parsEvo[evoTimeNow][4] = koff_T_p;
    parsEvo[evoTimeNow][5] = koff_D_p;
    parsEvo[evoTimeNow][6] = koff_T_m;
    parsEvo[evoTimeNow][7] = koff_D_m;
    parsEvo[evoTimeNow][8] = w_de;
    parsEvo[evoTimeNow][9] = w_re;
    parsEvo[evoTimeNow][10] = meanNow;
    parsEvo[evoTimeNow][11] = seNow;
    parsEvo[evoTimeNow][12] = tEvo*1.0;
    
    printf("tEvo: %.f  \n", tEvo*1.0);
}

void Evolution::storeCurrentFilament(int filamentStates[M][maxN], int filament[maxN], int count){
    int i=0;
    
    for (i=0; i<maxN; i++) {
        filamentStates[count][i] = filament[i];
    }
    
}



double Evolution::getMaxRatio(double mean1, double mean2, double mean3, double mean4){
    double R[6];
    double max=0.0;
    int i=0;
    
    for(i=0; i<6; i++) R[i]=0.0;
    
    R[0] = Absolute(mean1/mean2, 1.0);
    R[1] = Absolute(mean1/mean3, 1.0);
    R[2] = Absolute(mean1/mean4, 1.0);
    R[3] = Absolute(mean2/mean3, 1.0);
    R[4]= Absolute(mean2/mean4, 1.0);
    R[5] = Absolute(mean3/mean4, 1.0);
    
    for (i=0; i<6; i++) if (R[i]>max) max = 1.0 + R[i];
    
    return(max);
}

