#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <cmath>
#include <iostream>
#include <ctime>
#include <algorithm>
#include <vector>
#include <random>

#include "Global.hpp"
#include "Functions.hpp"
#include "InitiationFunctions.hpp"
#include "WriteFiles.hpp"
#include "Filaments.hpp"
#include "Evolution.hpp"


using namespace std;

// declare parameters
// system parameters
// D - 0 block, T - 1 block
double kon_T_p0 = 0.1;//0.8;
double kon_D_p0 = 0.0;//0.008;
double kon_T_m0 = 0.1;//0.008;//0.00002;
double kon_D_m0 = 0.0;//0.008;//0.00002;
double koff_T_p0 = 0.1;//0.25;
double koff_D_p0 = 0.0;//0.25;//2.5;
double koff_T_m0 = 0.1;//0.0025;
double koff_D_m0 = 0.0;//0.0025;//5.0;
double w_de0 = 1.0;     //T to D
double w_re0 = 0.1;   //D to T
//
double kon_T_p = kon_T_p0;
double kon_D_p = kon_D_p0;
double kon_T_m = kon_T_m0;
double kon_D_m = kon_D_m0;
double koff_T_p = koff_T_p0;
double koff_D_p = koff_D_p0;
double koff_T_m = koff_T_m0;
double koff_D_m = koff_D_m0;
double w_de = w_de0;
double w_re = w_re0;
// max and minimum rates
double kmin = 0.00001;
double kmax = 1.0;
// evolution parameters
double mu = 0.1;
double targetMean = 200.0;
double win = 10.0;
double targetRatio = 10.0;
int evoTimeNow = 1;
double mutationStep = 2.0;//2.0;
int maxEvo = 100;
double semIndex = 0.0;      // 1 if sem to be included in optimization and 0 otherwise
double rMean = 0.01;
// mutation indeces
int mu_konTp = 1;
int mu_konTm = 1;
int mu_koffTp = 1;
int mu_konDp = 1;
int mu_konDm = 1;
int mu_koffDp = 1;
int mu_w_de = 0;
int mu_w_re = 0;
int mu_konkoffT = 1;
int mu_konkoffD = 0;
// initial conditions
int D_bound_0 = 1;
int T_bound_0 = 1;
// concentrations
int D_bound_t = D_bound_0;
int T_bound_t = T_bound_0;

double tMin = 1000000.0;
double tMax = 2000000.0;

// simulation parameters (discretization)
double t_samp = 5.0/w_re;
double dt = 0.1;
double t_now=0.0;

int minus_end=0;
int plus_end=0;
int filamentLength = 0;

int samplingSteps = 20;

//Add objects
Filaments myFil;
Evolution   myEvo;

int indOut = 0;
int tEvo=0;

int SampleSize = 500000;

double maxRatioDelta = 0.0;

//////////////////////////////////////////////////////////////////////////////
////////////////////////// M A I N   P R O G R A M //////////////////////////
////////////////////////////////////////////////////////////////////////////

int main(int argc, char ** argv) {
	
//    srand(time(NULL));
    
    
    
    indOut = atoi(argv[1]);
    double randomize = 1.0*indOut;
    
    //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count() + randomize;
    //mt19937 gen(seed);
    
    //std::random_device rd;
    //std::mt19937 mt(rd());
    
    //mt19937 mt_rand((unsigned)time(NULL)+randomize);
    
    srand((unsigned)time(NULL)+randomize);
    
    //genrand_int32(time(NULL));
    
    // initiate evolution matrices
    Initiate_2x2(myEvo.meanSE);
    Initiate_ParsEvo(myEvo.paramsEvo);
    
    // keep omegas fixed
    w_de = w_de0;
    w_re = w_re0;
    
    // run dynamics once for initial conditions
    ///////////////////////////////////////////////////////////////////////
    // initial conditions to be set manually //
    // initiate filament to be made of one D and one T element
    Initiate_maxN(myFil.filament);
    myFil.filament[0]=1;
    myFil.filament[1]=0;
    // initiate position index matrices for bound elements
    Initiate_maxN(myFil.pos_0);
    Initiate_maxN(myFil.pos_1);
    myFil.pos_1[0] = 0;
    myFil.pos_0[0] = 1;
    D_bound_t = 1;
    T_bound_t = 1;
    // set indexes for plus and minus end positions
    plus_end = 0;
    minus_end = 1;
    // set initial length
    filamentLength = 2;
    
    // initiate matrix that can be used for debugging
    // Initiate_debug(myEvo.debugArray);
    ///////////////////////////////////////////////////////////////////////
    
    // counting ints
    int count=0;
    int count2 = 0;
    int count3 = 0;
    int count4 = 0;
    
    // randomise initial conditions //
    //InitialConditionsRandomize();
    
    // define t_samp for current set of parameters
    t_samp = 2.0/minOf10(kon_T_p, kon_D_p, kon_T_m, kon_D_m, koff_T_p, koff_D_p, koff_T_m, koff_D_m, w_de, w_re);
    
    printf("Initial value: w_re': %.10f  w_de%.10f\n", w_re, w_de);
    // mutate
    //myEvo.mutation();
    // define tMin for reaching SS for current set of parameters
    double tMin = 10000.0;
    tMin = 1000000.0;//1.0/minOf10( kon_T_p, kon_D_p, kon_T_m, kon_D_m, koff_T_p, koff_D_p, koff_T_m, koff_D_m, w_de, w_re);
    tMax =4.0*tMin;
    ////////////////////////////
    // filament dynamics loop //
    ////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    ////////////////////////////
    // filament dynamics loop //
    ////////////////////////////
    t_now=0.0;
    // filament dynamics until late tMin without storing
    while(t_now < tMin){
        // update reaction rates vector
        myFil.update_Ri(myFil.filament, myFil.Ri, filamentLength);
        // sample dt interval
        dt = myFil.sample_dt(myFil.Ri);
       // printf("dt: %.5f\n", dt);
        t_now = t_now+dt;
        // update filament
        myFil.update_f(myFil.filament, myFil.Ri, myFil.pos_0, myFil.pos_1);
     //   printf("time: %.5f\n", t_now);
    }
   
    // filament dynamics from tMin onwards until countMax recordings
    count=0;
    double meanNow=0.0;
    double mean2Now=0.0;
    double meanNow1=0.0;
    double mean2Now1=0.0;
    double meanNow2=0.0;
    double mean2Now2=0.0;
    double meanNow3=0.0;
    double mean2Now3=0.0;
    double meanNow4=0.0;
    double mean2Now4=0.0;
    double tPassed=0.0;
    
    count=0;
    count2 = 0;
    count3 = 0;
    count4 = 0;
    meanNow=0.0;
    mean2Now=0.0;
    meanNow1=0.0;
    mean2Now1=0.0;
    meanNow2=0.0;
    mean2Now2=0.0;
    meanNow3=0.0;
    mean2Now3=0.0;
    meanNow4=0.0;
    mean2Now4=0.0;
    maxRatioDelta = 0.0;
    tPassed=0.0;
    int countMax = 1200000; // number of samples to take average over
    int countAll = 0;       // counter of number of samples taken
    int counterSample = 0;
    while(countAll<countMax){
        // update reaction rates vector
        myFil.update_Ri(myFil.filament, myFil.Ri, filamentLength);
        // sample dt interval
        dt = myFil.sample_dt(myFil.Ri);
        t_now = t_now+dt;
        tPassed = tPassed+dt;
        counterSample = counterSample + 1;
        // update filament
        myFil.update_f(myFil.filament, myFil.Ri, myFil.pos_0, myFil.pos_1);
        // sampling
        // sample every t_samp second
        if(counterSample>samplingSteps){
            // get meanNow1
            if(countAll<=countMax/4.0){
                meanNow1=meanNow1+1.0*filamentLength;
                mean2Now1=mean2Now1+1.0*filamentLength*filamentLength;
                count++;
            }
            // get meanNow2
            else if(countAll<=2.0*countMax/4.0){
                meanNow2=meanNow2+1.0*filamentLength;
                mean2Now2=mean2Now2+1.0*filamentLength*filamentLength;
                count2++;
            }
            // get meanNow3
            else if(countAll<=3.0*countMax/4.0){
                meanNow3=meanNow3+1.0*filamentLength;
                mean2Now3=mean2Now3+1.0*filamentLength*filamentLength;
                count3++;
            }
            else{
                meanNow4=meanNow4+1.0*filamentLength;
                mean2Now4=mean2Now4+1.0*filamentLength*filamentLength;
                count4++;
            }
            counterSample = 0;
            countAll++;
        }
    }// while
    meanNow = (meanNow1/(1.0*count) + meanNow2/(1.0*count2) + meanNow3/(1.0*count3) + meanNow4/(1.0*count4))/4.0;
    mean2Now = (mean2Now1/(1.0*count) + mean2Now2/(1.0*count2) + mean2Now4/(1.0*count3) + mean2Now4/(1.0*count4));
    
    
//    while(t_now<tMax){
//        // update reaction rates vector
//        myFil.update_Ri(myFil.filament, myFil.Ri, filamentLength);
//        // sample dt interval
//        dt = myFil.sample_dt(myFil.Ri);
//        t_now = t_now+dt;
//        tPassed = tPassed+dt;
//        // update filament
//        myFil.update_f(myFil.filament, myFil.Ri, myFil.pos_0, myFil.pos_1);
//        // sampling
//        // sample every  second
//        if(tPassed>t_samp){
//            meanNow=meanNow+1.0*filamentLength;
//            mean2Now=mean2Now+1.0*filamentLength*filamentLength;
//            count++;
//            tPassed=0.0;
//        }
//    }// while
//    meanNow = meanNow/(1.0*count);
//    mean2Now = mean2Now/(1.0*count);
    
    // inittiate myEvo.meanSE
    Initiate_2x2(myEvo.meanSE);
    // store filamena mean length and length s.d. in myEvo.meanSE
    myEvo.updateMeanSE(meanNow, mean2Now, myEvo.meanSE, 0);
    
    ///////////////////////////////////
    /// E V O L U T I O N   L O O P ///
    ///////////////////////////////////
    tEvo=1 ;
    // store initial parameters
    myEvo.updateParamsEvo(myEvo.paramsEvo, myEvo.meanSE[0][0], myEvo.meanSE[0][1]);
    //while(pow(myEvo.meanSE[0][0]/myEvo.meanSE[0][1]-targetRatio, 2.0)>0.5 || tEvo<10){// && myEvo.meanSE[0][0]/myEvo.meanSE[0][1]>1.5){ //20.0){
    while(pow(myEvo.meanSE[0][0]-targetMean, 2.0)>win*win && tEvo < (numt-1)){
        ///////////////////////////////////////////////////////////////////////
        // initial conditions to be set manually //
        // initiate filament to be made of one D and one T element
        Initiate_maxN(myFil.filament);
        myFil.filament[0]=1;
        myFil.filament[1]=0;
        // initiate position index matrices for bound elements
        Initiate_maxN(myFil.pos_0);
        Initiate_maxN(myFil.pos_1);
        myFil.pos_1[0] = 0;
        myFil.pos_0[0] = 1;
        D_bound_t = 1;
        T_bound_t = 1;
        // set indeces for plus and minus end positions
        plus_end = 0;
        minus_end = 1;
        // set initial length
        filamentLength = 2;
        ///////////////////////////////////////////////////////////////////////
    
        // MUTATION
        myEvo.mutation();
        // define t_samp for current set of parameters
        t_samp = 0.1/minOf10(kon_T_p, kon_D_p, kon_T_m, kon_D_m, koff_T_p, koff_D_p, koff_T_m, koff_D_m, w_de, w_re);
        
        ////////////////////////////
        // filament dynamics loop //
        ////////////////////////////
        t_now = 0.0;
        // filament dynamics until late tMin without storing
        while(t_now < tMin){
            // update reaction rates vector
            myFil.update_Ri(myFil.filament, myFil.Ri, filamentLength);
            // sample dt interval
            dt = myFil.sample_dt(myFil.Ri);
            t_now = t_now+dt;
          //  printf("tnow: %.5f\nf", t_now);
            // update filament
            myFil.update_f(myFil.filament, myFil.Ri, myFil.pos_0, myFil.pos_1);
        }
        // filament dynamics for countMax steps
        count=0;
        count2 = 0;
        count3 = 0;
        count4 = 0;
        counterSample = 0;
        meanNow=0.0;
        mean2Now=0.0;
        meanNow1=0.0;
        mean2Now1=0.0;
        meanNow2=0.0;
        mean2Now2=0.0;
        meanNow3=0.0;
        mean2Now3=0.0;
        meanNow4=0.0;
        mean2Now4=0.0;
        maxRatioDelta = 0.0;
        tPassed=0.0;
        countMax = 1200000; // number of samples to take average over
        countAll = 0;       // counter of number of samples taken
        
        while(countAll<countMax){
            // update reaction rates vector
            myFil.update_Ri(myFil.filament, myFil.Ri, filamentLength);
            // sample dt interval
            dt = myFil.sample_dt(myFil.Ri);
            t_now = t_now+dt;
            tPassed = tPassed+dt;
            counterSample = counterSample + 1;
            // update filament
            myFil.update_f(myFil.filament, myFil.Ri, myFil.pos_0, myFil.pos_1);
            // sampling
            if(counterSample>samplingSteps){
                // get meanNow1
                if(countAll<=countMax/4){
                    meanNow1=meanNow1+1.0*filamentLength;
                    mean2Now1=mean2Now1+1.0*filamentLength*filamentLength;
                    count++;
                }
                // get meanNow2
                else if(countAll<=2.0*countMax/4.0){
                    meanNow2=meanNow2+1.0*filamentLength;
                    mean2Now2=mean2Now2+1.0*filamentLength*filamentLength;
                    count2++;
                }
                // get meanNow3
                else if(countAll<=3.0*countMax/4.0){
                    meanNow3=meanNow3+1.0*filamentLength;
                    mean2Now3=mean2Now3+1.0*filamentLength*filamentLength;
                    count3++;
                }
                else{
                    meanNow4=meanNow4+1.0*filamentLength;
                    mean2Now4=mean2Now4+1.0*filamentLength*filamentLength;
                    count4++;
                }
                counterSample = 0;
                countAll++;
            }
        }// while
        meanNow = (meanNow1/(1.0*count) + meanNow2/(1.0*count2) + meanNow3/(1.0*count3) + meanNow4/(1.0*count4))/4.0;
        mean2Now = (mean2Now1/(1.0*count) + mean2Now2/(1.0*count2) + mean2Now4/(1.0*count3) + mean2Now4/(1.0*count4))/4.0;
        // filament dynamics loop end
        
        // get mean filament length and variance in a distirbution sampled after 10^6 steps
        myEvo.updateMeanSE(meanNow, mean2Now, myEvo.meanSE, 1);
        // if evoStep > 0 select
        maxRatioDelta = myEvo.getMaxRatio(meanNow1/(1.0*count), meanNow2/(1.0*count2), meanNow3/(1.0*count3), meanNow4/(1.0*count4));
        
        myEvo.selection(myEvo.meanSE, myEvo.paramsEvo, maxRatioDelta);
        
        StoreParametersEvo(myEvo.paramsEvo, evoTimeNow, indOut);
        
        tEvo++;
       // if(tEvo > (numt-1)) break;
    }

    //StoreFilamentDynamics(myFil.FilamentLengthOverTime, tMin, t_samp);
    StoreParametersEvo(myEvo.paramsEvo, evoTimeNow, indOut);
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // repeat analysis for filament with final selected rates to export state of filament at M time points //
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

    // comment below if no need for filament composition sample //
    Initiate_MmaxN(myFil.filamentStates);
    // initiate filament to be made of one D and one T element
    Initiate_maxN(myFil.filament);
    myFil.filament[0]=1;
    myFil.filament[1]=0;
    // initiate position index matrices for bound elements
    Initiate_maxN(myFil.pos_0);
    Initiate_maxN(myFil.pos_1);
    myFil.pos_1[0] = 0;
    myFil.pos_0[0] = 1;
    D_bound_t = 1;
    T_bound_t = 1;
    // set indeces for plus and minus end positions
    plus_end = 0;
    minus_end = 1;
    // set initial length
    filamentLength = 2;
    ///////////////////////////////////////////////////////////////////////
    // define t_samp for current set of parameters
    t_samp = 5.0/minOf10( kon_T_p, kon_D_p, kon_T_m, kon_D_m, koff_T_p, koff_D_p, koff_T_m, koff_D_m, w_de, w_re);
    ////////////////////////////
    // filament dynamics loop //
    ////////////////////////////
    // filament dynamics until late tMin without storing
    t_now=0.0;
    while(t_now < tMin){
        // update reaction rates vector
        myFil.update_Ri(myFil.filament, myFil.Ri, filamentLength);
        // sample dt interval
        dt = myFil.sample_dt(myFil.Ri);
        t_now = t_now+dt;
        //printf("tnow: %.5f\nf", t_now);
        // update filament
        myFil.update_f(myFil.filament, myFil.Ri, myFil.pos_0, myFil.pos_1);
    }
    int indexCount=0;
    int indexSave = 0;
    count =0;
    // filament dynamics from tMin to tMax with registration
    while(count<numt){// && t_now<tMax){
        // update reaction rates vector
        myFil.update_Ri(myFil.filament, myFil.Ri, filamentLength);
        // sample dt interval
        dt = myFil.sample_dt(myFil.Ri);
        t_now = t_now+dt;
        // update filament
        myFil.update_f(myFil.filament, myFil.Ri, myFil.pos_0, myFil.pos_1);
        count++;
        if(indexSave==saveEnd){
            myEvo.storeCurrentFilament(myFil.filamentStates, myFil.filament, indexCount);
            indexCount++;
            indexSave =0;
        }
        indexSave++;
    }// while
    StoreMfilaments(myFil.filamentStates, indOut);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////

	return 0;
}



