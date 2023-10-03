#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include "Global.hpp"
#include "Functions.hpp"
#include "InitiationFunctions.hpp"

//////////////
//FUNTIONS //
/////////////

// auxilliary functions //

//rand_Normal, sample a number from a normal distribution with mean mean and s.e. sigma
double rand_Normal (double mean, double sigma){
    double U1, U2, W, mult;
    static double X1, X2;
    static int call = 0;
    
    if (call == 1)
    {
        call = !call;
        return (mean + sigma * (double) X2);
    }
    
    do
    {
        U1 = -1 + ((double) rand () / RAND_MAX) * 2;
        U2 = -1 + ((double) rand () / RAND_MAX) * 2;
        W = pow (U1, 2) + pow (U2, 2);
    }
    while (W >= 1 || W == 0);
    
    mult = sqrt ((-2 * log (W)) / W);
    X1 = U1 * mult;
    X2 = U2 * mult;
    
    call = !call;
    
    return (mean + sigma * (double) X1);
}//rand_Normal

//create random number from 0 to 1
double unifRand()
{
    return rand() / double(RAND_MAX);
}//unifRand

//random integer from rangeLow to rangeHight
int uniform_distribution(int rangeLow, int rangeHigh)
{
    int myRand = (int)rand();
    int range = rangeHigh - rangeLow + 1; //+1 makes it [rangeLow, rangeHigh], inclusive.
    int myRand_scaled = (myRand % range) + rangeLow;
    return myRand_scaled;
}//uniform_distribution

//absolute value of difference between two doubles
double Absolute(double a, double b){
    if(a>b) return(a-b);
    else return(b-a);
}//Absolute

double min(double a, double b){
    if(a>b) return(b);
    else return(a);
}//min

double max(double a, double b){
    if(a<b) return(b);
    else return(a);
}//max


double minOf10(double a1, double a2, double a3, double a4, double a5, double a6, double a7, double a8, double a9, double a10){
    double min = a1;
    
    if(a2<min) min = a2;
    if(a3<min) min = a3;
    if(a4<min) min = a4;
    if(a5<min) min = a5;
    if(a6<min) min = a6;
    if(a7<min) min = a7;
    if(a8<min) min = a8;
    if(a9<min) min = a9;
    if(a10<min) min = a10;
    
    return(min);
}













