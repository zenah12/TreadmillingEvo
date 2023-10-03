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
#include "Filaments.hpp"

using namespace std;

// local functions //
void cummulativeR(double R[numR], double cumul[numR], double R_T){
    int i = 0;
    int j = 0;
    double tot = 0.0;
    
    for(i=0; i<(numR); i++){
        tot = 0.0;
        for(j=0; j<(i+1); j++) tot = tot + R[j];
        cumul[i] = tot/R_T;
       // printf("i:%.d    Ci:%.5f\n", i, cumul[i]);
    }
}//cummulativeR

// This function goes through the matrix positions and finds the entry that is equal to j
// It then substitutes positions[j] with positions[total] and sets positions[total] to -1
void removeEntry_j(int positions[maxN], int j, int total, int INDEX){
 
    int i = 0;
    int k = 0;
    
    while(positions[i]!=j){
        i++;
    }
    positions[i] = positions[total];
    positions[total] = -1;
}


// update reaction terms
void Filaments::update_Ri(int fil[maxN], double Ri[numR], int length){
    
    // Ri[0]: Db->Tb, Ri[1]: Tb->Db,
    
    Ri[0] = (1.0*D_bound_t)*w_re;
    Ri[1] = (1.0*T_bound_t)*w_de;
    // binding process
    Ri[2] = kon_D_m;
    Ri[3] = kon_D_p;
    Ri[4] = kon_T_m;
    Ri[5] = kon_T_p;
    // if filament is larger than maxN-2 set bidning rates to 0
    if(length>maxN-5){
        Ri[2] = 0.0;
        Ri[3] = 0.0;
        Ri[4] = 0.0;
        Ri[5] = 0.0;
    }
    // unbinding process at minus end
    if(fil[minus_end] == 0){
        Ri[6] = koff_D_m;
        Ri[8] = 0.0;
    }
    else{
        Ri[6] = 0.0;
        Ri[8] = koff_T_m;
    }
    // unbinding process at plus end
    if(fil[plus_end] == 0){
        Ri[7] = koff_D_p;
        Ri[9] = 0.0;
    }
    else{
        Ri[7] = 0.0;
        Ri[9] = koff_T_p;
    }
    // if filament length below 3 set unbinding rates to 0
    if(length<3){
        Ri[6] = 0.0;
        Ri[7] = 0.0;
        Ri[8] = 0.0;
        Ri[9] = 0.0;
    }
}

double Filaments::sample_dt(double Ri[numR]){
    double r = 0.0;
    double R_T = 0.0;
    int j = 0;
    
    for (j=0; j<numR; j++) R_T = R_T + Ri[j];
    
    r = unifRand();
   // printf("RT:%.15f\n", R_T);
    
    return (-log(r)/R_T);
}

// D (NDP bound) state is 0, T (NTP bound) state is 1. A -1 in the filament indicates and empty spot
// pos0 holds index of where 0s lie, pos1 holds index of where 1s lie
void Filaments::update_f(int fil[maxN], double Ri[numR], int pos0[maxN], int pos1[maxN]){
    double r = unifRand(); // sample random number from 0 to 1
    int N_r = 0;
    int i=0;
    double R_T = 0.0;
    double cum_R[numR];
    
    for(i=0; i<numR; i++) R_T = R_T+Ri[i]; // get total reaction
    cummulativeR(Ri, cum_R, R_T);
    
    
    if(r<=cum_R[0] && cum_R[0]!=0.0){
       // printf("1. Here\n");
        // randomly choose Db and turn it to Tb
        N_r = uniform_distribution(0, D_bound_t-1);
       // printf("Dbound %.d N_r:%.d pos0:%.d    entry:%.d\n\n", D_bound_t, N_r, pos0[N_r], fil[pos0[N_r]]);
        fil[pos0[N_r]] = 1;
        // update pos0 and pos1
        pos1[int(T_bound_t)] = pos0[N_r];
        pos0[N_r] = pos0[int(D_bound_t)-1];
        pos0[int(D_bound_t)-1] = -1;
        // update concentrations
        D_bound_t = D_bound_t - 1;
        T_bound_t = T_bound_t + 1;
    }
    else if (r<=cum_R[1] && cum_R[1]!=0.0){
//        // randomly choose Tb and turn it to Db
        N_r = uniform_distribution(0, T_bound_t-1);
        fil[pos1[N_r]] = 0;
        // update pos0 and pos1
        pos0[int(D_bound_t)] = pos1[N_r];
        pos1[N_r] = pos1[int(T_bound_t)-1];
        pos1[int(T_bound_t)-1] = -1;
        // update concentrations
        D_bound_t = D_bound_t + 1;
        T_bound_t = T_bound_t - 1;
    }
    ///////////////////
    // kon processes //
    //////////////////
    else if (r<=cum_R[2] && cum_R[2]!=0.0){
//        printf("3. Here\n");
        // update - end position
        minus_end = minus_end + 1;
        if(minus_end>maxN-1) minus_end = 0;
        // turn an unbound to a bound D at - end
        fil[minus_end] = 0;
        // add to bound D (0) position index
        pos0[int(D_bound_t)] = minus_end;
        // update concentrations of bound D
        D_bound_t = D_bound_t + 1;
        // update filament length
        filamentLength = filamentLength + 1;
    }
    else if (r<=cum_R[3] && cum_R[3]!=0.0){
//        printf("4. Here\n");
        // update + end position
        plus_end = plus_end - 1;
        if(plus_end<0) plus_end = maxN-1;
        // turn an unbound to a bound D at + end
        fil[plus_end] = 0;
        // add to bound D (0) position index
        pos0[int(D_bound_t)] = plus_end;
        // update concentrations of bound D
        D_bound_t = D_bound_t + 1;
        // update filament length
        filamentLength = filamentLength + 1;
    }
    else if (r<=cum_R[4] && cum_R[4]!=0.0){
//        printf("5. Here\n");
        // update - end position
        minus_end = minus_end + 1;
        if(minus_end>maxN-1) minus_end = 0;
        // turn an unbound to a bound T at - end
        fil[minus_end] = 1;
        // add to bound T (1) position index
        pos1[int(T_bound_t)] = minus_end;
        // update concentrations of bound/unbound T
        T_bound_t = T_bound_t + 1;
        // update filament length
        filamentLength = filamentLength + 1;
    }
    else if (r<=cum_R[5] && cum_R[5]!=0.0){
//        printf("6. Here\n");
        // update + end position
        plus_end = plus_end - 1;
        if(plus_end<0) plus_end = maxN-1;
        // turn an unbound to a bound T at + end
        fil[plus_end] = 1;
        // add to bound D (0) position index
        pos1[int(T_bound_t)] = plus_end;
        // update concentrations of bound/unbound T
        T_bound_t = T_bound_t + 1;
        // update filament length
        filamentLength = filamentLength + 1;
    }
    ///////////////////
    // koff processes //
    //////////////////
    else if (r<=cum_R[6] && cum_R[6]!=0.0){
       // printf("7. Here filament at minus end: %.d\n", fil[minus_end]);
        // turn bound to an unbound D at - end
        fil[minus_end] = -1;
        // remove minus_end index from bound D (0) position index
        removeEntry_j(pos0, minus_end, D_bound_t-1, 7);
        // update concentrations of bound/unbound D
        D_bound_t = D_bound_t - 1;
        // update - end position
        minus_end = minus_end - 1;
        if(minus_end<0) minus_end = maxN-1;
        //update filament length
        filamentLength = filamentLength - 1;
    }
    else if (r<=cum_R[7] && cum_R[7]!=0.0){
        // turn a bound to an unbound D at + end
        fil[plus_end] = -1;
        // remove plus_end index from bound D (0) position index
        removeEntry_j(pos0, plus_end, D_bound_t-1, 8);
        // update concentrations of bound/unbound D
        D_bound_t = D_bound_t - 1;
        // update + end position
        plus_end = plus_end + 1;
        if(plus_end>maxN-1) plus_end = 0;
        //update filament length
        filamentLength = filamentLength - 1;
    }
    else if (r<=cum_R[8] && cum_R[8]!=0.0){
       // printf("9. Here filament at minus end: %.d\n", fil[minus_end]);
        // turn a bound to an unbound T at - end
        fil[minus_end] = -1;
        // remove minus_end index from bound T (1) position index
        removeEntry_j(pos1, minus_end, T_bound_t-1, 9);
        // update concentrations of bound/unbound T
        T_bound_t = T_bound_t - 1;
        // update - end position
        minus_end = minus_end - 1;
        if(minus_end<0) minus_end = maxN-1;
        //update filament length
        filamentLength = filamentLength - 1;
    }
    else {
        // turn a bound to an unbound T at + end
        fil[plus_end] = -1;
        // remove plus_end index from bound T (1) position index
        removeEntry_j(pos1, plus_end, T_bound_t-1, 10);
        // update concentrations of bound/unbound T
        T_bound_t = T_bound_t - 1;
        // update + end position
        plus_end = plus_end + 1;
        if(plus_end>maxN-1) plus_end = 0;
        //update filament length
        filamentLength = filamentLength - 1;
    }

}


















