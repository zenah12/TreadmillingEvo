#ifndef GLOBAL
#define GLOBAL //Global.hpp


//Define global fixed parameters
#define numt (5000000)//(300000)                       //maximum time steps for experiment
#define maxN (2000)                        // maximum filement size
#define M (500)                             // number of filaments to store at the end of evolution#define evoSteps (100)                      // maximum number of evolution steps
#define saveEnd (10000)                       // every how many steps to save filament state for struscture assessment at the end
#define numR (10)                           // number of posible reactions
#define numSave (10)                        // period for saving data
#define paramsN (14)                        // total number of system parameters
#define Pi (3.141592653589793238462643383279502884197169399375)


//define global parameters that can be modified in main
// parameters//
// rates
extern double kon_D_p;
extern double koff_D_p;
extern double kon_T_p;
extern double koff_T_p;
extern double kon_D_m;
extern double koff_D_m;
extern double kon_T_m;
extern double koff_T_m;
extern double w_de;
extern double w_re;
//
extern double kon_D_p0;
extern double koff_D_p0;
extern double kon_T_p0;
extern double koff_T_p0;
extern double kon_D_m0;
extern double koff_D_m0;
extern double kon_T_m0;
extern double koff_T_m0;
extern double w_de0;
extern double w_re0;
// evolution
extern double mu;               // mutation rate
extern double w_re_max;
extern double w_re_min;
extern double targetMean;
extern double targetRatio;
extern int evoTimeNow;          // current selection step
extern double tMin;
extern double meanNow;
extern double sdNow;
extern double mutationStep;
extern int maxEvo;
extern int tEvo;
// mutation indeces, the parameter mutates in simulation if = 1 and does not mutate if 0.
extern int mu_konTp;
extern int mu_konTm;
extern int mu_koffTp;
extern int mu_konDp;
extern int mu_konDm;
extern int mu_koffDp;
extern int mu_w_de;
extern int mu_w_re;
extern int mu_konkoffT;
extern int mu_konkoffD;

extern double semIndex;

// initial conditions parameters
extern int D_bound_0;
extern int T_bound_0;
// variables
extern int D_bound_t;
extern int T_bound_t;

// max and minimum rates
extern double kmin;
extern double kmax;

// indices for minus and plus ends for filament
extern int minus_end;
extern int plus_end;
extern int filamentLength;

extern int indOut;

extern int tEvo;

extern double rMean;

extern int samplingSteps; //number of Gillepsie algorithm steps for sampling intervals

extern double win;        // window far from target mean allowed


// definitions //
// D (NDP bound) state is 0, T (NTP bound) state is 1, empty entry is at -1
// + end is on the "left" (beginning of matrix)
// - end is on the "right" (end of matrix)

#endif
