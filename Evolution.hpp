class Evolution
{
private:
	int x;
public:
    
    double paramsEvo[numt][paramsN];                        // evolution of parameters
    double meanSE[2][2];                                    // stores the mean and SE at time t and t-1
    
    double debugArray[numt][6];
    
    void selection(double meansSEs[2][2], double parsEvo[numt][paramsN], double ratio);
    void mutation();
    void updateMeanSE(double mean_now, double mean2_now, double meanSE[2][2], int index);
    void updateParamsEvo(double parsEvo[numt][paramsN], double meanNow, double seNow);
    void storeCurrentFilament(int filamentStates[M][maxN], int filament[maxN], int count);
    double getMaxRatio(double mean1, double mean2, double mean3, double mean4);
    
};


