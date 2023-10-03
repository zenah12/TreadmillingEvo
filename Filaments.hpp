class Filaments
{
private:
	int x;
public:
    
    double Ri[numR];                                // reaction rates
    int filament[maxN];                             // filament
    int filamentStates[numt/saveEnd][maxN];                    // different fillaments at M time points following steady state to explroe filament composison
    int pos_0[maxN];                                // position of bound 0s (D blocks)
    int pos_1[maxN];                                // position of bound 1s (T blocks)
    
    //double FilamentLengthOverTime[numt][2];         // stores the length of the filament and the time
    
    void update_Ri(int fil[maxN], double Ri[numR], int length);
    double sample_dt(double Ri[numR]);
    void update_f(int fil[maxN], double Ri[numR], int pos0[maxN], int pos1[maxN]);
};


