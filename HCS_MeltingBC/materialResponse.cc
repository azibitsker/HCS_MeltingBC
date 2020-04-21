#include "materialResponse.h"


void materialResponse::Init_mr(int numRays, int numPtsPerRay_in, double T0_in)
{
    rays.resize(numRays);

    //Initialize temperature and grid state vectors
    for (int r = 0; r < numRays; r++)
    {
        rays[r].init(T0_in, numPtsPerRay_in);
    }

}

void materialResponse::SolveCondRays(vector<double> qdot_in, int numPtsPerRay_in, double dt_in, ofstream & IterationsFile)
{

    int Nx = numPtsPerRay_in;
    size_t numRays = qdot_in.size();
    double timeInc = dt_in;

    // ofstream RecessionRateFile("Sdot.dat");
     //RecessionRateFile << 0 << endl;

    for (int r = 0; r < numRays; ++r)
    {
        rays[r].HeatCondSolver(qdot_in[r], Nx, timeInc, IterationsFile);

    }

}