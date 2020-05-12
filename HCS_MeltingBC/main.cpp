// Main with constant number of cells 

#include "materialResponse.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>


// This code solves transient heat conduction equation using implicit scheme ( dT/dt=alpha(d^2T/dx^2).
//The main function is HeatCondSolver with inputs: dt - time step, qdot- heat flux, f0 - temperature state at current time step. Output: Tn - Temperature at the next time step.

int main(){   

    // ---------To be replaced by CFD simulation---------------------
    double t = 50;// simulation time [sec]
    double dt =0.01; // t / Nt; // s
    double qdot_i=6e6;
    double ti = 0;
    

    // Define incident heat flux vector "qdot_mr"
    const int numRays_mr = 1;
    vector<double> qdot_mr;
    for (int r = 0; r < numRays_mr; r++) { qdot_mr.push_back(qdot_i); }
    //-------------------------------------------------------------------

    int numPtsPerRay_mr = 20;
    double T0_mr = 800;
    int i = 1;

    //---------------------Solid conduction----------------------------------------

        // Initialize MR vectors
    materialResponse mr;
    mr.Init_mr(numRays_mr, numPtsPerRay_mr, T0_mr);

     // Write results into data file       
    int w = 15;

    ofstream ResultsFile("Results.dat");
    ResultsFile << left<< setw(w)<< "t" << setw(w) << "T" << setw(w) << "sdot_tot"<< setw(w) << "Xfront" << endl;
    ResultsFile << left << setw(w) << ti << setw(w) << T0_mr << setw(w) << 0 << setw(w) << mr.rays[0].x0[1]* 1e3 << endl;

    ofstream RatesFile("Rates.dat");
    RatesFile << left << setw(w) <<"t"<< setw(w) << "sdot_tot" << setw(w) << "dodt" << setw(w) << "dcodt" << setw(w) << "dco2dt" << setw(w) <<"dcdt"<< setw(w) <<"dc2dt"<< setw(w) <<"dc3dt"<< endl;
    RatesFile << left << setw(w) << ti << setw(w) << 0 << setw(w) << 0 << setw(w) << 0 << setw(w) << 0 << setw(w) << 0 << setw(w) << 0 << setw(w) << 0 << endl;
    
    // -----------------------------------------------------------------------------
       


    //---------------------------
    //           TIME LOOP
    //-----------------------------
    while (ti <= t-dt) {        

        ti = i * dt;       
        
         // Solve heat conduction equation and assess recession rate 
        mr.SolveCondRays(qdot_mr, numPtsPerRay_mr, dt);
        ///////////////////////////////////////////////////////////////             

        // Write temperature data into file
        ResultsFile << left << setw(15) << ti << " " << setw(15) << mr.rays[0].f0[1] << " " << setw(15) << mr.rays[0].sdot_out * 1e3<< setw(15) << mr.rays[0].x0[1] * 1e3 << endl;
        RatesFile << left << setw(w) << ti << setw(15) << mr.rays[0].sdot_out * 1e3 << setw(w) << mr.rays[0].mdot_sp[0] << setw(w) << mr.rays[0].mdot_sp[1] << setw(w) << mr.rays[0].mdot_sp[2] << setw(w) << mr.rays[0].mdot_sp[3] << setw(w) << mr.rays[0].mdot_sp[4] << setw(w) << mr.rays[0].mdot_sp[5] << endl;

        i += 1;


    }

    //for (int j = 1; j < numPtsPerRay_mr + 1; j++) { TempFile << setprecision(20)<< mr.rays[0].f0[j] << " "; }
    //TempFile << endl;
      
    ResultsFile.close();

    std::cout << "Hello World!" << endl;


    return 0;
}
