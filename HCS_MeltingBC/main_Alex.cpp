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
    double t = 40;// simulation time [sec]
    double dt =0.125; // t / Nt; // s
    double qdot_i;
    double ti = 0;
    qdot_i = 7.5e5; // W/m^2

    // Define incident heat flux vector "qdot_mr"
    const int numRays_mr = 1;
    vector<double> qdot_mr;
    for (int r = 0; r < numRays_mr; r++) { qdot_mr.push_back(qdot_i); }
    //-------------------------------------------------------------------

    int numPtsPerRay_mr = 20;
    double T0_mr = 300;
    int i = 1;

    //---------------------Solid conduction----------------------------------------

        // Initialize MR vectors
    materialResponse mr;
    mr.Init_mr(numRays_mr, numPtsPerRay_mr, T0_mr);

     // Write initial temperature vector into .dat file    
    ofstream NcellsFile("N_cells.dat");
    NcellsFile << "N_cells";
    NcellsFile << endl;
    ofstream IterationsFile("N_iter.dat");
    ofstream BoundaryFile("X_C.dat");
    for (int j = 1; j < 2; j++) { BoundaryFile << mr.rays[0].x0[j]* 1e3<<" "; }
    BoundaryFile << endl;
    ofstream RecessionRateFile("Sdot.dat");
    RecessionRateFile << 0 << endl;
    ofstream TimeFile("Time_C.dat");
    TimeFile << ti << endl;
    ofstream TempFile("Temp_C.dat");
    //for (int j = 1; j <2; j++) { TempFile << T0_mr << " "; }
    //TempFile << endl;
    // -----------------------------------------------------------------------------
       


    //---------------------------
    //           TIME LOOP
    //-----------------------------
    while (ti <= t-dt) {        

        ti = i * dt;       
        
         // Solve heat conduction equation and assess recession rate 
        mr.SolveCondRays(qdot_mr, numPtsPerRay_mr, dt, IterationsFile);
        ///////////////////////////////////////////////////////////////             

        // Write temperature data into file
        //for (int j = 1; j < numPtsPerRay_mr+1; j++) { TempFile << mr.rays[0].f0[j] << " "; }
        //TempFile << endl;
        // Write grid boundary locations into file
        for (int j = 1; j < 2; j++) { BoundaryFile << mr.rays[0].x0[j] * 1e3 << " "; }
        BoundaryFile << endl;
        // Write the recession rate into file
        RecessionRateFile << mr.rays[0].sdot_out* 1e3 << endl;
        TimeFile << ti << endl;
   
        i += 1;


    }

    for (int j = 1; j < numPtsPerRay_mr + 1; j++) { TempFile << setprecision(20)<< mr.rays[0].f0[j] << " "; }
    //TempFile << endl;

    TempFile.close();
    BoundaryFile.close();
    RecessionRateFile.close();
    TimeFile.close();
    IterationsFile.close();
    NcellsFile.close();


    cout << "Hello World!" << endl;


    return 0;
}
