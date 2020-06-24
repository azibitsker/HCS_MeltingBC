// Main with constant number of cells 

#include "materialResponse.h"
#include "solidCondRay.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>


// This code solves transient heat conduction equation using implicit scheme ( dT/dt=alpha(d^2T/dx^2).
//The main function is HeatCondSolver with inputs: dt - time step, qdot- heat flux, f0 - temperature state at current time step. Output: Tn - Temperature at the next time step.

int main(){   

    // ---------To be replaced by CFD simulation---------------------
    double t = 30;// simulation time [sec]
    double dt =0.001; // t / Nt; // s
    double qdot_i=0;
    double ti = 0;
    

    // Define incident heat flux vector "qdot_mr"
    const int numRays_mr = 1;
    vector<double> qdot_mr;
    for (int r = 0; r < numRays_mr; r++) { qdot_mr.push_back(qdot_i); }
    //-------------------------------------------------------------------

    int numPtsPerRay_mr = 20;
    double T0_mr = 1000;
    int i = 1;
    int N_iter = 50;


    //---------------------Solid conduction----------------------------------------

        // Initialize MR vectors
    materialResponse mr;    
    mr.Init_mr(numRays_mr, numPtsPerRay_mr, T0_mr);

     // Write initial temperature vector into .dat file        
    
    ofstream BoundaryFile("x.dat");
    for (int j = 1; j < mr.rays[0].x0.size() - 1; j++) { BoundaryFile << mr.rays[0].x0[j] * 1e3 << " "; }
    BoundaryFile << endl;
    
    ofstream TempFile("temp.dat");
    for (int j = 1; j < mr.rays[0].f0.size() - 1; j++) { TempFile << mr.rays[0].f0[j] << " "; }
    TempFile << endl;
    // -----------------------------------------------------------------------------        
        
         // Solve heat conduction equation on original mesh
        mr.SolveCondRays(qdot_mr, numPtsPerRay_mr, dt);        
        // Write the temperature into file
        for (int j = 1; j < mr.rays[0].f.size()-1; j++) { TempFile << setprecision(20) << mr.rays[0].f[j]<<" "; }
        TempFile << endl;  
        // Write grid into file
        for (int j = 1; j < mr.rays[0].x.size() - 1; j++) { BoundaryFile << mr.rays[0].x[j] * 1e3 << " "; }
        BoundaryFile << endl;

        // Move mesh
        mr.rays[0].moveInternalMesh(N_iter); // moves mesh of vector x
        // Solve heat conduction equation on moved mesh
        mr.SolveCondRays(qdot_mr, numPtsPerRay_mr, dt);
        // Write the temperature into file
        for (int j = 1; j < mr.rays[0].f.size() - 1; j++) { TempFile << setprecision(20) << mr.rays[0].f[j] << " "; }
        TempFile << endl;
       // Write grid into file
       for (int j = 1; j < mr.rays[0].x.size() - 1; j++) { BoundaryFile << setprecision(20) << mr.rays[0].x[j] * 1e3 << " "; }
       BoundaryFile << endl;      
        

    TempFile.close();
    BoundaryFile.close();    
   
    cout << "Hello World!" << endl;


    return 0;
}
