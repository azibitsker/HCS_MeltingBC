// Main with constant number of cells 

#include "materialResponse.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

//#define pi 3.14159265359


// This code solves transient heat conduction equation using implicit scheme ( dT/dt=alpha(d^2T/dx^2).
//The main function is HeatCondSolver with inputs: dt - time step, qdot- heat flux, f0 - temperature state at current time step. Output: Tn - Temperature at the next time step.

int main(){   

    // ---------To be replaced by CFD simulation---------------------
    double t = 0.077;// simulation time [sec]
    double dt =1e-4; // t / Nt; // s
    double qdot_i=2e6;
    double ti = 0;
    int ntSkip = 1;
    

    // Define incident heat flux vector "qdot_mr"
    const int numRays_mr = 1;
    vector<double> qdot_mr;
    for (int r = 0; r < numRays_mr; r++) { qdot_mr.push_back(qdot_i); }
    //-------------------------------------------------------------------

    int numPtsPerRay_mr = 100;
    double T0_mr = 300;
    int i = 1;
    int p = 6;
    //---------------------Solid conduction----------------------------------------

        // Initialize MR vectors
    materialResponse mr;
    mr.Init_mr(numRays_mr, numPtsPerRay_mr, T0_mr);

    string m ="0";
     // Write initial temperature vector into .dat file        
    
    ofstream BoundaryFile("x" + m +".dat");
    BoundaryFile << mr.rays[0].x0[1] ;
    for (int j = 1; j < mr.rays[0].x0.size() - 2; j++) { BoundaryFile << " " << (mr.rays[0].x0[j] + mr.rays[0].x0[j + 1]) / 2 ; }
    BoundaryFile << "  " << mr.rays[0].x0[numPtsPerRay_mr + 1];
    BoundaryFile << endl;

    ofstream RecessionRateFile("sdot" + m + ".dat");
    RecessionRateFile << 0 << endl;
    ofstream TimeFile("time" + m + ".dat");
    TimeFile << ti << endl;
    ofstream TempFile("Temp" + m + ".dat");
    TempFile << setprecision(p) << (mr.rays[0].f0[0] + mr.rays[0].f0[1]) / 2; // front surface temperature
    for (int j = 1; j < mr.rays[0].f0.size() - 1; j++) { TempFile << " " << setprecision(p) << mr.rays[0].f0[j]; } // cell centered values
    TempFile << " " << setprecision(p) << (mr.rays[0].f0[numPtsPerRay_mr] + mr.rays[0].f0[numPtsPerRay_mr + 1]) / 2;
    TempFile << endl;
    // -----------------------------------------------------------------------------      


    //---------------------------
    //           TIME LOOP
    //-----------------------------
    while (ti <= t-dt) {        

        ti = i * dt;       
        
         // Solve heat conduction equation and assess recession rate 
        mr.SolveCondRays(qdot_mr, numPtsPerRay_mr, dt);
        ///////////////////////////////////////////////////////////////          

        if (i % ntSkip == 0)
        {       

            // Write temperature data into file    
            BoundaryFile << mr.rays[0].x0[1];
            for (int j = 1; j < mr.rays[0].x0.size()-2; j++) { BoundaryFile <<" "<< (mr.rays[0].x0[j]+ mr.rays[0].x0[j+1])/2 ; }
            BoundaryFile <<"  "<< mr.rays[0].x0[numPtsPerRay_mr+1] ;
            BoundaryFile << endl;

            // Write the recession rate into file
            RecessionRateFile << mr.rays[0].sdot_out<< endl;
            TimeFile << ti << endl;

            TempFile << setprecision(p) << (mr.rays[0].f0[0]+ mr.rays[0].f0[1] )/2; // front surface temperature
            for (int j = 1; j < mr.rays[0].f0.size()-1; j++) { TempFile <<" " << setprecision(p) << mr.rays[0].f0[j]; } // cell centered values
            TempFile << " " << setprecision(p) << (mr.rays[0].f0[numPtsPerRay_mr]+ mr.rays[0].f0[numPtsPerRay_mr+1])/2;
            TempFile << endl;
        }

        i++;

    }

    //for (int j = 1; j < numPtsPerRay_mr + 1; j++) { TempFile << setprecision(20)<< mr.rays[0].f0[j] << " "; }
    //TempFile << endl;

    TempFile.close();
    BoundaryFile.close();
    RecessionRateFile.close();
    TimeFile.close();    


    cout << "Hello World!" << endl;


    return 0;
}
