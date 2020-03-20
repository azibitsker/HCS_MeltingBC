// Main with constant number of cells 

#include "materialResponse.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>


// This code solves transient heat conduction equation using implicit scheme ( dT/dt=alpha(d^2T/dx^2).
//The main function is HeatCondSolver with inputs: dt - time step, qdot- heat flux, f0 - temperature state at current time step. Output: Tn - Temperature at the next time step.

int main(){   

    // ---------To be replaced by CFD simulation---------------------
    double t = 60; // simulation time [sec]
    double dt = 0.01; // t / Nt; // s
    double qdot_i;
    double ti = 0;
    qdot_i = 5e4; // W/m^2

    // Define incident heat flux vector "qdot_mr"
    const int numRays_mr = 1;
    vector<double> qdot_mr;
    for (int r = 0; r < numRays_mr; r++) { qdot_mr.push_back(qdot_i); }
    //-------------------------------------------------------------------

    int numPtsPerRay_mr = 200;
    double T0_mr = 298;
    int i = 1;

    //---------------------Solid conduction----------------------------------------

     // Write initial temperature vector into .dat file
    ofstream IterationsFile("N_iter.dat");
    IterationsFile << "i N_iter";
    IterationsFile << endl;
    ofstream NcellsFile("N_cells.dat");
    NcellsFile << "N_cells";
    NcellsFile << endl;
    ofstream BoundaryFile("X_C.dat");
    //for (int j = 1; j < 2; j++) { BoundaryFile << (x0[j]) * 1e3; }
    BoundaryFile << endl;
    ofstream RecessionRateFile("Sdot.dat");
    RecessionRateFile << 0 << endl;
    ofstream TimeFile("Time_C.dat");
    TimeFile << ti << endl;
    ofstream TempFile("Temp_C.dat");
    for (int j = 1; j < 2; j++) { TempFile << T0_mr; }
    TempFile << endl;
    // -----------------------------------------------------------------------------
       
    // Initialize MR vectors
    materialResponse mr;
    mr.Init_mr(numRays_mr, numPtsPerRay_mr, T0_mr);

    //---------------------------
    //           TIME LOOP
    //-----------------------------
    while (ti <= t) {        

        ti = i * dt;       
        
         // Solve heat conduction equation and assess recession rate 
        mr.SolveCondRays(qdot_mr, numPtsPerRay_mr, dt);
        ///////////////////////////////////////////////////////////////             

        // Write temperature data into file
        for (int j = 1; j < 2; j++) { TempFile << mr.rays[0].f0[j]; }
        TempFile << endl;
        // Write grid boundary locations into file
        for (int j = 1; j < 2; j++) { BoundaryFile << (mr.rays[0].x0[j]) * 1e3; }
        BoundaryFile << endl;
        // Write the recession rate into file
        RecessionRateFile << mr.rays[0].sdot_out* 1e3 << endl;
        TimeFile << ti << endl;
   
        i += 1;


    }
    TempFile.close();
    BoundaryFile.close();
    RecessionRateFile.close();
    TimeFile.close();
    IterationsFile.close();
    NcellsFile.close();


    cout << "Hello World!" << endl;


    return 0;
}



/* 

// main with constant increment size

#include "SolidCondRay.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

// This code solves transient heat conduction equation using implicit scheme ( dT/dt=alpha(d^2T/dx^2).
//The main function is HeatCondSolver with inputs: dt - time step, qdot- heat flux, f0 - temperature state at current time step. Output: Tn - Temperature at the next time step.

int main()
{

    // ---------To be replaced by CFD simulation---------------------
    double t = 30; // simulation time [sec]
    double dt = 1e-2; // t / Nt; // s
    double qdot_i=5e4;
    double ti=0;        

    //---------------------Solid conduction----------------------------------------
     // Define
    SolidCond_Ray SolidTemp;

    double dx = 1e-5;
    const int Nx = round(SolidTemp.L0/dx); // number of material cells
    double T0 = 300; //initial material temperature [K]
    double sdot0;
    vector<double> f0(Nx + 2);  // Temperature solution array with Nx + 2 ghost cells
    for (int j = 0; j < Nx + 2; j++) { f0[j] = T0; }  // initialize f0

    vector<double> x0; // t=t0, location of the cells boundaries Nx+1 + 2 boundaries for the ghost cells
    vector<double> x; // t=t, array for the updated grid boundaries
    //SolidTemp.GenerateGeomGrid(x0, Nx + 3); // Generate initial grid using geometric series
    //SolidTemp.GenerateGeomGrid(x, Nx + 3); // Generate initial grid using geometric series
    SolidTemp.GenerateUniformGrid(x0, dx); // Generate initial grid with uniform spacing
    SolidTemp.GenerateUniformGrid(x, dx); // Generate initial grid with uniform spacing

    // Write initial temperature vector into .dat file
     ofstream IterationsFile("N_iter.dat");
    IterationsFile << "N_iter";
    IterationsFile << endl;

    ofstream NcellsFile("N_cells.dat");
    NcellsFile << "N_cells";
    NcellsFile << endl;

    ofstream sdotIter("sdotIter.dat");
    
    ofstream TempFile("Temp_C.dat");
    for (int j = 1; j < 2; j++) { TempFile << f0[j]; }
    TempFile << endl;

    ofstream BoundaryFile("X_C.dat");
    for (int j = 1; j < 2; j++) { BoundaryFile << (x0[j]) * 1e3; }
    BoundaryFile << endl;

    ofstream RecessionRateFile("Sdot.dat");
    RecessionRateFile << 0 << endl;
    
    ofstream TimeFile("Time_C.dat");
    TimeFile << ti << endl;

    //----- CFD Simulation Loop------
    int i = 1;    

    while(ti<=t) {
    //for(int i=1;i<Nt;i++){

        ti = i * dt;
        // Call function GetTemp to compute and update the solution vector f0
        ///////////////////////////////////////////////////////////////
        sdot0 = SolidTemp.HeatCondSolver(f0, x0, x, dt, qdot_i, Nx, i, IterationsFile, NcellsFile,sdotIter);
        ///////////////////////////////////////////////////////////////             

        // Write temperature data into file
        for (int j = 1; j < 2; j++) { TempFile << f0[j]; }
        TempFile << endl;
        // Write grid boundary locations into file
        for (int j = 1; j < 2; j++) { BoundaryFile << (x0[j]) * 1e3; }
        BoundaryFile << endl;
        // Write the recession rate into file
        RecessionRateFile << sdot0 * 1e3 << endl;
        TimeFile << ti << endl;    

    
        i += 1;
        

    }
    TempFile.close();
    BoundaryFile.close();
    RecessionRateFile.close();
    TimeFile.close();
    IterationsFile.close();
    NcellsFile.close();


    cout << "Hello World!" << endl;


    return 0;
}

*/