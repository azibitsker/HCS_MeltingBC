// Main with constant number of cells
/* 


#include "SolidCondRay.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

// This code solves transient heat conduction equation using implicit scheme ( dT/dt=alpha(d^2T/dx^2).
//The main function is HeatCondSolver with inputs: dt - time step, qdot- heat flux, f0 - temperature state at current time step. Output: Tn - Temperature at the next time step.

int main(){

    // ---------To be replaced by CFD simulation---------------------
    double t = 60; // simulation time [sec]
    //const int Nt = 2600; // number of time increments
    double dt = 0.01; // t / Nt; // s
    double qdot_i;
    double ti = 0;
    double k = 0.402, rho = 280, Cp = 983.9; // material properties
    double alpha = k / (rho * Cp); // thermal diffusivity

    qdot_i = 5e4; // W/m^2

    //---------------------Solid conduction----------------------------------------
    const int Nx = 2400; // number of material cells
    double T0 = 298; //initial material temperature [K]
    double sdot0;
    double f0[Nx + 2];  // Temperature solution array with Nx + 2 ghost cells
    for (int j = 0; j < Nx + 2; j++) { f0[j] = T0; }  // initialize f0

    double x0[Nx + 3]; // t=t0, location of the cells boundaries Nx+1 + 2 boundaries for the ghost cells
    double x[Nx + 3]; // t=t, array for the updated grid boundaries

    // Define
    SolidCond_Ray SolidTemp;

    //SolidTemp.GenerateGeomGrid(x0, Nx + 3); // Generate initial grid using geometric series
    //SolidTemp.GenerateGeomGrid(x, Nx + 3); // Generate initial grid using geometric series

    SolidTemp.GenerateUniformGrid(x0, Nx + 3); // Generate initial grid with uniform spacing
    SolidTemp.GenerateUniformGrid(x, Nx + 3); // Generate initial grid with uniform spacing

    // Write initial temperature vector into .dat file
    ofstream IterationsFile("N_iter.dat");
    IterationsFile << "i N_iter";
    IterationsFile << endl;

    ofstream NcellsFile("N_cells.dat");
    NcellsFile << "N_cells";
    NcellsFile << endl;

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
    double dx0;// = x[2] - x[1];
    //int fact = 20;
    //dt = fact *0.5 * pow(dx0, 2) / alpha;
    int i = 1;

    while (ti <= t) {
        //for(int i=1;i<Nt;i++){

        ti = i * dt;
        // Call function GetTemp to compute and update the solution vector f0
        ///////////////////////////////////////////////////////////////
        sdot0 = SolidTemp.HeatCondSolver(f0, x0, x, dt, qdot_i, Nx, i, IterationsFile, NcellsFile);
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

        //dx0 = x[2] - x[1];
        //dt = dx0 / sdot0;

        //if (sdot0 == 0 || dt>0.1) {            
            //dt = 0.1;
        //}


        //dx0 = x[2] - x[1];
        //dt = fact *0.5*pow(dx0,2) / alpha;

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
    double dt = 0.01; // t / Nt; // s
    double qdot_i=5e4;
    double ti=0;        

    //---------------------Solid conduction----------------------------------------
     // Define
    SolidCond_Ray SolidTemp;

    double dx = 1e-5;
    const int Nx = round(SolidTemp.L0/dx); // number of material cells
    double T0 = 298; //initial material temperature [K]
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
    IterationsFile << "i N_iter";
    IterationsFile << endl;

    ofstream NcellsFile("N_cells.dat");
    NcellsFile << "N_cells";
    NcellsFile << endl;

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
        sdot0 = SolidTemp.HeatCondSolver(f0, x0, x, dt, qdot_i, Nx, i, IterationsFile, NcellsFile);
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

        //dx0 = x[2] - x[1];
        //dt = dx0 / sdot0;

        //if (sdot0 == 0 || dt>0.1) {            
            //dt = 0.1;
        //}
        
        
        //dx0 = x[2] - x[1];
        //dt = fact *0.5*pow(dx0,2) / alpha;

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




