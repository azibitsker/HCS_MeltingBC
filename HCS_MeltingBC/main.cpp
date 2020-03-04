#include "SolidCondRay.h"


// This code solves transient heat conduction equation using implicit scheme ( dT/dt=alpha(d^2T/dx^2).
//The main function is HeatCondSolver with inputs: dt - time step, qdot- heat flux, f0 - temperature state at current time step. Output: Tn - Temperature at the next time step.

int main()
{

    // ---------To be replaced by CFD simulation---------------------
    double t = 30; // simulation time [sec]
    const int Nt = 1200; // number of time increments
    double dt, qdot_i;

    dt = t / Nt; // s
    qdot_i = 5e4; // W/m^2

    //---------------------Solid conduction----------------------------------------
    const int Nx = 1600; // number of material cells
    double T0 = 298; //initial material temperature [K]

    double f0[Nx + 2];  // Temperature solution array with Nx + 2 ghost cells
    for (int j = 0; j < Nx + 2; j++) { f0[j] = T0; }  // initialize f0

    double x0[Nx + 3]; // t=t0, location of the cells boundaries Nx+1 + 2 boundaries for the ghost cells
    double x[Nx + 3]; // t=t, array for the updated grid boundaries

    // Define
    SolidCond_Ray SolidTemp;

    SolidTemp.GenerateGrid(x0, Nx + 3); // Generate initial grid using geometric series
    SolidTemp.GenerateGrid(x, Nx + 3); // Generate initial grid using geometric series

    // Write initial temperature vector into .dat file

    ofstream TempFile("Temp_C.dat");
    for (int j = 1; j < 2; j++) { TempFile << f0[j] << " "; }
    TempFile << endl;

    ofstream BoundaryFile("X_C.dat");
    for (int j = 1; j < 2; j++) { BoundaryFile << (x0[j]) * 1e3 << " "; }
    BoundaryFile << endl;

    ofstream RecessionRateFile("Sdot.dat");
    RecessionRateFile << 0;
    RecessionRateFile << endl;

    //----- CFD Simulation Loop------

    for (int i = 1; i < Nt + 1; i++) {

        double sdot0 = 0; // [m/s] the initial guess for recession rate
        // Call function GetTemp to compute and update the solution vector f0
        sdot0 = SolidTemp.HeatCondSolver(f0, x0, x, dt, qdot_i, sdot0, Nx, i);

        // Write temperature data into file
        for (int j = 1; j < 2; j++) { TempFile << f0[j] << " "; }
        TempFile << endl;
        // Write grid boundary locations into file
        for (int j = 1; j < 2; j++) { BoundaryFile << (x0[j]) * 1e3 << " "; }
        BoundaryFile << endl;

        // Write the recession rate into file
        RecessionRateFile << sdot0 * 1e3 << endl;

    }
    TempFile.close();
    BoundaryFile.close();
    RecessionRateFile.close();


    //  Create simulation time .dat file
    ofstream TimeFile("Time_C.dat");
    if (TimeFile.is_open())
    {
        for (int j = 0; j < Nt + 1; j++) {
            TimeFile << j * dt;
            TimeFile << endl;
        }
        TimeFile.close();
    }
    else cout << "Unable to open file";


    cout << "Hello World - Hello!" << endl;
    cout << "Hello World - Hello!" << endl;
    cout << "Hello World - Hello!" << endl;

    return 0;
}




