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
    double t = 1000e-6;// simulation time [sec]
    double dt =1e-6; // t / Nt; // s
    double qdot_i=2167000;
    double ti = 0;
    

    // Define incident heat flux vector "qdot_mr"
    const int numRays_mr = 1;
    vector<double> qdot_mr;
    for (int r = 0; r < numRays_mr; r++) { qdot_mr.push_back(qdot_i); }
    //-------------------------------------------------------------------

    int numPtsPerRay_mr = 60;
    double T0_mr =297;
    int ns = 1;

    //---------------------Solid conduction----------------------------------------

        // Initialize MR vectors
    materialResponse mr;
    mr.Init_mr(numRays_mr, numPtsPerRay_mr, T0_mr);

     // Write results into data file       
    int w = 20; 
    
    ofstream mFractionsFile("mFrac.dat");
    mFractionsFile << left << setw(w) << "ns" << setw(w) << "t [s]" << setw(w) << "T [K]"  << setw(w) << "y_c10h160" << setw(w) <<  "y_n2"  << setw(w) << "y_o2" << setw(w) << "y_ar" <<endl;
    mFractionsFile << left << setw(w) <<ns-1<<setw(w)<< ti << setw(w) << T0_mr << setw(w) << 0 << setw(w) << 78./100 << setw(w) << 21./100 << setw(w) << 1./100 << endl;

    ofstream RatesFile("Rates.dat");
    RatesFile << left << setw(w) << "ns" << setw(w) <<"t [s]" << setw(w) << "T [K]" << setw(w) << "Xfront [mm]" << setw(w) << "sdot [mm/s]" << setw(w) << "omega_dot [kg/m^2s]" << setw(w) << "mdot_w [kg/m^2s]" << endl;
    RatesFile << left << setw(w) << ns << setw(w) << ti << setw(w) << T0_mr << setw(w) << mr.rays[0].x0[1]*1e3 << setw(w) << 0 << setw(w) << 0 << setw(w) << 0 << endl;
    // -----------------------------------------------------------------------------       


    //---------------------------
    //           TIME LOOP
    //-----------------------------
    while (ti <= t-dt) {        

        ti = ns * dt;       
        
         // Solve heat conduction equation and assess recession rate 
        mr.SolveCondRays(qdot_mr, numPtsPerRay_mr, dt, ns);
        ///////////////////////////////////////////////////////////////             

        // Write temperature data into file     
        mFractionsFile << left << setw(w) << ns << setw(w) << ti << setw(w) << (mr.rays[0].f0[0]+mr.rays[0].f0[1])/2 << setw(w) << mr.rays[0].yk_w_out[0] << setw(w) << mr.rays[0].yk_w_out[1] << setw(w) << mr.rays[0].yk_w_out[2] << setw(w) << mr.rays[0].yk_w_out[3] << endl;
        RatesFile << left << setw(w) <<ns << setw(w) << ti << setw(w) << (mr.rays[0].f0[0] + mr.rays[0].f0[1]) / 2 << setw(w) << mr.rays[0].x0[1]*1e3 << setw(w) << mr.rays[0].sdot_out*1e3 << setw(w) << mr.rays[0].mdot_k_out[0]  << setw(w) << mr.rays[0].mdot_w_out << endl;

        ns += 1;


    }

    //for (int j = 1; j < numPtsPerRay_mr + 1; j++) { TempFile << setprecision(20)<< mr.rays[0].f0[j] << " "; }
    //TempFile << endl;
      
    RatesFile.close();
    mFractionsFile.close();

    std::cout << "Hello World!" << endl;


    return 0;
}
