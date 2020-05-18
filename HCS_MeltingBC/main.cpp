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
    double t = 70;// simulation time [sec]
    double dt =0.005; // t / Nt; // s
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
    int w = 20; 
    
    ofstream mFractionsFile("mFrac21.dat");
    mFractionsFile << left << setw(w) << "t [s]" << setw(w) << "T [K]"  << setw(w) << "y_o" << setw(w) << "y_co" << setw(w) << "y_co2" << setw(w) << "y_c" << setw(w) << "y_c2" << setw(w) << "y_c3" << setw(w) << "y_o2" << setw(w) << "y_n2" << setw(w) << "y_n" << setw(w) << "y_no" << endl;
    mFractionsFile << left << setw(w) << ti << setw(w) << T0_mr << setw(w) << 4.34e-2/2 << setw(w) << 0 << setw(w) << 0 << setw(w) << 0 << setw(w) <<0 << setw(w) << 0 << setw(w) << 1.6e-1/2 << setw(w) << 7.5e-1/2 << setw(w) << 2.6e-3/2 << setw(w) << 4.4e-2/2 << endl;

    ofstream RatesFile("Rates21.dat");
    RatesFile << left << setw(w) <<"t [s]" << setw(w) << "T [K]" << setw(w) << "Xfront [mm]" << setw(w) << "sdot [mm/s]" << setw(w) << "dodt [kg/m^2s]" << setw(w) << "dcodt [kg/m^2s]" << setw(w) << "dco2dt [kg/m^2s]" << setw(w) <<"dcdt [kg/m^2s]"<< setw(w) <<"dc2dt [kg/m^2s]"<< setw(w) <<"dc3dt [kg/m^2s]"<< endl;
    RatesFile << left << setw(w) << ti << setw(w) << T0_mr << setw(w) << mr.rays[0].x0[1] * 1e3 << setw(w) << 0 << setw(w) << 0 << setw(w) << 0 << setw(w) << 0 << setw(w) << 0 << setw(w) << 0 << setw(w) << 0 << endl;
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
        RatesFile << left << setw(w) << ti << setw(w) << mr.rays[0].f0[1] << setw(w) << mr.rays[0].x0[1] * 1e3 << setw(w) << mr.rays[0].sdot_out * 1e3 << setw(w) << mr.rays[0].mdot_k_out[0] << setw(w) << mr.rays[0].mdot_k_out[1] << setw(w) << mr.rays[0].mdot_k_out[2] << setw(w) << mr.rays[0].mdot_k_out[3] << setw(w) << mr.rays[0].mdot_k_out[4] << setw(w) << mr.rays[0].mdot_k_out[5] << endl;
        mFractionsFile << left << setw(w) << ti << setw(w) << mr.rays[0].f0[1] << setw(w) << mr.rays[0].yk_w_out[0] << setw(w) << mr.rays[0].yk_w_out[1] << setw(w) << mr.rays[0].yk_w_out[2] << setw(w) << mr.rays[0].yk_w_out[3] << setw(w) << mr.rays[0].yk_w_out[4] << setw(w) << mr.rays[0].yk_w_out[5] << setw(w) << mr.rays[0].yk_w_out[6] << setw(w) << mr.rays[0].yk_w_out[7] << setw(w) << mr.rays[0].yk_w_out[8] << setw(w) << mr.rays[0].yk_w_out[9] << endl;

        i += 1;


    }

    //for (int j = 1; j < numPtsPerRay_mr + 1; j++) { TempFile << setprecision(20)<< mr.rays[0].f0[j] << " "; }
    //TempFile << endl;
      
    RatesFile.close();
    mFractionsFile.close();

    std::cout << "Hello World!" << endl;


    return 0;
}
