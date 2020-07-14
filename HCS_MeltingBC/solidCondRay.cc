#include "solidCondRay.h"
#include "materialResponse.h"

#include "chemkin.h"
#include "momentumBalance.h"
#include "massBalance.h"


void solidCond_Ray::HeatCondSolver(double qdot_in, int Nx, double dt) {        
    
    double sdot; // computed recession rate      
    double mdot_c; // carbon flux from the surface
    double Tw; // wall temperature
    

    vector<double> f(Nx + 2); // current temperature state vector
    for (int j = 0; j < Nx + 2; j++) { f[j] = f0[j]; }  
        
    double p_eta = 1e-2*101350.; // stagnation pressure

    // convergence parameters
    double norm_yk0 = 0;
    double conv_sdot;
    double conv_yk_w = 1;   
    double norm_yk;

    double conv_Tw;
    double Tw0=0;

    t_check = dt * i_check;
    i_check++;

    chemkin SR;
    massBalance MasB;
    momentumBalance MomB;

    // initial guesses
    double sdot0 = 0;
    double yk0_w =  1. / SR.Ns; // mass fractions at the wall
    double ps_c = 0; // partial pressure of camphor species for sublimation   

    // initialize 
    SR.init_chemParam();        
    MasB.init_MassBal(SR, yk0_w); // initialize initial species concentrations at the wall with a guess
    mdot_w_out = 0;

    while (conv_yk_w > 1e-7) {
        conv_sdot = 1;
        norm_yk = 0;
        // convergence on sdot for the computed surface oxidation and sublimation
        while (conv_sdot > 1e-7) {

            // Evaluate temperature at the next time step, vector f is updated
            EvaluateTemp(f, qdot_in, sdot0, dt, Nx);
            Tw = f[1]; //wall temperature            
            SR.getSublimRates(Tw, ps_c); // update values in vector mdot_k
            SR.getSpeciesFlux(); // evaluates mdot_w    
            mdot_w_out = SR.mdot_w;

            // Compute Camphor mass flux from the surface
            mdot_c = SR.mdot_k[0]; //[kg/m^2sec] 
            sdot = mdot_c / rho_s;

            ContractGridBoundaries(sdot, dt, Nx + 3); // array x is updated

            conv_sdot = abs(sdot - sdot0) / sdot;
            conv_Tw = abs(Tw - Tw0) / Tw;

            Tw0 = Tw;
            sdot0 = sdot; // update the guess value (sdot + sdot0) / 2;            

        }        

        //cout << SR.mdot_sub[0] << " " << SR.mdot_sub[1] << " " << SR.mdot_sub[2] << " " << endl;

        
        MomB.solveMomentumBalance(SR, MasB, p_eta, Tw); // compute velocity,density and pressure at the wall
        MasB.solveMassBal(SR, MomB.rho_w, x); // solve for mass fractions in the ghost cell yk_ghost and yk_w

        // update partial pressures of camphor      

            ps_c = MasB.yk_w[0] * MomB.rho_w* SR.R0 / SR.Mw_s[0] * Tw;
        

        // compute norm of yk_w to check the convergence
        for (int i = 0; i < SR.Ns; i++) {
            norm_yk += MasB.yk_w[i] * MasB.yk_w[i];
        }

        conv_yk_w = abs(norm_yk - norm_yk0) / norm_yk;
        norm_yk0 = norm_yk;

    }

       //cout << "conv_Tw = " << conv_Tw << endl;

       // Update  x0
       for (int j = 0; j < Nx + 3; j++) { x0[j] = x[j]; }
       sdot_out = sdot0;

       // Update  f0
       for (int j = 0; j < Nx + 2; j++) { f0[j] = f[j]; }

       // Update  species mass flux
       mdot_k_out.resize(SR.Ns);
       yk_w_out.resize(SR.Ns);

       mdot_w_out = SR.mdot_w;
       for (int j = 0; j < SR.Ns; j++) { mdot_k_out[j] = SR.mdot_k[j]; }
       for (int j = 0; j < SR.Ns; j++) { yk_w_out[j] = MasB.yk_w[j]; }

   
}

// Generate initial material grid including nodes for the two ghost cells using geometric series
void solidCond_Ray::GenerateGeomGrid(const int Nx) {
    // Generate the initial material grid. Fill out the array x0 with location of the
    // cells boundaries using geometric series

     
    double dx0;
    double dxj;    

    dx0 = L0 * (1 - F) / (1 - pow(F, Nx - 3)); // initial grid cell size

    x0.push_back(-dx0);
    x0.push_back(0); // 1st ghost cell - its size is equal to the size of the first material cell

    for (int j = 2; j < Nx - 1; j++) {

        dxj = pow(F, j - 2) * dx0;
        x0.push_back(x0[j - 1] + dxj); // space the cell boundaries using geometric series
    }

    x0.push_back(x0[Nx - 2] + dxj); // last ghost cell - its size is equal to the size of the last material cell
    

    /* 10 geometric series grid with first cell divided into additional 10 geometric series grid
    double dx01, dx02,dxj;

    dx02 = L0 * (1 - F) / (1 - pow(F,10)); // initial grid cell size
    dx01= dx02 * (1 - F) / (1 - pow(F, 10)); // initial grid cell size

    x0.push_back(-dx01);
    x0.push_back(0); // 1st ghost cell - its size is equal to the size of the first material cell

    int j;
    for (j = 2; j <12; j++) {

        dxj = pow(F, j - 2) * dx01;
        x0.push_back(x0[j - 1] + dxj); // space the cell boundaries using geometric series
    }

    int jj;
    for (jj = j; jj < 21; jj++) {

        dxj = pow(F, jj) * dx02;
        x0.push_back(x0[jj - 1] + dxj); // space the cell boundaries using geometric series
    }

    x0.push_back(x0[jj-1] + dxj); // last ghost cell - its size is equal to the size of the last material cell

    */

    for (int j = 0; j < Nx; j++) { x.push_back(x0[j]); }
    //cout<< "L0 = "<< x0[Nx-2]-x0[1]<< endl;

}

void solidCond_Ray::GenerateUniformGrid(int Nx) {
    // Generate the initial material grid. Fill out the array x0 with location of the
    // cells boundaries using uniform series
    
    
    double dx0;    
    dx0 = L0/ (Nx - 3); // initial grid cell size
        
    x0.push_back(-dx0);

    for (int j = 1; j < Nx; j++) {

        x0.push_back(x0[j - 1] + dx0); // space the cell boundaries using uniform spacing
    }

    //x0.push_back(x0[Nx - 2] + dx0); // last ghost cell - its size is equal to the size of the last material cell    
       
    for (int j = 0; j < Nx; j++) { x.push_back(x0[j]); }

}

// Update grid boundaries based on the computed recession rate sdot0
void solidCond_Ray::ContractGridBoundaries(double sdot0, double dt, const int Nx) {

    double dx_melt = sdot0 * dt; // total length of the melted cells
    double L;

    L = x0[Nx - 2] - x0[1]; // current material length

//------- Update location of each cell boundary relative to the initial location at t=0------------------
    double dx0, dxj;

    x[1] = x0[1] + dx_melt; // update location of the recession front
    double Recession = x[1];

    if (Recession >= L0) {
        cout << "The whole material has melted away" << endl;
        return;
    }

    for (int j = 2; j < Nx - 1; j++) {

        dx0 = x0[j] - x0[j - 1]; // compute material initial cell length
        dxj = dx0 - dx0 / L * dx_melt; // compute the new cell length

        x[j] = x[j - 1] + dxj; // update location of each cell boundary

    }
    x[0] = x[1] - (x[2] - x[1]); // location of the 1st ghost cell boundary
    x[Nx - 1] = x[Nx - 2] + dxj;  // location of the last ghost cell boundary

}

// Get a,b,c coefficients for Thomas algorithm
void solidCond_Ray::Get_abcd_coeff( vector<double>& f,double qdot0, double sdot0, double dt, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d, const int Nx) {

    double aW, aP, aE, aP0;
    double dxW, dxP, dxP0, dxE;
    double Uw, Ue;    
    double qdot_in;
    double qdot_rad;
    double Tw;

    // Compute "a" coefficient
    // --------------------------------------
    
    for (int j = 1; j < Nx - 1; j++) {

        // Precompute material properties k,Cp,alpha
        //---------------------------------------------------------------------------------------------------            
        
        alpha = k_s / (rho_s * Cp_s); 
        //-------------------------------------------------------------------------------------------------------------------------

        // Compute "a" coefficient
        // --------------------------------------
        dxW = x[j] - x[j - 1]; // west cell length
        dxP = x[j + 1] - x[j]; // central cell length
        Uw = (x[j] - x0[j]) / dt; // velocity of the west boundary of the P cell
        aW = Uw / (2 * alpha) - 2 / (dxW + dxP);
        a[j] = aW;

        // Compute "b" coefficient
        // ------------------------------------
        dxW = x[j] - x[j - 1]; // west cell length
        dxP = x[j + 1] - x[j]; // central cell length
        dxE = x[j + 2] - x[j + 1]; // east cell length
        Uw = (x[j] - x0[j]) / dt; // velocity of the west boundary of the P cell
        Ue = (x[j + 1] - x0[j + 1]) / dt; // velocity of the east boundary of the P cell
        aP = 2 / (dxE + dxP) + 2 / (dxW + dxP) - Ue / (2 * alpha) + Uw / (2 * alpha) + dxP / (dt * alpha);
        b[j] = aP;

        // Compute "c" coefficient
        // ---------------------------------
        dxE = x[j + 2] - x[j + 1]; //  east cell length
        dxP = x[j + 1] - x[j]; // central cell length
        Ue = (x[j + 1] - x0[j + 1]) / dt; // velocity of the east boundary of the cell
        aE = -2 / (dxE + dxP) - Ue / (2 * alpha);
        c[j] = aE;

        // Compute "d" coefficient (d=B*f0+C)
        // --------------------------------------
        dxP0 = x0[j + 1] - x0[j]; // central cell length at previous time
        aP0 = dxP0 / (dt * alpha);
        d[j] = aP0 * f0[j];

    }

    a[0] = 0;
    a[Nx - 1] = 1;
    // -------------------------------------
    
    b[0] = -1;    
    b[Nx - 1] = -1;
    // -----------------------------------------
        
    c[0] = 1;
    c[Nx - 1] = 0;
    //------------------------------------------
     
    
    Tw = f[1]; // current wall temperature
    qdot_rad = Epsilon * sigma * (pow(Tw,4) - pow(T_inf,4));
    double dh_sub = (2.5e-6 * Tw * Tw - 2.2e-3 * Tw + 0.773)*1e6;
    double qdot_sub = mdot_w_out * dh_sub;
    qdot_in = -(qdot0 - qdot_sub - qdot_rad) * (x[2] - x[0]) / (2 * k_s);   

    d[0] = qdot_in;    
    d[Nx - 1] = 0;
    //--------------------------------------

}
// Solver linear system of equations using Thomas algorithm
void  solidCond_Ray::SolveTDMA(vector<double>& f, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d, const int Nx) {

    // Solve system of eq. A*T=B*T_t0+C

    // Compute c_star
    vector<double> c_star(Nx - 1);
    c_star[0] = c[0] / b[0];

    for (int i = 1; i < Nx - 1; i++) {
        c_star[i] = c[i] / (b[i] - c_star[i - 1] * a[i]);
    }

    // Compute d_star
    vector<double> d_star(Nx);
    d_star[0] = d[0] / b[0];

    for (int i = 1; i < Nx; i++) {
        d_star[i] = (d[i] - d_star[i - 1] * a[i]) / (b[i] - c_star[i - 1] * a[i]);
    }

    // Compute solution vector f
    f[Nx - 1] = d_star[Nx - 1];
    for (int i = Nx - 2; i >= 0; i--) {
        f[i] = d_star[i] - c_star[i] * f[i + 1];  // update the solution vector f
    }


}


void solidCond_Ray::EvaluateTemp( vector<double>& f, double qdot0, double sdot0, double dt, const int Nx) {

    vector<double> a(Nx + 2), b(Nx + 2), c(Nx + 2), d(Nx + 2);

    //Generate a,b,c,d vectors for TDMA Algorithm //
    Get_abcd_coeff(f, qdot0, sdot0, dt, a, b, c, d, Nx + 2);
    //--------------------------------------------------

    //Solve linear system of equations and update temperature solution f0
    SolveTDMA(f, a, b, c, d, Nx + 2);
    //-----------------------------------------------------------------------


}


void solidCond_Ray::init(double T0, int numPts) {

    GenerateGeomGrid(numPts + 3); // Generate initial grid with geometric spacing
    //GenerateUniformGrid(numPts + 3);    

    for (int ni = 0; ni < numPts+2; ++ni) {        

        f0.push_back(T0);
    }
  
    
}


