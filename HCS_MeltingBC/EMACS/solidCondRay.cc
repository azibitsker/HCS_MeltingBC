#include "SolidCondRay.h"

double SolidCond_Ray::HeatCondSolver(double* f0, double* x0, double* x, double dt, double qdot0, const int Nx, int ind_t, ofstream& IterationsFile, ofstream& NcellsFile) {


    double sdot; // computed recession rate  
    //double dT;
    double sdot0 = 0;
    vector<double> f(Nx + 2);



    // Evaluate temperature at the next time step, vector f is updated
    EvaluateTemp(x0, x, f0, f, qdot0, sdot0, dt, ind_t, Nx);

    // Check if at least one of the material cells has reached melting temperature
    if (CheckMelting(f)) {

        int N_iter = 0;

        // Loop until the guess of sdot0 unsures enough recession to happen to reduce the q_in and cause Twall=Tmelting
        while (f[1] < Tm - Eps_T || f[1]>Tm) {

            N_iter++;

            // Compute recession rate based on the number of material cells that have reached the melting temperature       

            sdot = GetRecessionRate(f, x0, dt);
            sdot0 = (sdot0 + sdot) / 2; // guess sdot0 for the next evaluation of temperature profile

            // Contract grid x0 boundaries based on the comouted sdot0
            ContractGridBoundaries(x0, x, sdot0, dt, Nx + 3, ind_t); // array x is updated

            EvaluateTemp(x0, x, f0, f, qdot0, sdot0, dt, ind_t, Nx); // reevaluate temperature with the updated sdot0 and x; vector f is updated


        }

        IterationsFile << ind_t << " " << N_iter;
        IterationsFile << endl;


        // count how many cells are melted at each time step
        double N_ratio;
        double dx_melt = sdot0 * dt;

        N_ratio = dx_melt / (x0[2] - x0[1]);
        NcellsFile << N_ratio;
        NcellsFile << endl;

        // Update  x0
        for (int j = 0; j < Nx + 3; j++) { x0[j] = x[j]; }

    }

    // Update  f0
    for (int j = 0; j < Nx + 2; j++) { f0[j] = f[j]; }


    return sdot0;
}

// Generate initial material grid including nodes for the two ghost cells using geometric series
void SolidCond_Ray::GenerateGeomGrid(double* x0, const int Nx) {
    // Generate the initial material grid. Fill out the array x0 with location of the
    // cells boundaries using geometric series

    double dx0;
    double dxj;

    dx0 = L0 * (1 - F) / (1 - pow(F, Nx - 3)); // initial grid cell size
    x0[0] = -dx0;
    x0[1] = 0; // 1st ghost cell - its size is equal to the size of the first material cell

    for (int j = 2; j < Nx - 1; j++) {

        dxj = pow(F, j - 2) * dx0;
        x0[j] = x0[j - 1] + dxj; // space the cell boundaries using geometric series
    }

    x0[Nx - 1] = x0[Nx - 2] + dxj; // last ghost cell - its size is equal to the size of the last material cell

    //cout<< "L0 = "<< x0[Nx-2]-x0[1]<< endl;

}

void SolidCond_Ray::GenerateUniformGrid(double* x0, const int Nx) {
    // Generate the initial material grid. Fill out the array x0 with location of the
    // cells boundaries using uniform series

    double dx0;
    //double dxj;

    dx0 = L0 / (Nx - 3); // initial grid cell size
    x0[0] = -dx0;

    for (int j = 1; j < Nx - 1; j++) {

        x0[j] = x0[j - 1] + dx0; // space the cell boundaries using uniform spacing
    }

    x0[Nx - 1] = x0[Nx - 2] + dx0; // last ghost cell - its size is equal to the size of the last material cell

    //cout<< "L0 = "<< x0[Nx-2]-x0[1]<< endl;

}

// Update grid boundaries based on the computed recession rate sdot0
void SolidCond_Ray::ContractGridBoundaries(double* x0, double* x, double sdot0, double dt, const int Nx, int ind_t) {

    double dx_melt = sdot0 * dt; // total length of the melted cells
    double L;

    L = x0[Nx - 2] - x0[1]; // current material length

//------- Update location of each cell boundary relative to the initial location at t=0------------------
    double dx0, dxj;

    x[1] = x0[1] + dx_melt; // update location of the recession front
    double Recession = x[1];

    if (Recession >= L0) {
        cout << "The whole material was melted away" << "time = " << ind_t * dt << "[sec]" << endl;
        //return;
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
void SolidCond_Ray::Get_abcd_coeff(double* x0, double* x, double* f0, double qdot0, double sdot0, double dt, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d, const int Nx, int ind_t) {

    double aW, aP, aE, aP0;
    double dxW, dxP, dxP0, dxE;
    double Uw, Ue;

    // Compute "a" coefficient
    // --------------------------------------
    a[0] = 0;
    for (int j = 1; j < Nx - 1; j++) {

        dxW = x[j] - x[j - 1]; // west cell length
        dxP = x[j + 1] - x[j]; // central cell length
        Uw = (x[j] - x0[j]) / dt; // velocity of the west boundary of the P cell

        aW = Uw / (2 * alpha) - 2 / (dxW + dxP);

        a[j] = aW;

    }
    a[Nx - 1] = 1;
    // -------------------------------------

    // Compute "b" coefficient
    // ------------------------------------
    b[0] = -1;

    for (int j = 1; j < Nx - 1; j++) {

        dxW = x[j] - x[j - 1]; // west cell length
        dxP = x[j + 1] - x[j]; // central cell length
        dxE = x[j + 2] - x[j + 1]; // east cell length
        Uw = (x[j] - x0[j]) / dt; // velocity of the west boundary of the P cell
        Ue = (x[j + 1] - x0[j + 1]) / dt; // velocity of the east boundary of the P cell

        aP = 2 / (dxE + dxP) + 2 / (dxW + dxP) - Ue / (2 * alpha) + Uw / (2 * alpha) + dxP / (dt * alpha);

        b[j] = aP;

    }
    b[Nx - 1] = -1;
    // -----------------------------------------

    // Compute "c" coefficient
    // ---------------------------------
    c[0] = 1;

    for (int j = 1; j < Nx - 1; j++) {

        dxE = x[j + 2] - x[j + 1]; //  east cell length
        dxP = x[j + 1] - x[j]; // central cell length
        Ue = (x[j + 1] - x0[j + 1]) / dt; // velocity of the east boundary of the cell

        aE = -2 / (dxE + dxP) - Ue / (2 * alpha);

        c[j] = aE;

    }
    c[Nx - 1] = 0;

    // Compute "d" coefficient (d=B*f0+C)
    // --------------------------------------
    double qdot_in;

    qdot_in = -(qdot0 - rho * Qstar * sdot0) * (x[2] - x[0]) / (2 * k);
    //if (qdot_in > 0) {
       // cout << "int_t = " << ind_t << ": " "qdot0 < rho*Qstar*sdot0" << endl;
        //return;
    //}
    d[0] = qdot_in;

    for (int j = 1; j < Nx - 1; j++) {

        dxP0 = x0[j + 1] - x0[j]; // central cell length at previous time

        aP0 = dxP0 / (dt * alpha);

        d[j] = aP0 * f0[j];

    }
    d[Nx - 1] = 0;
    //--------------------------------------

}
// Solver linear system of equations using Thomas algorithm
void  SolidCond_Ray::SolveTDMA(vector<double>& f, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d, const int Nx) {

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

double SolidCond_Ray::GetMeltedLength(vector<double>& f, double* x0) {

    int jm = 0;
    int j = 1;
    double dx_melt = 0;

    while (f[jm + 1] > Tm) {

        jm += 1;

    }

    // Compute melted length using interpolation to find location of Tm
    if (jm > 0) {

        double xc1 = (x0[jm] + x0[jm + 1]) / 2;
        double xc2 = (x0[jm + 1] + x0[jm + 2]) / 2;
        double m = (f[jm] - f[jm + 1]) / (xc1 - xc2); // linear slope between Tj<Tm and Tj>Tm
        double x_m = (Tm - f[jm + 1]) / m + xc2; // location of Tm
        double dx_inc = x_m - xc1; // incremental distance between Tm and Tj

        dx_melt = x_m - x0[1]; // total melted length

    }

    return dx_melt;
}


double SolidCond_Ray::GetRecessionRate(vector<double>& f, double* x0, double dt) {

    double sdot;
    double dx_melt; // number of cells that have reached the melting temperature

    // Count the number of cells that have reached the melting temperature
    dx_melt = GetMeltedLength(f, x0);

    // Compute the amount of recession
    sdot = dx_melt / dt;

    return sdot;
}


void SolidCond_Ray::EvaluateTemp(double* x0, double* x, double* f0, vector<double>& f, double qdot0, double sdot0, double dt, int ind_t, const int Nx) {

    vector<double> a(Nx + 2), b(Nx + 2), c(Nx + 2), d(Nx + 2);

    //Generate a,b,c,d vectors for TDMA Algorithm //
    Get_abcd_coeff(x0, x, f0, qdot0, sdot0, dt, a, b, c, d, Nx + 2, ind_t);
    //--------------------------------------------------

    //Solve linear system of equations and update temperature solution f0
    SolveTDMA(f, a, b, c, d, Nx + 2);
    //-----------------------------------------------------------------------


}

bool SolidCond_Ray::CheckMelting(vector<double>& f) {

    int jm = 0;

    while (f[jm + 1] > Tm) {

        jm += 1;

    }

    if (jm != 0) {
        return true;
    }
    else {
        return false;
    }

}
