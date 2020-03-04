#include "SolidCondRay.h"

double SolidCond_Ray::HeatCondSolver(double* f0, double* x0, double* x, double dt, double qdot0, double sdot0, const int Nx, int ind_t) {

    double alpha; // thermal diffusivity
    int Nmelt; // number of cells that have reached the melting temperature
    double sdot = 0; // computed recession rate
    double sdot_diff = 1;
    vector<double> a(Nx + 2), b(Nx + 2), c(Nx + 2), d(Nx + 2);
    double Tw;
    double sdot_saver[100] = { 0 };
    int jj = 2;

    alpha = k / (rho * Cp);

    // f0 and sdot are coupled. The loop converges when abs(sdot0-sdot)<Epsilon
    while (sdot_diff > Eps_sdot) {

        //Generate a,b,c,d vectors for TDMA Algorithm //
        Get_abcd_coeff(x0, x, f0, qdot0, sdot0, alpha, dt, a, b, c, d, Nx + 2, ind_t);
        //--------------------------------------------------

        //Solve linear system of equations and update temperature solution f0
        SolveTDMA(f0, a, b, c, d, Nx + 2);
        //-----------------------------------------------------------------------

        // Count the number of cells that have reached the melting temperature
        Nmelt = GetMeltedCells(f0, Nx + 2);

        // Compute recession rate based on the number of material cells that have reached the melting temperature
        sdot = GetRecessionRate(x0, dt, Nmelt);  // [m/s] updated recession rate

        sdot_diff = abs(sdot - sdot0);


        if (sdot_diff > Eps_sdot) {

            sdot_saver[jj] = sdot;
            jj += 1;
            Tw = f0[1];
            //cout<<sdot0<<" ";
            sdot0 = sdot; // update the guess for recession rate
            ContractGridBoundaries(x0, x, sdot0, dt, Nx + 3, ind_t); // array x is updated
        }
    }

    // cout<<endl;

    sdot0 = sdot_saver[jj - 2];
    ContractGridBoundaries(x0, x, sdot0, dt, Nx + 3, ind_t); // array x is updated
    // Update  x0
    for (int j = 0; j < Nx + 3; j++) { x0[j] = x[j]; }

    return sdot0;
}

// Generate initial material grid including nodes for the two ghost cells
void SolidCond_Ray::GenerateGrid(double* x0, const int Nx) {
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

// Update grid boundaries based on the computed recession rate sdot0
void SolidCond_Ray::ContractGridBoundaries(double* x0, double* x, double sdot0, double dt, const int Nx, int ind_t) {

    double dx_m = sdot0 * dt; // total length of the melted cells
    double L;

    L = x0[Nx - 2] - x0[1]; // current material length

//------- Update location of each cell boundary relative to the initial location at t=0------------------
    double dx0, dxj;

    x[1] = x0[1] + dx_m; // update location of the recession front
    double Recession = x[1];

    if (Recession >= L0) {
        cout << "The whole material was melted away" << "time = " << ind_t * dt << "[sec]" << endl;
        //return;
    }

    for (int j = 2; j < Nx - 1; j++) {

        dx0 = x0[j] - x0[j - 1]; // compute material initial cell length
        dxj = dx0 - dx0 / L * dx_m; // compute the new cell length

        x[j] = x[j - 1] + dxj; // update location of each cell boundary

    }
    x[0] = x[1] - (x[2] - x[1]); // location of the 1st ghost cell boundary
    x[Nx - 1] = x[Nx - 2] + dxj;  // location of the last ghost cell boundary

}

// Get a,b,c coefficients for Thomas algorithm
void SolidCond_Ray::Get_abcd_coeff(double* x0, double* x, double* f0, double qdot0, double sdot0, double alpha, double dt, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d, const int Nx, int ind_t) {

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
    if (qdot_in > 0) {
        cout << "int_t = " << ind_t << ": " "qdot0 < rho*Qstar*sdot0" << endl;
        //return;
    }
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
void  SolidCond_Ray::SolveTDMA(double* f0, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d, const int Nx) {

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

    f0[Nx - 1] = d_star[Nx - 1];
    for (int i = Nx - 2; i >= 0; i--) {
        f0[i] = d_star[i] - c_star[i] * f0[i + 1];  // update the solution vector f0
    }


}

int SolidCond_Ray::GetMeltedCells(double* f0, const int Nx) {

    double Tj;
    int jm = 0;

    // Find the cells that have reached the melting temperature and count the number of such cells
    for (int j = 1; j < Nx - 1; j++) {

        Tj = f0[j];

        if (Tj >= Tm) {

            jm += 1;
        }
        else {

            break;
        }

    }

    return jm;
}


double SolidCond_Ray::GetRecessionRate(double* x0, double dt, int Nmelt) {


    double sdot;


    // Compute the amount of recession
    double dx_m;

    dx_m = x0[Nmelt + 1] - x0[1];
    sdot = dx_m / dt;

    return sdot;
}



