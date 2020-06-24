#include "solidCondRay.h"
#include "materialResponse.h"

void solidCond_Ray::HeatCondSolver(double qdot_in, int Nx, double dt) {

    double sdot0 = 0;   
    

    // Evaluate temperature at the next time step, vector f is updated
    EvaluateTemp(qdot_in, sdot0, dt, Nx);   
    
    // Update  x0
    //for (int j = 0; j < Nx + 3; j++) { x0[j] = x[j]; }    
   
    // Update  f0
    //for (int j = 0; j < Nx + 2; j++) { f0[j] = f[j]; }
    
        
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
void solidCond_Ray::ContractGridBoundaries(double sdot0, double dt, const int Nx)
{

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

void solidCond_Ray::moveInternalMesh(int N_iter)
{

    int Nx = x.size();
    double dx_f;

    /* initialize random seed: */
    srand(time(NULL));    

//------- Move randomly internal face of the mesh------------------    

    for (int i = 0; i < N_iter; ++i)
    {
        for (int j = 2; j < Nx - 2; j++)
        {
            dx_f = (((double)rand() / (RAND_MAX)) * 2 - 1) / 2e4;
            x[j] = x[j] + dx_f;
        }
    }


    

}

// Get a,b,c coefficients for Thomas algorithm
void solidCond_Ray::Get_abcd_coeff(double qdot0, double sdot0, double dt, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d, const int Nx) {

    double aW, aP, aE, aP0;
    double dxW, dxP, dxP0, dxE;
    double Uw, Ue;
    double fw, fe;

    // Compute "a" coefficient
    // --------------------------------------
    a[0] = 0;
    for (int j = 1; j < Nx - 1; j++) {

        dxW = x[j] - x[j - 1]; // west cell length
        dxP = x[j + 1] - x[j]; // central cell length
        Uw = (x[j] - x0[j]) / dt; // velocity of the west boundary of the P cell
        fw = 1 / (1 + dxW / dxP);
        aW = Uw*fw - 2*alpha / (dxW + dxP);

        a[j] = aW;

    }
    a[Nx - 1] = 1;
    // -------------------------------------

    // Compute "b" coefficient
    // ------------------------------------
    double dx0 = x[1] - x[0]; // left ghost cell length
    double dx1 = x[2] - x[1]; // adjacent material cell length
    b[0] = 1;

    for (int j = 1; j < Nx - 1; j++) {

        dxW = x[j] - x[j - 1]; // west cell length
        dxP = x[j + 1] - x[j]; // central cell length
        dxE = x[j + 2] - x[j + 1]; // east cell length
        Uw = (x[j] - x0[j]) / dt; // velocity of the west boundary of the P cell
        Ue = (x[j + 1] - x0[j + 1]) / dt; // velocity of the east boundary of the P cell
        fe = 1 / (1 + dxP / dxE);
        fw = 1 / (1 + dxW / dxP);

        aP = 2*alpha / (dxE + dxP) + 2*alpha / (dxW + dxP) - (Ue*fe-(1-fw)*Uw) + dxP / dt;

        b[j] = aP;

    }
    b[Nx - 1] = -1;
    // -----------------------------------------

    // Compute "c" coefficient
    // ---------------------------------
    c[0] = -1;

    for (int j = 1; j < Nx - 1; j++) 
    {

        dxE = x[j + 2] - x[j + 1]; //  east cell length
        dxP = x[j + 1] - x[j]; // central cell length
        Ue = (x[j + 1] - x0[j + 1]) / dt; // velocity of the east boundary of the cell
        fe = 1 / (1 + dxP / dxE);

        aE = -2*alpha / (dxE + dxP) - (1-fe)*Ue;

        c[j] = aE;

    }
    c[Nx - 1] = 0;

    // Compute "d" coefficient (d=B*f0+C)
    // --------------------------------------
    double qdot_in;

    qdot_in = 0.;
    
    d[0] = qdot_in;

    for (int j = 1; j < Nx - 1; j++) {

        dxP0 = x0[j + 1] - x0[j]; // central cell length at previous time

        aP0 = dxP0 / (dt);

        d[j] = aP0 * f0[j];

    }
    d[Nx - 1] = 0;
    //--------------------------------------

}
// Solver linear system of equations using Thomas algorithm
void  solidCond_Ray::SolveTDMA(vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d, const int Nx) {

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

double solidCond_Ray::GetMeltedLength(vector<double>& f) {

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


double solidCond_Ray::GetRecessionRate(vector<double>& f, double dt) {

    double sdot;
    double dx_melt; // number of cells that have reached the melting temperature

    // Count the number of cells that have reached the melting temperature
    dx_melt = GetMeltedLength(f);

    // Compute the amount of recession
    sdot = dx_melt / dt;

    return sdot;
}


void solidCond_Ray::EvaluateTemp(double qdot0, double sdot0, double dt, const int Nx) {

    vector<double> a(Nx + 2), b(Nx + 2), c(Nx + 2), d(Nx + 2);

    //Generate a,b,c,d vectors for TDMA Algorithm //
    Get_abcd_coeff(qdot0, sdot0, dt, a, b, c, d, Nx + 2);
    //--------------------------------------------------

    //Solve linear system of equations and update temperature solution f0
    SolveTDMA(a, b, c, d, Nx + 2);
    //-----------------------------------------------------------------------


}

bool solidCond_Ray::CheckMelting(vector<double>& f) {

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


void solidCond_Ray::init(double T0, int numPts) {

    //GenerateGeomGrid(numPts + 3); // Generate initial grid with uniform spacing
    GenerateUniformGrid(numPts + 3);

    f.resize(numPts + 2);

    for (int ni = 0; ni < numPts+2; ++ni) {     

        f0.push_back(T0);
    }
  
    
}


