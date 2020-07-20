#include "solidCondRay.h"
#include "materialResponse.h"


void solidCond_Ray::HeatCondSolver(double qdot_in, int Nx, double dt) {
    
    
    double sdot=0; // computed recession rate  
    double sdot0 = 0;
    double sdot_g = 0;
    double Tl = 0;
    double Tr = 0;
    double Twall = 0;
    double TT = 0;
    vector<double> f(Nx + 2);
    vector<double> grid(Nx + 3);    
        

    if (IterSdot)
    {

        for (int j = 0; j < Nx + 3; j++) { grid[j] = x[j]; }
        EvaluateTemp(f, grid, qdot_in, sdot0, dt, Nx);
        //-----------------------------------------------------------------
        // Convergence on Tabl by halving the guess of sdot        

        // Check if at least one of the material cells has reached melting temperature
        if (CheckMelting(f)) {

            // Interval Halving method to find sdot such that Tw==Tm

            Tl = (f[0] + f[1]) / 2;
            sdot = GetRecessionRate(f, dt);

            // Iterative search for the true recession that assures Twall=Tmelt             

            ContractGridBoundaries(sdot, dt, Nx + 3); // array x is updated
            for (int j = 0; j < Nx + 3; j++) { grid[j] = x[j]; }
            EvaluateTemp(f,grid, qdot_in, sdot, dt, Nx); // reevaluate temperature with the updated sdot0 and x; vector f is updated
            Tr = (f[0] + f[1]) / 2;
            TT = (Tl - Tm) * (Tr - Tm);

            // Perform a check if sdot0 and sdot are from two sides of the true sdot solution//
            // If not, halve sdot until the requirement of interval halving method is fulfilled 
            if (TT > 0)
            {

                while (TT > 0)
                {
                    sdot = sdot / 2; // halve sdot until requirement of the method is fulfilled
                    ContractGridBoundaries(sdot, dt, Nx + 3); // array x is updated
                    for (int j = 0; j < Nx + 3; j++) { grid[j] = x[j]; }
                    EvaluateTemp(f,grid, qdot_in, sdot, dt, Nx); // reevaluate temperature with the updated sdot0 and x; vector f is updated
                    Tr = (f[0] + f[1]) / 2;
                    TT = (Tl - Tm) * (Tr - Tm);
                }

            }

            sdot_g = (sdot0 + sdot) / 2; // initial guess, sdot_true is between sdot0 and sdot
            ContractGridBoundaries(sdot_g, dt, Nx + 3); // array x is updated
            for (int j = 0; j < Nx + 3; j++) { grid[j] = x[j]; }
            EvaluateTemp(f,grid, qdot_in, sdot_g, dt, Nx); // reevaluate temperature with the updated sdot0 and x; vector f is updated
            Tr = (f[0] + f[1]) / 2;
            Twall = Tr;

            // Loop until the guess of sdot0 unsures enough recession to happen to reduce the q_in and cause Twall=Tmelting
            while (Twall < Tm - Eps_T || Twall>Tm) {

                TT = (Tl - Tm) * (Tr - Tm);

                if (TT < 0) {
                    sdot = sdot_g;
                }
                else {
                    sdot0 = sdot_g;
                    Tl = (f[0] + f[1]) / 2;
                }

                sdot_g = (sdot0 + sdot) / 2; // guess sdot0 for the next evaluation of temperature profile
                ContractGridBoundaries(sdot_g, dt, Nx + 3); // array x is updated
                for (int j = 0; j < Nx + 3; j++) { grid[j] = x[j]; }
                EvaluateTemp(f,grid, qdot_in, sdot_g, dt, Nx); // reevaluate temperature with the updated sdot0 and x; vector f is updated
                Tr = (f[0] + f[1]) / 2;
                Twall = Tr;
            }

            // Update  x0
            for (int j = 0; j < Nx + 3; j++) { x0[j] = x[j]; }

        }

        // Update  f0
        for (int j = 0; j < Nx + 2; j++) { f0[j] = f[j]; }
        sdot_out = sdot_g;
        

    }

    else if (SteadyStateSdot) // Imitate recession at steady state by passing constant sdot from the steady state solution
    {
        // Evaluate temperature at the next time step, vector f is updated
        if (CheckMelting(f0))
        {
            sdot0 = knownRecessionRate;
            ContractGridBoundaries(sdot0, dt, Nx + 3); // array x is updated            
        }
        
        for (int j = 0; j < Nx + 3; j++) { grid[j] = x0[j]; }
        EvaluateTemp(f,grid, qdot_in, sdot0, dt, Nx);
        // Update  f0
        for (int j = 0; j < Nx + 2; j++) { f0[j] = f[j]; }
        //Update x0
        for (int j = 0; j < Nx + 3; j++) { x0[j] = x[j]; }  
        sdot_out = sdot0;
    }

    else if (PureCond) // pure conduction
    {

        for (int j = 0; j < Nx + 3; j++) { grid[j] = x0[j]; }
        EvaluateTemp(f, grid, qdot_in, sdot0, dt, Nx);
        // Update  f0
        for (int j = 0; j < Nx + 2; j++) { f0[j] = f[j]; }
        sdot_out = sdot0;
    }
    else
    {
        cout << "No simulation type has been defined" << endl;
    } 
    //--------------------------------------------------------------------------------
        
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

// Get a,b,c coefficients for Thomas algorithm
void solidCond_Ray::Get_abcd_coeff(vector<double>& grid, double qdot0, double sdot0, double dt, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d, const int Nx) {

    double aW, aP, aE, aP0;
    double dxW, dxP, dxE;
    double Uw, Ue;
    double fw, fe;
    double Aw, Ae;
    double Vp, Vp0;

    // Compute "a" coefficient
    // --------------------------------------
    a[0] = 0;
    for (int j = 1; j < Nx - 1; j++) {

        dxW = grid[j] - grid[j - 1]; // west cell length
        dxP = grid[j + 1] - grid[j]; // central cell length
        Uw = (x[j] - x0[j]) / dt; // velocity of the west boundary of the P cell
        fw = 1 / (1 + dxW / dxP);
        Aw = pow((2 * (Rout-x[j])), m) * pow(pi, 0.5 * m * (3 - m))/numRays;

        aW = fw*Uw*Aw - 2*alpha*Aw / (dxW + dxP);
        a[j] = aW;

    }
    a[Nx - 1] = 1;
    // -------------------------------------

    // Compute "b" coefficient
    // ------------------------------------
    double dx0 = grid[1] - grid[0]; // left ghost cell length
    double dx1 = grid[2] - grid[1]; // adjacent material cell length
    b[0] = 2 * k / (dx0 + dx1);

    for (int j = 1; j < Nx - 1; j++) {

        dxW = grid[j] - grid[j - 1]; // west cell length
        dxP = grid[j + 1] - grid[j]; // central cell length
        dxE = grid[j + 2] - grid[j + 1]; // east cell length
        Uw = (x[j] - x0[j]) / dt; // velocity of the west boundary of the P cell
        Ue = (x[j + 1] - x0[j + 1]) / dt; // velocity of the east boundary of the P cell
        fe = 1 / (1 + dxP / dxE);
        fw = 1 / (1 + dxW / dxP);
        Aw = pow((2 * (Rout - x[j])), m) * pow(pi, 0.5 * m * (3 - m)) / numRays;
        Ae = pow((2 * (Rout - x[j+1])), m) * pow(pi, 0.5 * m * (3 - m)) / numRays;
        Vp = pow(2, m) * pow(pi, 0.5 * m * (3 - m)) / (m + 1) * (pow((Rout - x[j]), m + 1) - pow((Rout - x[j + 1]), m + 1))/numRays;

        aP = 2*alpha*Ae / (dxE + dxP) + 2*alpha*Aw / (dxW + dxP) - (Ue*fe*Ae-(1-fw)*Uw*Aw) + Vp / dt;
        b[j] = aP;

    }
    b[Nx - 1] = -1;
    // -----------------------------------------

    // Compute "c" coefficient
    // ---------------------------------
    c[0] = -2 * k / (dx0 + dx1);

    for (int j = 1; j < Nx - 1; j++) 
    {

        dxE = grid[j + 2] - grid[j + 1]; //  east cell length
        dxP = grid[j + 1] - grid[j]; // central cell length
        Ue = (x[j + 1] - x0[j + 1]) / dt; // velocity of the east boundary of the cell
        fe = 1 / (1 + dxP / dxE);
        Ae = pow((2 * (Rout - x[j + 1])), m) * pow(pi, 0.5 * m * (3 - m)) / numRays;

        aE = -2*alpha*Ae / (dxE + dxP) - (1-fe)*Ue*Ae;
        c[j] = aE;

    }
    c[Nx - 1] = 0;

    // Compute "d" coefficient (d=B*f0+C)
    // --------------------------------------
    double qdot_in;
    double Awall;
    Awall = pow((2 * (Rout - x[1])), m) * pow(pi, 0.5 * m * (3 - m)) / numRays;

    qdot_in = qdot0/Awall - rho * Qstar * sdot0;
   
    d[0] = qdot_in;

    for (int j = 1; j < Nx - 1; j++) {
        
        Vp0 = pow(2, m) * pow(pi, 0.5 * m * (3 - m)) / (m + 1) * (pow((Rout - x0[j]), m + 1) - pow((Rout - x0[j + 1]), m + 1)) / numRays;

        aP0 = Vp0 / (dt);
        d[j] = aP0 * f0[j];

    }
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

double solidCond_Ray::GetMeltedLength(vector<double>& f) {

    int jm = 0;    
    double dx_melt = 0;
    double Ti = (f[0]+f[1])/2;
    double x_m;

    while (Ti > Tm) {

        jm += 1;
        Ti = f[jm];
    }

    // Compute melted length using interpolation to find location of Tm
    if (jm > 0) {

        if (jm == 1)
        {
            double xc1 = x0[jm]; // at the wall
            double xc2 = (x0[jm] + x0[jm + 1]) / 2;
            Ti = (f[0] + f[1]) / 2;
            double m = (Ti - f[jm]) / (xc1 - xc2); // linear slope between Tj<Tm and Tj>Tm
            x_m = (Tm - f[jm]) / m + xc2; // location of Tm
            double dx_inc = x_m - xc1; // incremental distance between Tm and Tj
        }
        else
        {
            double xc1 = (x0[jm] + x0[jm + 1]) / 2;
            double xc2 = (x0[jm + 1] + x0[jm + 2]) / 2;
            double m = (f[jm] - f[jm + 1]) / (xc1 - xc2); // linear slope between Tj<Tm and Tj>Tm
            x_m = (Tm - f[jm + 1]) / m + xc2; // location of Tm
            double dx_inc = x_m - xc1; // incremental distance between Tm and Tj
        }        

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


void solidCond_Ray::EvaluateTemp(vector<double>& f, vector<double>& grid, double qdot0, double sdot0, double dt, const int Nx) {

    vector<double> a(Nx + 2), b(Nx + 2), c(Nx + 2), d(Nx + 2);

    //Generate a,b,c,d vectors for TDMA Algorithm //
    Get_abcd_coeff(grid, qdot0, sdot0, dt, a, b, c, d, Nx + 2);
    //--------------------------------------------------

    //Solve linear system of equations and update temperature solution f0
    SolveTDMA(f, a, b, c, d, Nx + 2);
    //-----------------------------------------------------------------------


}

bool solidCond_Ray::CheckMelting(vector<double>& f) {

    double Twall = 0;

    Twall = (f[0] + f[1]) / 2;

    if (Twall>Tm) {
        return true;
    }
    else {
        return false;
    }

}


void solidCond_Ray::init(double T0, int numPts) {

    if (GeomGrid) GenerateGeomGrid(numPts + 3);  
    else GenerateUniformGrid(numPts + 3);
    

    for (int ni = 0; ni < numPts+2; ++ni) {       

        f0.push_back(T0);
    }
  
    
}


