// Constant array size

/* 
#ifndef SOLIDCONDRAY_H_INCLUDED
#define SOLIDCONDRAY_H_INCLUDED

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

class SolidCond_Ray {

    double L0 = 0.02; // initial material length
    double F = 1.001; // cell size increase ratio
    double k = 0.402, rho = 280, Cp = 983.9; // material properties
    double Tm = 800; // [K] material melting temperature
    double Qstar = 1e5; // [J/kgK] heat of ablation
    double Eps_T = 0.5;
    double alpha = k / (rho * Cp); // thermal diffusivity


public:
    double HeatCondSolver(double* f0, double* x0, double* x, double dt, double qdot0, const int Nx, int t, ofstream& IterationsFile, ofstream& NcellsFile);
    void GenerateGeomGrid(double* x0, const int Nx);
    void GenerateUniformGrid(double* x0, const int Nx);
    void ContractGridBoundaries(double* x0, double* x, double sdot0, double dt, const int Nx, int ind_t);
    void Get_abcd_coeff(double* x0, double* x, double* f0, double qdot0, double sdot0, double dt, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d, const int Nx, int ind_t);
    void SolveTDMA(vector<double>& f, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d, const int Nx);
    void EvaluateTemp(double* x0, double* x, double* f0, vector<double>& f, double qdot0, double sdot0, double dt, int ind_t, const int Nx);
    bool CheckMelting(vector<double>& f);
    double GetMeltedLength(vector<double>& f, double* x0);
    double GetRecessionRate(vector<double>& f, double* x0, double dt);


};
#endif // SOLIDCONDRAY_H_INCLUDED

*/

// Constant increment size

#ifndef SOLIDCONDRAY_H_INCLUDED
#define SOLIDCONDRAY_H_INCLUDED

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

class SolidCond_Ray {

public:

    double L0 = 0.02; // initial material length
    double F = 1.001; // cell size increase ratio
    double k = 0.402, rho = 280, Cp = 983.9; // material properties
    double Tm = 800; // [K] material melting temperature
    double Qstar = 1e5; // [J/kgK] heat of ablation
    double Eps_T = 1;
    double alpha= k / (rho * Cp); // thermal diffusivity
    

public:
    double HeatCondSolver(vector<double>& f0, vector<double>& x0, vector<double>& x, double dt, double qdot0, const int Nx, int t,ofstream & IterationsFile,ofstream & NcellsFile);
    void GenerateGeomGrid(vector<double>& x0, const int Nx);
    void GenerateUniformGrid(vector<double>& x0, double dx);
    void ContractGridBoundaries(vector<double>& x0, vector<double>& x, double sdot0, double dt, const int Nx, int ind_t);
    void Get_abcd_coeff(vector<double>& x0, vector<double>& x, vector<double>& f0, double qdot0, double sdot0, double dt, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d, const int Nx, int ind_t);
    void SolveTDMA(vector<double>& f, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d, const int Nx);
    void EvaluateTemp(vector<double>& x0, vector<double>& x, vector<double>& f0, vector<double> &f, double qdot0, double sdot0, double dt, int ind_t, const int Nx);
    bool CheckMelting(vector<double>& f);
    double GetMeltedLength(vector<double>& f, vector<double>& x0);
    double GetRecessionRate(vector<double>& f, vector<double>& x0,double dt);
	
	
};
#endif // SOLIDCONDRAY_H_INCLUDED
