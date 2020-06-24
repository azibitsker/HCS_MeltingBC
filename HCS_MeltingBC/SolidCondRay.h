#ifndef SOLIDCONDRAY_H_INCLUDED
#define SOLIDCONDRAY_H_INCLUDED

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

class solidCond_Ray {

public:

    double F = 1.3; // cell size increase ratio
    double k = 0.4, rho = 1500, Cp = 1500; // material properties
    double Tm = 1100; // [K] material melting temperature
    double Qstar = 1e6; // [J/kgK] heat of ablation
    double Eps_T = 1e-3;
    double alpha = k / (rho * Cp); // thermal diffusivity
    vector<double> f0, x0,x;
    double sdot_out;   
    double eta = 0.5;
    int m = 2; //  0: Cartesian, 1: Cylindrical, 2: Spherical
    double Rout=0.04; // outer radius of cylinder or sphere
    double Rin=0.01; // inner radius of ylinder or sphere
    double L0 = Rout - Rin; // initial material length
    int numRays = 1;
    double pi = 3.14159265359;
  

public:
    void HeatCondSolver(double qdot_in,int Nx, double dt, ofstream& IterationsFile);
    void EvaluateTemp(vector<double>& f, double qdot0, double sdot0, double dt, const int Nx);
    void GenerateGeomGrid(const int Nx);
    void GenerateUniformGrid(int Nx);
    void ContractGridBoundaries(double sdot0, double dt, const int Nx);
    void Get_abcd_coeff(double qdot0, double sdot0, double dt, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d, const int Nx);
    void SolveTDMA(vector<double>& f, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d, const int Nx);
    bool CheckMelting(vector<double>& f);
    double GetMeltedLength(vector<double>& f);
    double GetRecessionRate(vector<double>& f, double dt);
    void init(double T0, int numPts);
    
};


#endif //SOLIDCONDRAY_H_INCLUDED

