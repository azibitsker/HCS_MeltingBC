#ifndef SOLIDCONDRAY_H_INCLUDED
#define SOLIDCONDRAY_H_INCLUDED

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <time.h>       /* time */

using namespace std;

class solidCond_Ray {

public:

    double L0 = 0.01; // initial material length
    double F = 1.3; // cell size increase ratio
    double k = 0.4, rho = 1500, Cp = 1500; // material properties
    double Tm = 1100; // [K] material melting temperature
    double Qstar = 1e6; // [J/kgK] heat of ablation
    double Eps_T = 1e-3;
    double alpha = k / (rho * Cp); // thermal diffusivity
    vector<double> f0,f, x0,x;    
    double sdot_out;   
    double eta = 0.5;
  

public:
    void HeatCondSolver(double qdot_in,int Nx, double dt);
    void EvaluateTemp(double qdot0, double sdot0, double dt, const int Nx);
    void GenerateGeomGrid(const int Nx);
    void GenerateUniformGrid(int Nx);
    void ContractGridBoundaries(double sdot0, double dt, const int Nx);
    void moveInternalMesh(int N_iter);
    void Get_abcd_coeff(double qdot0, double sdot0, double dt, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d, const int Nx);
    void SolveTDMA(vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d, const int Nx);
    bool CheckMelting(vector<double>& f);
    double GetMeltedLength(vector<double>& f);
    double GetRecessionRate(vector<double>& f, double dt);
    void init(double T0, int numPts);
    
};


#endif //SOLIDCONDRAY_H_INCLUDED

