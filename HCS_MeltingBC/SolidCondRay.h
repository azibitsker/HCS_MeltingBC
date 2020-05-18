#ifndef SOLIDCONDRAY_H_INCLUDED
#define SOLIDCONDRAY_H_INCLUDED

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

class solidCond_Ray {

public:

    double L0 = 0.03; // initial material length
    double F = 1.1; // cell size increase ratio
    //double k = 0.2, rho = 2000, Cp = 1000; // material properties
    double Tavg = 1000;
    double k;
    double rho = 1775;
    double Cp;
    double Tm = 800; // [K] material melting temperature
    double Qstar = 2e6; // [J/kgK] heat of ablation
    double Eps_T = 1e-3;
    double alpha; // thermal diffusivity
    double sdot_out;   
    double eta = 0.5;
    vector<double> f0, x0, x, mdot_k_out, yk_w_out;
    

    // Radiation parameters
    double Epsilon = 0.9;
    double sigma = 5.670374e-8; // [W/m^2K^4]
    double T_inf = 300; // ambient temperature

    double t_check = 0;
    int i_check = 1; 

  

public:
    void HeatCondSolver(double qdot_in,int Nx, double dt);
    void EvaluateTemp(vector<double>& f, double qdot0, double sdot0, double dt, const int Nx);
    void GenerateGeomGrid(const int Nx);
    void GenerateUniformGrid(int Nx);
    void ContractGridBoundaries(double sdot0, double dt, const int Nx);
    void Get_abcd_coeff(vector<double>& f, double qdot0, double sdot0, double dt, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d, const int Nx);
    void SolveTDMA(vector<double>& f, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d, const int Nx);
    bool CheckMelting(vector<double>& f);
    double GetMeltedLength(vector<double>& f);
    double GetRecessionRate(vector<double>& f, double dt);
    void init(double T0, int numPts);
    
};


#endif //SOLIDCONDRAY_H_INCLUDED

