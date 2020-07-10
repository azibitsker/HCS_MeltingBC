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
    
    // Camphor properties
    double k_s = 0.2, rho_s = 990, Cp_s = 1.781e3;  
    double Tmelt = 453.3; 
    double Eps_T = 1e-3;
    double alpha; // thermal diffusivity
    double sdot_out;   
    double eta = 0.5;
    vector<double> f0, x0, x, mdot_k_out, yk_w_out;
    double mdot_w_out; 


    // Radiation parameters
    double Epsilon = 0.84;
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
    void init(double T0, int numPts);
    
};


#endif //SOLIDCONDRAY_H_INCLUDED

