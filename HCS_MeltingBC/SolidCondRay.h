#ifndef SOLIDCONDRAY_H_INCLUDED
#define SOLIDCONDRAY_H_INCLUDED

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

class SolidCond_Ray {

    double L0 = 0.04; // initial material length
    double F = 1.001; // cell size increase ratio
    double k = 0.402, rho = 280, Cp = 983.9; // material properties
    double Tm = 800; // [K] material melting temperature
    double Qstar = 1e5; // [J/kgK] heat of ablation
    double Eps_sdot = 1e-6;

public:
    double HeatCondSolver(double* f0, double* x0, double* x, double dt, double qdot0, double sdot0, const int Nx, int t); //
    void GenerateGrid(double* x0, const int Nx);
    void ContractGridBoundaries(double* x0, double* x, double sdot0, double dt, const int Nx, int ind_t);
    //void Get_abcd_coeff(double *x0,double *x,double *f0,double qdot0,double sdot0,double alpha, double dt,double *a,double *b,double *c,double *d,const int Nx,int ind_t);
    void Get_abcd_coeff(double* x0, double* x, double* f0, double qdot0, double sdot0, double alpha, double dt, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d, const int Nx, int ind_t);
    void SolveTDMA(double* f0, vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d, const int Nx);
    int GetMeltedCells(double* f0, const int Nx);
    double GetRecessionRate(double* x0, double dt, int Nmelt);
};
#endif // SOLIDCONDRAY_H_INCLUDED
