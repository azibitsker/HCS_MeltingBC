#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

class chemkin
{

public:

	// General parameters	
	int const Ns = 4; // number of species	
	double Sc = 0.5; //Schmidt number
	double Ds; //species diffusion coefficient
	double mu_air = 47.88e-6; //[Pa*s] air dynamic viscosity at 1000 C	
	double rho_c = 990; //[kg/m^3] camphor density

	double mdot_w; // // overall flux of species out of the control volume
	vector<double> mdot_k; // vector of individual species mass flux
	vector<double> Mw_s; //molar weight of each species
	vector<double> yk_f = {0., 78.0 / 100, 21. / 100, 1./100 };  // mass fraction of species in the fluid cell {C10H12O N2 O2 Ar}
	
	// Molar masses of species [kg/mol]	
	double M_c10h16o = 152.2334e-3, M_n2 = 28.e-3, M_o2 = 32.e-3,  M_ar = 39.948e-3;

	//Constants 		
	double NA = 6.0221409e+23; //[1/mol] Avogadro Number	
	double kb = 1.381e-23; //boltzman constant in [m^2kg/s^2K]	
	double R0 = 8.3145; //[J/molK]
	double pi = 3.14159; //pi	

	// Simulation conditions
	

	//---------- Sublimation parameters---------------------
	int N_sub = 1;
	double alpha_sub, A_sub,B_sub, C_sub;
	double ps_c;	

public:

	void init_chemParam();    
	void getSublimRates(double Tw, double ps_c);
	void getSpeciesFlux();


};

