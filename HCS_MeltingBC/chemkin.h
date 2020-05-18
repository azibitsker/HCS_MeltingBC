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
	int const Ns = 10; // number of species	
	double Sc = 0.5; //Schmidt number
	double Ds; //species diffusion coefficient
	double mu_air = 47.88e-6; //[Pa*s] air dynamic viscosity at 1000 C	
	double rho_c = 1775; //[kg/m^3] graphite density

	double mdot_w; // // overall flux of species out of the control volume
	vector<double> mdot_k; // vector of individual species mass flux
	vector<double> Mw_s; //molar weight of each species
	vector<double> yk_f = { 4.34e-2, 0, 0 , 0, 0, 0, 1.6e-1, 7.5e-1, 2.6e-3 ,4.4e-2 };  // mass fraction of species in the fluid cell {O,CO,CO2,C,C2,C3,O2,N2,N,NO}
	
	// Molar masses of species [kg/mol]
	double M_o = 16.e-3, M_co = 28.e-3, M_co2 = 44.e-3;
	double M_c = 12.e-3, M_c2 = 24.e-3, M_c3 = 36.e-3;
	double M_o2 = 32.e-3, M_n2 = 28.e-3, M_n = 14.e-3, M_no = 30.e-3;	

	//--------- Oxidation parameters--------------		
	//constants 	
	double NA = 6.0221409e+23; //[1/mol] Avogadro Number
	double mo = 2.6567e-26; //mass of oxygen atom in Kg
	double kb = 1.381e-23; //boltzman constant in [m^2kg/s^2K]
	double h = 6.626e-34; //planck's const in [m^2kg/s]
	double R0 = 8.3145; //[J/molK]
	double pi = 3.14159; //pi	

	// Simulation conditions
	double B = 1e-5; //total site molar concentration [mol]/m^2	
	double S0 = 0.5e-5; // initial concentration of available sites [mol/m^2]
	double S_ss; // steady state surface concentration of availables sites	 
	double C_O0; // [mol/m^2] constant oxygen concentration as the surface
	double Fox = 1e21/NA;

	//---------- Sublimation parameters---------------------
	int N_sub = 3;
	vector<double> alpha_sub, Ps_sub, Qs_sub;
	vector<double> ps_c;	

public:

	void init_chemParam();
    void getOxidRates(double Tw);
	void getSublimRates(double Tw, vector<double>& ps_c);
	void getSpeciesFlux();


};

