#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

class chemkin
{

public:

	//--------- Oxidation parameters--------------
	const int N_ox = 4;
	vector<double> surfcon_oxid;
	vector<double> mdot_oxid;
	double dt_sim = 1e-7;

	//constants needed	
	double NA = 6.0221409e+23; //[1/mol] Avogadro Number
	double mo = 2.6567e-26; //mass of oxygen atom in Kg
	double kb = 1.381e-23; //boltzman constant in [m^2kg/s^2K]
	double h = 6.626e-34; //planck's const in [m^2kg/s]
	double R0 = 8.3145; //[J/molK]
	double pi = 3.14159; //pi	

	// Molar masses of species [kg/mol]
	double M_o = 16.e-3;
	double M_co = 28.e-3;
	double M_co2 = 44.e-3;

	// Simulation conditions
	double B = 1e-5; //total site molar concentration [mol]/m^2	
	double S = 0.5e-5; // initial concentration of available sites [mol/m^2]
	double O_flux = 1e23 / NA; // [mol/m^2s] constant oxygen atoms flux	
	double C_O0 = (O_flux) * 4 / sqrt(8 * kb * 1000 / pi / mo); // [mol/m^2] constant oxygen concentration as the surface
	

	//---------- Sublimation parameters---------------------
	int N_sub = 3;
	vector<double> alpha_sub, Ps_sub, Qs_sub;
	vector<double> rhoi_sub, masscon_sub, ps_c;
	vector<double> mdot_sub;

	//Molar masses of carbon species [kg/mol]
	double M_c = 12.e-3;
	double M_c2 = 24.e-3;
	double M_c3 = 36.e-3;
	double rho_c = 1775; //[kg/m^3]


	// ---------Overall--------
	double mdot_w; // // overall flux of species out of the control volume
	vector<double> mdot_k; // vector of individual species mass flux
	vector<double> Mw_s; //molar weight of each species
	int Ns; // number of species
	double Sc = 0.5; //Schmidt number
	double mu_air = 47.88e-6; //[Pa*s] air dynamic viscosity at 1000 C
	double Ds; //species diffusion coefficient

public:

	void init_OxidParam();
	void init_SublimParam();

    void getOxidRates(double Tw);
	void getSublimRates(double Tw, vector<double>& ps_c);
	void getSpeciesFlux();


};

