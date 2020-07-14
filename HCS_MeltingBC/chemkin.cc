#include "chemkin.h"
#include "solidCondRay.h"
#include "materialResponse.h"

void chemkin::init_chemParam(){

	// General parameters	
	Mw_s.resize(Ns);
	mdot_k.resize(Ns);	
	
	Mw_s = {M_c10h16o, M_n2, M_o2, M_ar}; // molar mass vector
	
	for (int i = 0; i < Ns; i++) {
		mdot_k[i] = 0;
	}

	// Sublimation parameters
	alpha_sub = 0.18;
	A_sub = 8.52;
	B_sub = 2714.91;
	C_sub = 277.67;

}




void chemkin::getSublimRates(double Tw, double ps_c)
{	
	
	// initialize production rates
	double pvs_c; 	
		
	// camphor species vapor pressure
	pvs_c = 133.322 * pow(10, A_sub - B_sub / (Tw - 273.15 + C_sub)); 
	

	//production rates camphor gas [kg/m^2/s]

	mdot_k[0] = alpha_sub * sqrt(Mw_s[0] / (2 * pi * R0 * Tw))* (pvs_c - ps_c);
	

}

void chemkin::getSpeciesFlux() {

	mdot_w = 0;

	for (int i = 0; i < Ns; i++) {

		mdot_w += mdot_k[i]; // [kg/m^2s] total flux of species
		
	}	


}