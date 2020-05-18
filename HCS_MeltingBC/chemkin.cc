#include "chemkin.h"
#include "solidCondRay.h"
#include "materialResponse.h"

void chemkin::init_chemParam(){

	// General parameters	
	Mw_s.resize(Ns);
	mdot_k.resize(Ns);	
	
	Mw_s = {M_o, M_co, M_co2, M_c, M_c2, M_c3, M_o2, M_n2, M_n, M_no}; // molar mass vector
	
	for (int i = 0; i < Ns; i++) {
		mdot_k[i] = 0;
	}

	// Sublimation parameters
	alpha_sub.resize(N_sub);
	Ps_sub.resize(N_sub);
	Qs_sub.resize(N_sub);	

	alpha_sub = { 0.14,0.26,0.030 };
	Ps_sub = { -85715,-98363,-93227 };
	Qs_sub = { 18.69,22.2,23.93 };
	

}



void chemkin::getOxidRates(double Tw) {
	
	double dodt, dcodt, dco2dt;	
	//double Fox = 1e23 / NA;
	

	C_O0 = Fox * 4 / sqrt(8 * kb * 1000 / pi / mo);

	//rate constants
	double k1 = (1. / (4. * B)) * sqrt((8. * kb * Tw) / (pi * mo));
	double k2 = (2. * pi * mo * kb * kb * Tw * Tw) * exp(-44277. / Tw) / (B * NA * h * h * h);
	double k3 = (1. / (4. * B)) * sqrt((8. * kb * Tw) / (pi * mo)) * 57.37 * exp(-4667. / Tw);
	double k4 = (1. / (4. * B)) * sqrt((8. * kb * Tw) / (pi * mo)) * 8.529e-6 * exp(6958. / Tw);
	double k5 = (1. / (4. * B)) * sqrt((8. * kb * Tw) / (pi * mo)) * 0.1203 * exp(2287. / Tw);

	S_ss = (k2 + k4 * C_O0) / (k1 * C_O0 + k2 + k4 * C_O0) * B; // sites [mol/m^2] steady state number of available sites	

	//production rates [mol]/m^2/s
	dodt = k2 * (B - S_ss) + Fox - k1 * C_O0 * S_ss - ((k3 + k4) * C_O0 * (B - S_ss)) - (k5 * C_O0 * S_ss);
	dcodt = (k3 * C_O0 * (B - S_ss)) + (k5 * C_O0 * S_ss);
	dco2dt = (k4 * C_O0 * (B - S_ss));

	mdot_k[0] = dodt * Mw_s[0]; // [mol/m^2sec]
	mdot_k[1] = dcodt * Mw_s[1]; // [mol/m^2sec]
	mdot_k[2] = dco2dt * Mw_s[2]; // [mol/m^2sec]	

}

void chemkin::getSublimRates(double Tw, vector<double>& ps_c)
{	
	
	// initialize production rates
	vector <double> pvs_c(3); 
	double dcdt, dc2dt, dc3dt;	
		
	// carbon species vapor pressure
	pvs_c[0] = 101300 * exp(Ps_sub[0] / Tw + Qs_sub[0]);
	pvs_c[1] = 101300 * exp(Ps_sub[1] / Tw + Qs_sub[1]);
	pvs_c[2] = 101300 * exp(Ps_sub[2] / Tw + Qs_sub[2]);

	//production rates of c, c2, c3 [kg/m^2/s]
	dcdt = alpha_sub[0] * (pvs_c[0] - ps_c[0]) * sqrt(Mw_s[3] / (2 * pi * R0 * Tw));
	dc2dt= alpha_sub[1] * (pvs_c[1] - ps_c[1]) * sqrt(Mw_s[4] / (2 * pi * R0 * Tw));
	dc3dt= alpha_sub[2] * (pvs_c[2] - ps_c[2]) * sqrt(Mw_s[5] / (2 * pi * R0 * Tw));

	mdot_k[3] = dcdt;
	mdot_k[4] = dc2dt;
	mdot_k[5] = dc3dt;	

}

void chemkin::getSpeciesFlux() {

	mdot_w = 0;

	for (int i = 0; i < Ns; i++) {

		mdot_w += mdot_k[i]; // [kg/m^2s] total flux of species
		
	}	


}