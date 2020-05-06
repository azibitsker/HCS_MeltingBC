#include "chemkin.h"
#include "solidCondRay.h"
#include "materialResponse.h"

void chemkin::init_OxidParam()
{

	surfcon_oxid.resize(N_ox);
	mdot_oxid.resize(N_ox-1);
	Ns = N_ox - 1 + N_sub;
	Mw_s.resize(Ns);
	mdot_k.resize(Ns);

	//surfcon_oxid[0] = 0;// (475. / 8.314 * M_o); // O
	//surfcon_oxid[1] = 0;// 1.e-5 * M_co; // CO
	//surfcon_oxid[2] = 0;// 1.e-5 * M_co2; // CO2
	//surfcon_oxid[3] = S; // sites

	Mw_s = {M_o, M_co, M_co2, M_c, M_c2, M_c3}; // molar mass vector

}

void chemkin::init_SublimParam()
{
	alpha_sub.resize(N_sub);
	Ps_sub.resize(N_sub);
	Qs_sub.resize(N_sub);
	rhoi_sub.resize(N_sub);
	masscon_sub.resize(N_sub);	
	mdot_sub.resize(N_sub);

	alpha_sub = { 0.14,0.26,0.030 };
	Ps_sub = { -85715,-98363,-93227 };
	Qs_sub = { 18.69,22.2,23.93 };	

	for (int i = 0; i < N_sub; i++) {
		rhoi_sub[i] = 0;		
	}

}

void chemkin::getOxidRates(double Tw) {
	/*
	 * find values of O, CO, CO2 and available adsorbtion sites on surface after
	 * chemical kinetic model converges
	 *
	 * input:
	 * double &rhoi - mass densities of O, CO and CO2 on surface (for 5species is rhoi[4] in KATS!)
	 * double &S - double containing site concentration
	 * double &T - surface temperature
	 * double &timeStep - time step size of 1 simulation step dt[gid].cell(c)
	 *
	 * output:
	 * vector<double> - vector containing mass density values for O
	 *                  and molar concentration of CO, CO2, sites
	 */

	 // organize vector with molar concentrations	
	double timeStep = dt_sim;	

	//rate constants
	double k1 = (1. / (4. * B)) * sqrt((8. * kb * Tw) / (pi * mo));
	double k2 = (2. * pi * mo * kb * kb * Tw * Tw) * exp(-44277. / Tw) / (B * NA * h * h * h);
	double k3 = (1. / (4. * B)) * sqrt((8. * kb * Tw) / (pi * mo)) * 57.37 * exp(-4667. / Tw);
	double k4 = (1. / (4. * B)) * sqrt((8. * kb * Tw) / (pi * mo)) * 8.529e-6 * exp(6958. / Tw);
	double k5 = (1. / (4. * B)) * sqrt((8. * kb * Tw) / (pi * mo)) * 0.1203 * exp(2287. / Tw);

	surfcon_oxid[3] = (k2 + k4 * C_O0) / (k1 * C_O0 + k2 + k4 * C_O0) * B; // sites [mol/m^2] steady state number of available sites

	// total steps for netwon solver
	double dt_loc = 1e-10;
	//int i_tot = timeStep / dt_loc;
	
	
	// initialize production rates
	double dodt = 0;
	double dsdt = 0;
	double dcodt = 0;
	double dco2dt = 0;
	double surfcov = 1;//(B - surfcon_oxid[3]) / B;
	double surfcov_conv = surfcov;

	int j=1;
	//start explicit netwon method: molcon_(n+1)=molcon_(n)+(dt*d(molcon)/dt)
	//for (int i=0;i<i_tot;i++) {
	
		//
		

		//production rates [mol]/m^2/s
		dodt = k2 * (B - surfcon_oxid[3]) + O_flux - k1 * C_O0 * surfcon_oxid[3] - ((k3 + k4) * C_O0 * (B - surfcon_oxid[3])) - (k5 * C_O0 * surfcon_oxid[3]);
		//dsdt = (-k1 * C_O0 * surfcon_oxid[3]) + (k2 * (B - surfcon_oxid[3])) + (k4 * C_O0 * (B - surfcon_oxid[3]));
		dcodt = (k3 * C_O0 * (B - surfcon_oxid[3])) + (k5 * C_O0 * surfcon_oxid[3]);
		dco2dt = (k4 * C_O0 * (B - surfcon_oxid[3]));

		//update values
		//surfcon_oxid[0] += dt_loc * dodt;
		//surfcon_oxid[1] += dt_loc * dcodt;
		//surfcon_oxid[2] += dt_loc * dco2dt;
		//surfcon_oxid[3] += dt_loc * dsdt;

		//surf_cov_p = surfcon_oxid[0] * surfcon_oxid[0] + surfcon_oxid[1] * surfcon_oxid[1] + surfcon_oxid[2] * surfcon_oxid[2] + surfcon_oxid[3] * surfcon_oxid[3];

		//surfcov_conv = abs(surfcov - surfcon_oxid[3]) / surfcov;
		//surfcov = surfcon_oxid[3];
		//surfcov_conv = abs(surfcov - (B - surfcon_oxid[3]) / B);
		//surfcov = (B - surfcon_oxid[3]) / B;
		
		//if (surfcov_conv < 1.e-5) {
		//	cout << "iter: " << i << " out of " << i_tot << "\n";
		//	break;
		//}
		//j++;
	//}

	//cout << "j = " << j - 1 << endl;
	//cout << "dodt: " << dodt * Mw_s[0] <<" vs "<< surfcon_oxid[0]* Mw_s [0]/ timeStep << "[kg/m^2sec]" << endl;
	//cout << "dcodt: " << dcodt * Mw_s[1] << " vs " << surfcon_oxid[1] * Mw_s[1] / timeStep <<" [kg/m^2sec]" << endl;
	//cout << "dco2dt: " << dco2dt * Mw_s[2] << " vs " << surfcon_oxid[2] * Mw_s[2] / timeStep <<" [kg/m^2sec]"<< endl;
	//cout << "dsdt: " << dsdt << " vs " << (surfcon_oxid[3]-S) / timeStep <<" [mol/m^2sec]"<< endl;

	//surfcon_oxid[0] *= Mw_s[0]; // convert O molar concentration to mass density [kg/m^2]
	//surfcon_oxid[1] *= Mw_s[1]; // convert CO molar concentration to mass density [kg/m^2]
	//surfcon_oxid[2] *= Mw_s[2]; // convert CO2 molar concentration to mass density [kg/m^2]

	mdot_oxid[0] = dodt * Mw_s[0]; // [kg/m^2sec]
	mdot_oxid[1] = dcodt * Mw_s[1]; // [kg/m^2sec]
	mdot_oxid[2] = dco2dt * Mw_s[2]; // [kg/m^2sec]	

}

void chemkin::getSublimRates(double Tw,vector<double>& ps_c)
{
	// organize vector with molar concentrations	
	double timeStep = dt_sim;

	masscon_sub[0] = rhoi_sub[0];   // C density (kg/m^3) 
	masscon_sub[1] = rhoi_sub[1];  // C2 density  (kg/m^3)
	masscon_sub[2] = rhoi_sub[2];  // C3 density  (kg/m^3) 	

	// total steps for netwon solver
	double dt_loc = 1e-10;
	//int i_tot = timeStep / dt_loc;	
	

	// initialize production rates
	vector <double> pvs_c(3); 
	double dcdt, dc2dt, dc3dt;

	//start explicit netwon method: molcon_(n+1)=molcon_(n)+(dt*d(molcon)/dt)
	//for (int i = 0; i < i_tot; i++) {
		//
		
		// carbon species vapor pressure
		pvs_c[0] = 101300 * exp(Ps_sub[0] / Tw + Qs_sub[0]);
		pvs_c[1] = 101300 * exp(Ps_sub[1] / Tw + Qs_sub[1]);
		pvs_c[2] = 101300 * exp(Ps_sub[2] / Tw + Qs_sub[2]);

		//production rates of c, c2, c3 [kg/m^2/s]
		dcdt = alpha_sub[0] * (pvs_c[0] - ps_c[0]) * sqrt(Mw_s[3] / (2 * pi * R0 * Tw));
		dc2dt= alpha_sub[1] * (pvs_c[1] - ps_c[1]) * sqrt(Mw_s[4] / (2 * pi * R0 * Tw));
		dc3dt= alpha_sub[2] * (pvs_c[2] - ps_c[2]) * sqrt(Mw_s[5] / (2 * pi * R0 * Tw));

		//update values
		//masscon_sub[0] += dt_loc * dcdt;
		//masscon_sub[1] += dt_loc * dc2dt;
		//masscon_sub[2] += dt_loc * dc3dt;			

	//}
		mdot_sub[0] = dcdt;
		mdot_sub[1] = dc2dt;
		mdot_sub[2] = dc3dt;

	//for (int i = 0; i < N_sub; i++) {

		//mdot_sub[i] = (masscon_sub[i] - rhoi_sub[i]) / timeStep; // [kg/m^2sec]
	//}


}

void chemkin::getSpeciesFlux() {

	mdot_w = 0;

	for (int i = 0; i < N_ox-1; i++) {

		mdot_w += mdot_oxid[i]; // total flux of species
		mdot_k[i] = mdot_oxid[i]; // vector of each species flux
	}

	for (int i = 0; i < N_sub; i++) {

		mdot_w += mdot_sub[i];
		mdot_k[i+ N_ox - 1] = mdot_sub[i];
	}


}