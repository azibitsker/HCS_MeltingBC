#include "massBalance.h"
#include "chemkin.h"
#include "momentumBalance.h"


void massBalance::init_MassBal(int Ns,double yk0_w){	
	
	yk_w.resize(Ns);
	yk_gh.resize(Ns);	

	for (int i = 0; i < Ns; i++) {

		yk_w[i] = yk0_w;

	}	

}


void massBalance::solveMassBal(chemkin SR, double rho_w, vector<double> x)
{
	double dx = (x[2] + x[1]) / 2 - (x[0] + x[1]) / 2; // distance between centers of ghost cell and first material cell

	// compute diffusion coefficient based on rho_w, mu_air@1000C and Sc=0.5 [ref. Keenan 1994]
	double Ds = SR.mu_air / (SR.Sc * rho_w);
	double y_sum = 0;

	for (int i = 0; i < SR.Ns; i++) {

		yk_gh[i] = 1 / (1 + SR.mdot_w * dx / (2 * rho_w * Ds)) * (SR.yk_f[i] * (1 - SR.mdot_w * dx / (2 * rho_w * Ds)) + (SR.mdot_k[i] * dx) / (rho_w * Ds));
		y_sum += yk_gh[i];		

		yk_w[i] = (SR.yk_f[i] + yk_gh[i]) / 2; // update species mass fractions at the boundary

	}

}