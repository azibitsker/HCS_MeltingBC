#include "momentumBalance.h"
#include "massBalance.h"
#include "chemkin.h"


void momentumBalance::solveMomentumBalance(chemkin SR, massBalance MasB, double p_eta, double Tw)
{
	
	double Rs; // specific gas constant
	double Mw_ps; // average molar mass of produced species
	Mw_ps = 0; // average species molar weight
	
	for (int i = 0; i < SR.Ns-4; i++) {

		Mw_ps += MasB.yk_w[i] * SR.Mw_s[i]; // compute overall molar weight of produced species at the wall

	}

	Rs = SR.R0 / Mw_ps;

	// Compute species velocity, density and pressure at the wall
	U_w = (2 * Rs * Tw * SR.mdot_w * SR.mdot_w) / (sqrt(p_eta * p_eta - 4 * Rs * Tw * SR.mdot_w * SR.mdot_w));
	rho_w = (p_eta + sqrt(p_eta * p_eta - 4 * Rs * Tw * SR.mdot_w * SR.mdot_w)) / (2 * Rs * Tw);
	p_w= (p_eta + sqrt(p_eta * p_eta - 4 * Rs * Tw * SR.mdot_w * SR.mdot_w)) / 2;

}
