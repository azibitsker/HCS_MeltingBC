#include "momentumBalance.h"
#include "massBalance.h"
#include "chemkin.h"
#include "solidCondRay.h"


void momentumBalance::solveMomentumBalance(chemkin SR, massBalance MasB, double p_eta, double Tw)
{
	solidCond_Ray sC;

	double Rs; // specific gas constant
	double Mw_ps; // average molar mass of produced species
	Mw_ps = 0; // average species molar weight
	
	for (int i = 0; i < SR.Mw_s.size(); i++) {

		Mw_ps += MasB.yk_w[i] * SR.Mw_s[i]; // compute overall molar weight of blowing species at the wall

	}

	Rs = SR.R0 / Mw_ps;

	// Compute species velocity, density and pressure at the wall
	double denom = ((p_eta * p_eta - 4 * Rs * Tw * SR.mdot_w * SR.mdot_w));

	if (denom <= 0) { cout << "ts: "<<sC.timeStep<<" The term in moment balance (p_eta * p_eta - 4 * Rs * Tw * SR.mdot_w * SR.mdot_w) is negative!" << endl; exit(1); }

	U_w = (2 * Rs * Tw * SR.mdot_w * SR.mdot_w) / (sqrt(p_eta * p_eta - 4 * Rs * Tw * SR.mdot_w * SR.mdot_w));
	rho_w = (p_eta + sqrt(p_eta * p_eta - 4 * Rs * Tw * SR.mdot_w * SR.mdot_w)) / (2 * Rs * Tw);
	p_w= (p_eta + sqrt(p_eta * p_eta - 4 * Rs * Tw * SR.mdot_w * SR.mdot_w)) / 2;

}
