#include "momentumBalance.h"
#include "massBalance.h"
#include "chemkin.h"


void momentumBalance::solveMomentumBalance(chemkin SR, massBalance MSB, double p_eta, double Tw)
{
	
	double Rs; // specific gas constant
	Mw_all = 0;
	
	for (int i = 0; i < SR.Ns; i++) {

		Mw_all += MSB.yk_w[i] * SR.Mw_s[i]; // compute overall molar weight of species at the wall

	}

	Rs = SR.R0 / Mw_all;

	U_w = (2 * Rs * Tw * SR.mdot_w * SR.mdot_w) / (sqrt(p_eta * p_eta - 4 * Rs * Tw * SR.mdot_w * SR.mdot_w));
	rho_w = (p_eta + sqrt(p_eta * p_eta - 4 * Rs * Tw * SR.mdot_w * SR.mdot_w)) / (2 * Rs * Tw);
	p_w= (p_eta + sqrt(p_eta * p_eta - 4 * Rs * Tw * SR.mdot_w * SR.mdot_w)) / 2;


}
