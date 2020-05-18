#pragma once

#include "massBalance.h"
#include "chemkin.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

class momentumBalance
{
public:

	// conditions at the wall
	double U_w, rho_w, p_w;	


public:

	void solveMomentumBalance(chemkin SR, massBalance MSB,double p_eta,double Tw);



};

