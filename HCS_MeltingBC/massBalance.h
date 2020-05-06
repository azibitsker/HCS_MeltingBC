#pragma once


#include "chemkin.h"
//#include "momentumBalance.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

class massBalance
{
public:
	
	vector<double> yk_w; // mass fraction of species at the wall
	vector<double> yk_gh; // mass fraction of species in the ghost cell 

public:

	void init_MassBal(int Ns,double yk_0);
	void solveMassBal(chemkin SR, double rho_w, vector<double> yk_f, vector<double> x);



};

