#ifndef MATERIALRESPONSE_H_INCLUDED
#define MATERIALRESPONSE_H_INCLUDED

#include "solidCondRay.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

class materialResponse {
  
public:
    
    vector<solidCond_Ray>  rays;     
    
public:    
    
    void Init_mr(int numRays, int numPtsPerRay_in, double T0_in);
    void SolveCondRays(vector<double> qdot_in, int numPtsPerRay_in, double dt_in);

};

#endif MATERIALRESPONSE_H_INCLUDED