#include "SolidCondRay.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
using namespace std;

class materialResponse {
  
  vector<SolidCond_Ray> :: rays;
  int                   :: numRays; 

  
public:  

  materialReponse(int numRays_in,int numPtsPerRay_in,double T0_in)
  {
    numRays=numRays_in;
    rays.resize(numRays);
    for (int r=0;r<numRays;++r)
      {
        rays[r].init(T0_in,numPtsPerRay_in);
      }    
  }

  void materialResponse::updateRays()
  {
    for (int r=0;r<numRays;++r)
      {
        rays[r].HeatCondSolver();
        
      }
  }
  
};

