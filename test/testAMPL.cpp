/**
 * This file generates the AMPL files for testing purpose. 
 */
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <pthread.h>
#include "main.h" 
#include "GraphPPM.h"
#include "PottsModel.h"  
#include "PottsGraph.h" 
using namespace std;

int main() {
   int N1=3, N2=3, Q=2, seed=0; 
    PottsModel* mp = new PottsModel(N1, N2, Q, 0.1, 1.0, seed);
    mp->genFactors();
    PottsGraphPPM* graph = new PottsGraphPPM(mp, 2, 2, 1, 1); 
    graph->genSubRegions(); 
    // Revision 0.1 BEGIN CUT
    string modelFileName="ppm.mod";
    string dataFileName="ppm.dat";
    string amplFileName="ppm.ampl";
    graph->getAMPL(modelFileName, dataFileName, amplFileName);  
    // Revision 0.1 END CUT 
  return 0; 
}




