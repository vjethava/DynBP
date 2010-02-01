
#include <iostream>
#include <ctime>
#include "PottsGraph.h"
#include "VideoGraph.h"
using namespace std; 

int main() {
    int N1=5, N2=5, Q=2, seed=0; 
    //cout<<" Enter N1: " cin>>N1;
    //cout<<" Enter N1: " cin>>N1;
    PottsModel* mp = new PottsModel(N1, N2, Q, 0.1, 1.0, seed);
    mp->genFactors();
    PottsGraphPPM* graph = new PottsGraphPPM(mp, 2, 2, 1, 1); 
    graph->genSubRegions(); 
    // MyGraph* G = (MyGraph*) graph;
    //  cout<<*G; 
    
    // Revision 0.0 BEGIN CUT
    // graph->PPM_algo(0.5, 1.0);
    // Revision 0.0 END CUT
    
    // Revision 0.1 BEGIN CUT
    string modelFileName="ppm.mod";
    string dataFileName="ppm.dat";
    string amplFileName="ppm.ampl";
    graph->getAMPL(modelFileName, dataFileName, amplFileName);  
    // Revision 0.1 END CUT 
    
    return 0; 
}

/* REVISION LOG
 * ============
 *  
 * 0.1 July 17, 2008: Commented out PPM_algo testing, use for testing AMPL code generation
 * 
 */