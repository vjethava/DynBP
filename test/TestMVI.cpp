#include "PottsGraph.h"
#include "VR.h"
#include "MVI.h" 
using namespace std; 

int main(int argc, char* argv[]) {
	string fname; 
	if (argc != 2)
    {
        fname = "/home/vjethava/bw80x50.avi"; 
    } else 
        fname = string(argv[1]);  
	// string fname = "/home/vjethava/eclipse/GBP/data/jumping_girl_input.avi";
	
	MVI mvi(fname, 2, 0.99, 0.70); 
	mvi.init(); 
	//for(int i=0; i < 30; i++) { mvi.run(true, 0); } 
	//for(int i=30; i < 110; i++) { mvi.run(true, 70); }
	for(int i=0; i < 200; i++) { mvi.run(true, 150); } 
	return 0; 
}
