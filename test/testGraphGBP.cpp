#include <iostream>
#include <fstream> 
#include <string>
#include <cstdio> 
#include <sstream> 
#include "GraphGBP.h" 
using namespace std;
/// this results in the readme file 
void readme( string str=""); 

int main() {
	
	// model variables
	int seed = 0;  
	int N1 = 3; 
	int N2 = 3; 
	int Q = 2; 
	
	LD J = 0.1; 
	LD H = 0.1;
	PottsModel* model = new PottsModel(N1, N2, Q, J, H, seed);  
	model->genFactors();
	
	// initialize the graph 
	int rx = 2; 
	int ry = 2; 
	int cx = 1; 
	int cy = 1; 
	LD theta = 0.1; 
	LD delta_t = 1.0; 
	bool isNodeInitProbRandom = true; 
	GraphGBP* graph = new GraphGBP(model, rx, ry, cx, cy);
 	graph->construct_graph(isNodeInitProbRandom, theta, delta_t); 
	
	// view the graph; 
	cout<<*graph; 
}


void readme(string str) {
	FILE* file = fopen("README.txt", "w"); 
	string message = "The GBP experiment tests GBP with  Potts Model.\n" ;  
	fprintf(file, "%s\n", message.c_str()) ; 
	if(str != "") fprintf(file, "%\n", str.c_str()); 
	fclose(file); 
}
