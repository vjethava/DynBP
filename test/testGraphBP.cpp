#include "GraphBP.h"
#include "node.h" 
#include <string> 
#include <vector>
using namespace std; 


// tests a network of the form:
// f1->x1->f(x3,x4|x1,x2)<->x3
// f2->x2->				 <->x4

void updateBlf(Belief* b, FactorGraph* graph, int Q, int numVar, bool isRandom, int seed=0) {
	int NS = 1;
	// global vector of states
	vector<int> vals;
	for(int i=0; i < numVar; i++) { NS*= Q; vals.push_back(i); }
	srand(seed);
	LD pEq  = ((LD) (1.0))/((LD) NS);
	
	for(int i=0; i < NS; i++) {
		vector<INT> vvec = graph->getIthStVec(i, numVar, vals); 
		StatePtr sptr = graph->getStatePtr(vvec); 	
		if(isRandom == false) {
			b->setPr(sptr, pEq); 
		} else {
			int r = rand()%NS; 
			b->setPr(sptr, r);
		}
	}
	if(isRandom) {
		b->normalize(); 
	}
}

int testConnectivity() {
	GraphBP graph(0);
	VI v0, v1; v0.push_back(0); v1.push_back(1); 
	StatePtr s0 = new State(v0); 
	StatePtr s1 = new State(v1); 
	Node *n1 = new Node("X1", 2, new Belief(), VARIABLE);
	Node *n2 = new Node("X2", 2, new Belief(), VARIABLE);
	n1->getBeliefPtr()->setPr(s0, 0.5); 
	n1->getBeliefPtr()->setPr(s1, 0.5); 
	n2->getBeliefPtr()->setPr(s0, 0.5); 
	n2->getBeliefPtr()->setPr(s1, 0.5); 
	Node* f1 = new Node("F_X1", 2, new Belief(), FACTOR); 
	Node*f12 = new Node("F_X1_X2", 2, new Belief(), FACTOR);
	f1->getBeliefPtr()->setPr(s0, 0.9);
	f1->getBeliefPtr()->setPr(s1, 0.1);
	vector<int> v[4]; 
	vector<int> vals; vals.push_back(0); vals.push_back(1); 
	for(int i=0; i < 4; i++) {
		v[i] = graph.getIthStVec(i, 2, vals);
	}
	
	graph.addNode(n1); 
	graph.addNode(n2); 
	graph.addEdge(f1, n1, true); 
	graph.addEdge(f12, n1, false); 
	graph.addEdge(f12, n2, false); 
	graph.addEdge(f12, n2, true); 
	graph.bp_iter(); 
	cout<<graph; 
	 
	return 0; 
}

int main(){ 
	testConnectivity();
	return 0;
} 
