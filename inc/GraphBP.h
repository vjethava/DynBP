#ifndef GRAPHBP_H_
#define GRAPHBP_H_
#include "mygraph.h"  
#include "PottsModel.h" 
using namespace std; 


class GraphBP: public FactorGraph {
	PottsModel* model; 
public: 

	GraphBP(PottsModel* _model) {
		model = _model; 
	}
	
	pair<EdgePtr, EdgePtr> addEdge(const NodePtr& n1, NodePtr const& n2, bool isParent, bool checkExisting=false);
	pair<EdgePtr, EdgePtr> findEdge(NodePtr nps, NodePtr npt);
	LD updateM(NodePtr np, NodePtr nc);
	LD updateN(NodePtr np, NodePtr nc);
	
};



#endif /*GRAPHBP_H_*/
