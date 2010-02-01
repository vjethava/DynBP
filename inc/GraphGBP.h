/**
 * This class re-implements GBP message passing algorithm as applied to time 
 */
#ifndef GRAPHGBP_H_
#define GRAPHGBP_H_
#include "mygraph.h" 
#include "RegionPPM.h" 
#include "PottsGraph.h" 
class GraphGBP : public PottsGraphPPM {
protected:
	/// returns the prod_{parent} parent-> region contribution {DENOMINATOR}
	LD lambdaAlphaToParent(RegionPPM* alpha, StatePtr currentState, StatePtr nextState);
	/// returns the prod_{child} child -> region contribution {NUMERATOR}    	
	LD lambdaChildToAlpha(RegionPPM* alpha, StatePtr currentState, StatePtr nextState);
	/// updates the \f$m_\alpha\f$ for given region \f$\alpha\f$   	
	LD getMalphaUpdate(NodePtr vn);
	/// computes the joint belief for given node \f$\alpha\f$
	LD getBeliefUpdate(NodePtr vn, StatePtr cs, StatePtr ns);
	/// Computes the parent->child message update 
	LD getMpcUpdate(NodePtr vp, NodePtr vc, StatePtr jointSt);
public:
	/** initializes the potts graph structure (adds nodes/edges)
	 * and updates the pHat vector for each node 
	 */ 
	void construct_graph(bool isNodeInitProbRandom, LD theta, LD delta_t); 
	/// initializes all the messages to 1.0   
	void step_init(); 
	/// updates the joint belief vectors and stores in RegionPPM::bJointRev
	void step_updateJointB(); 
	/// updates all child->parent messages 
	void step_updateMpc(); 
	
	/// the constructs that allows construction of PottsModel based graphs using PottsGraphPPM code
	GraphGBP(PottsModel* _M, int _rx, int _ry, int _cx, int _cy) : 
      PottsGraphPPM::PottsGraphPPM(_M, _rx, _ry, _cx, _cy) {
    }
	
};
#endif /*GRAPHGBP_H_*/

/* REVISION LOG
 * 
 * 07 July 08 - started recoding GBP for MVI 
 *
 *
 */
