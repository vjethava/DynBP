/**
 * This class re-implements the GBP modification for PPM
 */
#include "GraphGBP.h" 

/// constructs the potts graph and initializes pHat for each node 
void GraphGBP::construct_graph(bool isNodeProbRandom, LD theta, LD delta_t) {
	PottsGraphPPM::genSubRegions();
	PottsGraphPPM::initAllNodesOrigP(isNodeProbRandom);
	PottsGraphPPM::PPM_step2(theta, delta_t); 	
}

/// initializes all the messages 
void GraphGBP::step_init() {
	PottsGraphPPM::PPM_step3();
	// initialize all region normalization beliefs
	foreach(nodeIter, adj_map) {
		RegionPPM* node = (RegionPPM*) nodeIter->first; 
		BidirEdgeList& edgeList = nodeIter->second;
		node->m_alpha = 1.0; 
	}
}

/**
 * Computes the message contributions from node from its parents
 * \arg vn: node ptr \f$\alpha\f$
 * \arg cs: node current state \f$ \vec{x}_\alpha^t\f$
 * \arg ns: node next state \f$ \vec{x}_\alpha^{t+\delta t} \f$
 * \returns \f$\prod_{\{\gamma:\alpha\subset\gamma\} } m_{\alpha\to\gamma}(\vec{x}^{t,\delta t}_\alpha)\f$
 */
LD GraphGBP::lambdaAlphaToParent(RegionPPM* vn, StatePtr cs, StatePtr ns) {	
	StatePtr jointSt = this->getJointState(cs, ns); 
	ADJ_MAP_ITER nodeIter = adj_map.find(vn);
	ASSERT(nodeIter != adj_map.end());  
	vector<MyEdge> revEdges = nodeIter->second.second;
	LD result = 1.0;  
	foreach(edgeIter, revEdges) { // each child->parent contribution
		MyEdge& edge = *edgeIter;  
		LD message = edge.getMsg(jointSt); 
		if(message == 0.0) { message = EPSILON_ZERO; } 
		result = result*message; 	
	}
	return result; 
} 

/**
 * Computes all children contribution for given node
 * \arg vn: node ptr \f$\alpha\f$
 * \arg cs: current node state \f$ \vec{x}_\alpha^t\f$
 * \arg ns: next node state \f$ \vec{x}_\alpha^{t+\delta t} \f$
 * \returns \f$\prod_{\{\beta:\beta\subset\alpha\} } m_{\beta\to\alpha}(\vec{x}^{t,\delta t}_\beta)\f$
 */
LD GraphGBP::lambdaChildToAlpha(RegionPPM* vn, StatePtr cs, State* ns) {
	ADJ_MAP_ITER pIter = adj_map.find(vn); 
	ASSERT(pIter != adj_map.end());
	LD result = 1.0;  
	vector<MyEdge>& forwardEdges = pIter->second.first; 
	foreach(edgeIter, forwardEdges) { 
		NodePtr& child = edgeIter->ends[1]; 
		pair<EdgePtr, EdgePtr> edgePair = findEdge(vn, child);
		MyEdge& revEdge = *(edgePair.second);
		StatePtr childCurrSt = getChildStatePtr(vn, child, cs);
		StatePtr childNextSt = getChildStatePtr(vn, child, ns);
		StatePtr jointSt = getJointState(childCurrSt, childNextSt); 
		LD message = revEdge.getMsg(jointSt);
		if(message == 0.0) message = EPSILON_ZERO;  
		result *= message;   
	}
	return result; 
}

/**
 * Computes the belief for region \f$\alpha\f$ 
 * \arg vn region \f$\alpha\f$ 
 * \arg cs current state \f$\vec{x}_\alpha^{t}\f$
 * \arg ns next state \f$\vec{x}_\alpha^{t}\f$
 * \returns \f$b_\alpha(x_\alpha^t)\times\frac{\lambda^{t+\delta t | t}(x_\alpha^{t,\delta t})}{e\times m_\alpha^{(i)} }\times \Big\{\frac{\lambda_{child\to\alpha}(x_\alpha^{t,\delta t})}{\lambda_{\alpha\to parent}(x_\alpha^{t,\delta t})}\Big\}^{1/c_\alpha}\f$ 
 */
 LD GraphGBP::getBeliefUpdate(NodePtr vn, StatePtr cs, StatePtr ns) {
 	RegionPPM* node = (RegionPPM*) node; 
 	LD bOrig = node->getBeliefPtr()->getPr(cs); 
 	LD lambda_cond = node->pHat->getPr(cs, ns);
 	LD lambda_par = lambdaAlphaToParent(node, cs, ns);
 	LD lambda_child = lambdaChildToAlpha(node, cs, ns); 
 	int c_alpha = node->getCR();
 	ASSERT((lambda_par > 0.0) && (lambda_child > 0.0));  	  
 	LD factor = lambda_child/lambda_par;
 	LD multiplier = (LD) pow(factor, ((1.0)/((double) c_alpha) )); 
 	LD bJoint = bOrig * lambda_cond * multiplier;   
 	return bJoint; 
 } 
 
 /**
  * NOTE: Assumes RegionPPM::bJointRev as \f$b_\alpha (x^{t+\delta t}, x^t)\f$ 
  * \arg vp parent node \f$\alpha\f$ 
  * \arg vc child node \f$\beta\f$
  * \arg jointSt child joint state \f$x_\beta^{t,\delta t}\f$ 
  * \returns \f[ m_{\beta\to\alpha}^{(i+1)}(x_\beta^{t,\delta t}) = 
  * m_{\beta\to\alpha}^{(i)}(x_\beta^{t,\delta t})\times \Big\{\frac{b_\beta^{(i)}(x_\beta^{t, \delta t})}{\sum_{x_\alpha\backslash\beta^{t,\delta t}} b_\alpha^{(i)}(x_\alpha^{t, \delta t})}\Big\}^{\frac{c_\alpha c_\beta}{c_\alpha + c_\beta}}\f]
  * 
  */
 LD GraphGBP::getMpcUpdate(NodePtr vp, NodePtr vc, StatePtr jointSt) {
 	RegionPPM* parent = (RegionPPM*) vp;
 	RegionPPM* child = (RegionPPM*) vc;
 	pair<StatePtr, StatePtr> childStPair = getChildFromMsgState(jointSt);
 	StatePtr& childPastSt = childStPair.first;
 	StatePtr& childNextSt = childStPair.second;  
 	// we iterate over the states present in vp->pHat allowing reduced state space search
 	int c_alpha = parent->getCR();
 	int c_beta = child->getCR(); 
 	// returns pointer to the joint beliefs for parent 
 	map<StatePtr, Belief, StatePtrCmp>* cbMpPtr = parent->bJointRev->getCbMpPtr();
 	LD sumParMarginalizedBlfs = 0.0;
 	foreach(mpIter, *cbMpPtr) { // iterate over the reversed joint beliefs
 		StatePtr parNextSt = mpIter->first;
 		StatePtr parNextStReducedToChild = getChildStatePtr(parent, child, parNextSt);
		if(State::SSD(*parNextStReducedToChild, *childNextSt) == 0.0) {	 
 			map<StatePtr, LD, StatePtrCmp> oldStateMp = mpIter->second.getMp(); 
 			foreach(innerIter, oldStateMp) {
 				StatePtr parPastSt = innerIter->first;
 				StatePtr parPastStReducedToChild = getChildStatePtr(parent, child, parPastSt);
 				if(State::SSD(*parPastStReducedToChild, *childPastSt) == 0.0) {
 					// add belief to margSum
 					sumParMarginalizedBlfs += innerIter->second; 	
 				}
 			}
 		}	
 	}
 	ASSERT(sumParMarginalizedBlfs > 0.0); 
 	// note reversed order for accessing bChild
 	LD bChild = child->bJointRev->getPr(childNextSt, childPastSt);
 	// old message value
 	LD mOld = getMsg(parent, child, jointSt, false);
 	LD factor = bChild/sumParMarginalizedBlfs;
 	LD powMult;  
 	if( (c_alpha + c_beta) != 0 ) {
 		powMult = ((double) (c_alpha*c_beta) )/ ((double) (c_alpha+c_beta));  
 	} else {
 		// note it goes to zero or -infinity limiting case depending
 		// on whether factor>1.0 or factor<1.0
 		powMult = -1.0; 
 	}	
 	LD mult = pow(factor, powMult); 
 	LD result = mOld*mult; 
 	return result; 
 }
 
 
 
 
/**
 * NOTE: uses bJointRev as update scheme \f$b_\alpha (x^{t+\delta t}, x^{t})\f$
 * \arg vn the region \f$\alpha\f$ to be updated
 */
LD GraphGBP::getMalphaUpdate(NodePtr vn) {
	RegionPPM* node = (RegionPPM*) vn;
	// at present un-implemented as just a normalizing factor - can be ignored  
	return 1.0; 
}

void GraphGBP::step_updateJointB() {
	// update joint beliefs for all nodes
	foreach(nodeIter, adj_map) {
		LD sumB = 0.0; // new joint beliefs normalization constant
		RegionPPM* node = (RegionPPM*) nodeIter->first; 
		BidirEdgeList& edgeList = nodeIter->second; 
		node->bJointRev->clear();
		// look at permissible changes as outlined in pHat
		map<StatePtr, Belief, StatePtrCmp>* cbMpPtr = node->pHat->getCbMpPtr();
		foreach(mpIter, *cbMpPtr) {
			StatePtr nodePastSt = (StatePtr) mpIter->first;
			// next state possibilities 
			BeliefPtr nStBlf = node->pHat->getBeliefForState(nodePastSt);
			map<StatePtr, LD, StatePtrCmp> newStateMp = nStBlf->getMp();
			foreach(newStIter, newStateMp) {
				StatePtr nodeNextSt = (StatePtr) newStIter->first; 
				LD bJoint = getBeliefUpdate(node, nodePastSt, nodeNextSt);
				sumB += bJoint; 
				node->bJointRev->setPr(nodeNextSt, nodePastSt, bJoint); 
			}
		}
		node->bNew->clear(); 
		ASSERT(sumB > 0.0); 
		// update bNew
		map<StatePtr, Belief, StatePtrCmp>* cbMpPtr2 = node->bJointRev->getCbMpPtr();
		foreach(mpIter, *cbMpPtr2) {
			StatePtr nodeNextSt = mpIter->first; 
			LD totalAlloc = node->bJointRev->getAllocPr(nodeNextSt); 
			node->bNew->setPr(nodeNextSt, (totalAlloc/sumB)); 	
		}
	}
} 

/**
 * NOTE: Assumes RegionPPM::bJointRev correctly initialized
 */
void GraphGBP::step_updateMpc() {
	foreach(mpIter, adj_map) {
		RegionPPM* child = (RegionPPM*) mpIter->first;
		if(mpIter->second.second.size() == 0) {
			// has no parents - top level nodes
			continue; 
		} else {
			vector<MyEdge>& edgeList = mpIter->second.second; 
			foreach(parEdge, edgeList) {
				// note: ends[1] == child
				RegionPPM* parent = (RegionPPM*) parEdge->ends[0]; 
				map<StatePtr, Belief, StatePtrCmp>* cbMpPtr2 = child->bJointRev->getCbMpPtr();
				foreach(mpIter, *cbMpPtr2) {
					StatePtr nodeNextSt = mpIter->first; 
					BeliefPtr nStBlf = child->bJointRev->getBeliefForState(nodeNextSt);
					map<StatePtr, LD, StatePtrCmp> newStateMp = nStBlf->getMp();
					foreach(stIter, newStateMp) {
						StatePtr nodePastSt = stIter->first; 
						StatePtr jointSt = getJointState(nodePastSt, nodeNextSt);
						LD mNew = getMpcUpdate(parent, child, jointSt);
						updateMsg(parent, child, jointSt, mNew, false);  
					}
				} 	
			}
		}
	}
}

/** Revision Log
 *
 * July 07 2008. Initialized class
 * July 08 2008. Added construct_graph(), step_init() 
 * July 09 2008. Added lambda() functions; \f$m_alpha\f$ initialization 
 * July 16 2008. 
 */
