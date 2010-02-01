/**
 * Modification of GraphPPM.h to a separate class 
 * -used by both GraphPPM and GraphGBP  
 */
#ifndef REGIONPPM_H_
#define REGIONPPM_H_

#include "mygraph.h"
#include "node.h" 

class RegionPPM : public Node {
   int cr; 
public:
   	 
    /// single normalizing factor \f$m_\alpha\f$ for region \f$\alpha\f$  
    LD m_alpha; 
    
    
    /// normalizing factor per initial state used by GraphPPM 
    BeliefPtr lambdaR; 
    /// the conditional evolution model used by GraphGBP & GraphPPM 
    ConditionalBeliefPtr pHat; 
    /// new values for \f$b_\alpha^{t, t+\delta t}\f$ saved as \f$t+\delta t, t\f$ (in reverse) 
    ConditionalBeliefPtr bJointRev; 
    /// the new beliefs which are committed to region beliefs at commit stage
    BeliefPtr bNew; 
    /// record holder for last beliefs
    BeliefPtr bLast;
    /// Holds all alpha->parent contributions per alpha state, used by GraphPPM
    ConditionalBeliefPtr lambdaPar;
    /// Holds all child->alpha per alpha state (computed using child reduced state messages)
    ConditionalBeliefPtr lambdaChild;
    /// average(over all states) deviation in belief compared with parent normalized beliefs
    LD avgDevFromPar; 
    /// average(over all states) change in the belief compared with last iteration beliefs
    LD avgChngInBlfs; 
    /// maximum deviation in belief compared with parent normalized beliefs
    LD maxChngInBlfs;
    /// maximum(over all states) change in the belief compared with last iteration beliefs
    LD maxDevFromPar;
    /// region entropy  
    LD SR; 
    RegionPPM(string s, int numVar, Belief* bp=0, int _type=-1) : Node(s, numVar, bp, _type) { 
        lambdaR=new Belief(); 
        pHat=new ConditionalBelief(); 
        bJointRev=new ConditionalBelief(); 
        bNew=new Belief();
        bLast = new Belief(); 
        lambdaPar=new ConditionalBelief();
        lambdaChild=new ConditionalBelief(); 
        avgDevFromPar = 0.0; 
        avgChngInBlfs = 0.0; 
        maxChngInBlfs = 0.0; 
        maxDevFromPar = 0.0; 
        SR = 0.0;
        m_alpha = 1.0;  
    }
    inline void setCR(int _cr) { cr = _cr; }
    inline int getCR() { return cr; }    
    void clear() {
        bLast->clear(); 
        if(getBeliefPtr() != 0) getBeliefPtr()->clear(); 
        lambdaR->clear();  
        pHat->clear(); 
        bJointRev->clear(); 
        bNew->clear(); 
        lambdaPar->clear(); 
        lambdaChild->clear();
    }
    ~RegionPPM() {
        delete(bLast); 
        delete(lambdaR);  
        delete(pHat); 
        delete(bJointRev); 
        delete(bNew); 
        delete(lambdaPar); 
        delete(lambdaChild);
    } 
};

#endif /*REGIONPPM_H_*/

/* REVISION LOG
 * ============
 * 09-08-2008 Added m_alpha which is used by GraphGBP as normalizing factor 
 * 
 */ 