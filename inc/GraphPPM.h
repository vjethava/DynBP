
#ifndef GRAPHPPM_H_
#define GRAPHPPM_H_
#include "mygraph.h"
#include "node.h" 
#include "RegionPPM.h" 


/**
 * Graph that implements PPM with partial state space
 * 
 */
class GraphPPM : public MyGraph {
    LD getMargPr_5c(RegionPPM* vp, RegionPPM* vc, StatePtr cStOld, StatePtr cStNew); 
public:
    LD SR; // overall dual variable
    LD minMsg;
    LD maxMsg; 
    LD avgDiffFromPars; 
    LD avgChngInBlfs; 
    LD maxChngInBlfs; 
    LD maxDevFromPars;
    LD avgMaxCaseDevFromPars;
    LD avgMaxCaseChngInBlfs; 
    /// initialize all messages
    StatePtr getJointState(StatePtr currSt, StatePtr nextSt) {
        
        vector<INT> vals;
        fi(0, currSt->num) {
            vals.push_back(currSt->getData(i));
        }
        fi(0, nextSt->num) {
            vals.push_back(nextSt->getData(i));
        }
        StatePtr jointSt = getStatePtr(vals);
        return jointSt; 
    }
    void PPM_inner_loop(); 
    void PPM_step3();
    void PPM_step4a();
    void PPM_step4b(); 
    void PPM_step5a(); 
    void PPM_step5b(); 
    void PPM_step5c(); 
    virtual StatePtr getChildStatePtr(NodePtr vp,  NodePtr vc, StatePtr parState)=0; 
     
   
    void updateNodesInfo(); 
    string get_status();
    /// does the log normalization on the messages
    void PPM_LN(LD thetaLN=1.0); 
};



/// class that provides sort based on nearest child match 
struct ANN_Cmp {
    NodePtr vp;
    NodePtr vc;
    StatePtr cSt; 
    GraphPPM* graph; 
    ANN_Cmp(GraphPPM* _graph, NodePtr _parId, NodePtr _cId, StatePtr _sptr) {
        vc = _cId; 
        vp =  _parId; 
        cSt = _sptr;  
        graph = _graph; 
    }
    bool operator()(const pair<StatePtr, LD>& p1, const pair<StatePtr, LD>& p2) const {
        StatePtr cp1=NULL, cp2=NULL; 
        if(vc != NULL) {
            cp1 = graph->getChildStatePtr(vp, vc, p1.first); 
            cp2 = graph->getChildStatePtr(vp, vc, p2.first);
        } else {
            cp1 = p1.first;
            cp2 = p2.first; 
        } 
        
        ASSERT( (cp1->num == cp2->num) && (cp1->num == cSt->num)); 
        LD d1 = State::SSD(*cp1, *cSt); 
        LD d2 = State::SSD(*cp2, *cSt); 
        if(d1 == d2) return (p1.second > p2.second);
        else return (d1 < d2); 
    }
};


#endif /*GRAPHPPM_H_*/
