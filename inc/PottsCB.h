#ifndef POTTSCB_H_
#define POTTSCB_H_
#include "belief.h"
#include "VideoNode.h"
#include "PottsModel.h" 
class PottsGraphPPM; 




class PottsCB: public ConditionalBelief {
    ConditionalBelief* cglobal;
    VideoNode* node; 
    PottsGraphPPM* graph; 
    static LD theta; 
    static LD dt; 
public:
    
    PottsCB(PottsGraphPPM* _graph, VideoNode* _node) {
        node = _node;
        graph = _graph; 
    } 
    static void setTheta(LD _theta) { theta = _theta; }
    static void setDt(LD _dt) { dt = _dt; }
    static LD getTheta() { return theta; }
    static LD getDt() { return dt; } 
    LD getPr(StatePtr s1, StatePtr s2); 
};

#endif /*POTTSCB_H_*/
