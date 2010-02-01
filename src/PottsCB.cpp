#include "PottsCB.h"
#include "PottsGraph.h"
#include "state.h"
#include "VideoNode.h"
using namespace std;
LD PottsCB::theta = 0.5; 
LD PottsCB::dt = 1.0; 

LD PottsCB::getPr(StatePtr s1, StatePtr s2) {
        cout<<"PottsCB("<<node->getNodeId()<<" "<<*s1<<", "<<*s2<<"\n" ;
        ASSERT(node->getNodeNumVar() == s1->num); 
        ASSERT(s1->num == s2->num); 
        LD val = graph->get_p_hat(node, s1, s2, theta, dt); 
        return val; 
}   
