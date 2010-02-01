#include "node.h"

Node::Node(string s, int _numVar,  Belief* bp, int type) {
    //    assert(uids.find(s) == uids.end()); 
    id = s;
    numVar = _numVar;
    
    bPtr = bp;
//    uids.insert(make_pair<string, NodePtr>(s, this)); 
}



DynNode::DynNode(string _s, int _numVar, Belief* bp, ConditionalBelief* cbp, int type ) : Node(_s, _numVar, bp, type) 
{
  cPtr = cbp; 
}

// int Node::numVal=0; 
