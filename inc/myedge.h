#ifndef MYEDGE_H_
#define MYEDGE_H_

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "main.h"
#include "node.h"
using namespace std;
 
const LD MSG_NOT_FOUND=-1.0;
struct MyEdge {
   
   NodePtr ends[2]; 
   map<StatePtr, LD, StatePtrCmp> msgMp; 
   MyEdge(NodePtr n1, NodePtr n2) {
        ends[0] = n1;
        ends[1] = n2;
       
   }
   LD getMsg(StatePtr val) {
        map<StatePtr, LD, StatePtrCmp>::iterator mIter = msgMp.find(val);
        if(mIter == msgMp.end()) {
            return MSG_NOT_FOUND; 
        } else 
            return mIter->second; 
   }
   
   inline void setMsg(StatePtr val, LD d) {
        map<StatePtr, LD, StatePtrCmp>::iterator miter = msgMp.find(val);
        if(miter != msgMp.end()) msgMp.erase(miter);  
        msgMp.insert(make_pair<StatePtr, LD>(val, d)); 
   }
   
   friend ostream& operator<<(ostream& os, MyEdge& e) {
        string ids = e.ends[0]->getNodeId(); 
        string idt = e.ends[1]->getNodeId(); 
        os<<"{"<<ids<<"->"<<idt;
        foreach(iter, e.msgMp) {
            os<<" "<<(*(iter->first))<<"="<<iter->second;
        }
        os<<" }";
        return os;
   }
};

typedef MyEdge* EdgePtr; 
class EdgeList {
    vector<MyEdge> edgeList; 
    
public:
    EdgeList(); 
    inline void addEdge(MyEdge e) { edgeList.push_back(e); }
    inline int getSize() { return edgeList.size(); }
    vector<MyEdge>::iterator begin() { return edgeList.begin(); }
    vector<MyEdge>::iterator end() { return edgeList.end(); }
    
};


#endif /*MYEDGE_H_*/
