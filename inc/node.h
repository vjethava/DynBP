#ifndef NODE_H_
#define NODE_H_

#include "main.h"
#include "belief.h"
#include <cmath>
#include <map>
/// node types
enum NodeType { VARIABLE, FACTOR, REGION };

class Node;
typedef Node* NodePtr; 
class DynNode;
typedef DynNode* DynNodePtr;
/**
 * The node class has the following convention regarding naming of the nodes
 * Variable nodes are named as VXXXX ( where XXXX represent numeric code) for example, V012
 * Factor nodes are named as F_VXXX_VXX ( F representing its a factor ) and  VXXX represent the 
 */    
class Node {
    //static int totalNodeCount; 
    
    /// set of nodes already present
    //static map<string, NodePtr> uids;
    
    /// the number of variables at this node
    int numVar;
    /// the number of values each variable can take
    // static int numVal;
    /// the unique string identifier for this node
    string id;
    /// set of beliefs about current state
    Belief* bPtr;
    /// the type of node
    NodeType type;
    /// used in graph algorithms 
    int color;

public:
   
    inline int getColor() const { return color; }
    inline void setColor(int _col) { color = _col; } 
    Node() { 
    //    totalNodeCount++;
        numVar = 1;
    //    numVal = 1; 
    //    int n = rand();
        id = "10000"; 
        
        bPtr = 0; 
    }
    inline void setNodeNumVar(int _n) { numVar = _n; }
  //  inline void setNodeNumVal(int _n) { numVal = _n; }
    inline void setNodeId(string _s) { id = _s; }
    inline void setBelief(Belief* b) { bPtr = b; }
 
    Node(string s, int numVar, Belief* bp=0, int _type=-1);
    inline int getNodeNumVar() const { return numVar; }
   // static int getNodeNumVal() { return numVal; }
   // static void setNodeNumVal(int _numVal) { numVal = _numVal; } 
    //static int getTotalNodeCount()  { return totalNodeCount; }
    inline string getNodeId() const { return id; }
    inline Belief* getBeliefPtr() { return bPtr; } 
    ~Node() { delete(bPtr); } 
    friend ostream& operator<<(ostream& os, Node& n) {
        os<<"<N: "<<n.id;
        if(n.bPtr != 0) {
            os<<"\nB: ";
            os<<(*(n.bPtr)); 
        } 
        os<<">\n"; 
        return os; 
    }
     
};


/// node having conditional forward time belief information
class DynNode : public Node {
    /// set of conditional beliefs about next state given current state
    ConditionalBelief* cPtr;
public:
    DynNode() { ; }
    inline void setConditionalBelief(ConditionalBelief* c) { cPtr = c; }
    DynNode(string s, int numVar, Belief* bp=0, ConditionalBelief* cbp=0, int _type=-1);
    inline ConditionalBelief* getConditionalBeliefPtr() { return cPtr; }
    friend ostream& operator<<(ostream& os, DynNode& dn) {
        Node* n = ((Node*)(&dn)); 
        os<<(*n); 
        /*
        if(dn.cPtr != 0) {
            os<<"CB: ";
            os<<(*(dn.cPtr));
        }*/ 
        return os; 
    }
    ~DynNode() { delete(cPtr); }
};

struct NodeCmp {
    bool operator()(const NodePtr& n1, const NodePtr& n2) const {
        return (n1->getNodeId() < n2->getNodeId()); 
    }
};

#endif /*NODE_H_*/
