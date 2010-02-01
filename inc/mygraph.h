#ifndef MYGRAPH_H_
#define MYGRAPH_H_
#include <iostream>
#include <vector>
#include <map>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "node.h"
#include "myedge.h"

#include "main.h"
using namespace std; 
const LD SPIN_FLIP_RATE=0.5;
const LD DELTA_T=1.0; 
/// used to initialize messages in case has value zero 
const LD EPSILON_ZERO=1e-9; 
typedef vector<MyEdge> VE; 
typedef pair<vector<MyEdge>, vector<MyEdge> > BidirEdgeList; 
typedef map<NodePtr, BidirEdgeList, NodeCmp>::iterator ADJ_MAP_ITER; 

/**
 * Generic graph data structure allowing both factor and region graphs
 */

class MyGraph {
    int numNodes;
  
 
    vector<INT> values; // values each variable takes; 
    vector<string> vars; // the variable nodes in the graph
    
public:
	/// returns the pair <currentState, nextState> based on the jointState
	pair<StatePtr, StatePtr> getChildFromMsgState(StatePtr sptr);
    map<StatePtr, int, StatePtrCmp> sMp;
    int numEdges;
    /// list of associated edges which begin at the current node
    map<NodePtr, BidirEdgeList, NodeCmp> adj_map;   
    StatePtr getStatePtr(vector<INT> vals); 
    StatePtr getStatePtr(INT v) {
        vector<INT> vals; 
        vals.push_back(v); 
        return getStatePtr(vals);    
    }
    MyGraph() {
        numNodes = 0;
        numEdges = 0; 
    }
    int getDegree(NodePtr np) {
     map<NodePtr, BidirEdgeList, NodeCmp>::iterator mIt =  adj_map.find(np);   
        ASSERT(mIt != adj_map.end()); 
        int df = mIt->second.first.size();
        int db = mIt->second.second.size(); 
        return (df + db); 
    
    }
    
    pair<EdgePtr, EdgePtr> addEdge(const NodePtr& nps, const NodePtr& npt, bool checkExisting = false); 
    pair<EdgePtr, EdgePtr> findEdge(NodePtr nps, NodePtr npt);
    //void addEdge(Node nps, Node npt); 
    NodePtr addNode(NodePtr np); 
    inline void setVarVals(vector<INT> v) { values = v; }
    inline vector<INT> getVarVals() const { return values; } 
    inline void setVarIds(vector<string> v) { vars = v; }
    inline vector<string> getVarIds() const { return vars; } 
    void updateMsg(NodePtr parent, NodePtr child, StatePtr num, LD value, bool forwardMsg=true);
    
    
    vector<NodePtr> getParents(NodePtr node); 
    vector<NodePtr> getChildren(NodePtr node); 
    /**
     * function that returns the message from child->parent when childToParentMsg is true 
     * else returns parent->child reverse edge message. jointState gives the child joint state 
     */
    LD getMsg(NodePtr parent, NodePtr child, StatePtr jointState, bool forwardMsg=true); 
    
    NodePtr getNode(string nid);  
    NodePtr getNode(NodePtr np); // actually goes and searches the adj_map
    ~MyGraph() {
        // foreach(iter, sMp) delete(iter->first); 
      //  foreach(iter, adj_map) delete(iter->first); 
    }
    friend ostream& operator<<(ostream& os, MyGraph& g) {
        map<NodePtr, BidirEdgeList, NodeCmp>::iterator itB = g.adj_map.begin();
        map<NodePtr, BidirEdgeList, NodeCmp>::iterator itE = g.adj_map.end();
         
        while(itB != itE) {
           os<<"Node: ";
         Node n = (Node) *(itB->first);
         os<<n;
            int locVar = 0; 
            VE vCurr = itB->second.first;
    TOP:
            VE::iterator vItCurr = vCurr.begin();  
          
            while(vItCurr != vCurr.end()) {
                MyEdge e = *(vItCurr);
                if(locVar == 1) os<<"*RE*";
                os<<"\t"<<e<<"\n";  
                vItCurr++;
            }      
            if(locVar == 0) { 
                locVar++;   
                vCurr = itB->second.second;
               // cout<<"Hi YTh";  
                goto TOP;
            }  
            itB++;
            
                
        }
        os<<"GRAPH: numNodes = "<<g.numNodes<<" numEdges = "<<g.numEdges<<"\n"; 
        return os;  
    }
};


/**
 * Graph structure corresponding to factor graphs
 *  - supports BP algorithm 
 *  - constrains factor node id to F_(list of var_ids separated by _) eg, F_X0.0_X0.1
 *  - constrains variable node id to X(nodeInt) eg, X0.0
 */
class FactorGraph: public MyGraph { 
    bool dirtyZ; 
    LD trueZ; 
  
    LD theta; // the bp constant
    map<long long, LD> trueProbMp;
public:
  bool dirtyF; 
    LD Fapprox;
    bool bp_iter(bool ch); 
    int getStPtrIdx(StatePtr sp) {
        long long idx=0;
        int n=0; 
        map<INT, INT> valMp; 
        vector<INT> varVals = getVarVals();
        fi(0, varVals.size()) valMp[varVals[i] ] = i;
        int Q=varVals.size(); 
        long long mlt=1;   
        for(int i=0; i < sp->num; i++) {
            idx += mlt*valMp[sp->getData(i)];
            mlt *= Q; 
        }
        return idx; 
    }
    FactorGraph() { 
        trueZ = 0.0; 
        dirtyZ = true;  
        dirtyF = true; 
        theta = 0.01;    
    } 
    LD getTrueZ();
    LD getTrueP(StatePtr s) ; 
    StatePtr getFactorStatePtr(StatePtr sComplete, string fctrNodeId); 
    inline bool isFactor(const NodePtr np) { 
        string id = np->getNodeId();
        return (id.at(0) == 'F'); 
    }  
    vector<INT> getIthStVec(int i, int N, vector<INT> vals); 
    /// updates the child to parent message
    LD updateN(NodePtr np, NodePtr nc);
    /// updates the parent to child message
    LD updateM(NodePtr np, NodePtr nc); 
    /// updates the belief at variable
    
    void updateMN(NodePtr np, NodePtr nc) {
        updateN(np, nc); 
        updateM(np, nc); 
     
        getchar(); 
    }
    bool updateB(NodePtr nv);   
    /// additionally initializes the msgMp for each edge
    pair<EdgePtr, EdgePtr> addEdge(NodePtr ns, NodePtr nt, bool checkExisting=false);
    // synchronous bp 
    void bp_synch(bool ch); 
    /// computes the marginal probability for given state
    LD getMarginalProb(vector<pair<string, int> > VVP); 
    
};
class RegionGraph; 
typedef FactorGraph* FactorGraphPtr;
typedef RegionGraph* RegionGraphPtr;
/**
 * Implements the region 
 */
class Region : public DynNode {
    
    /// the individual factor nodes
    Cluster cl;
    RegionGraphPtr rg; 
    FactorGraphPtr fg; 
 
    
public:
    /// counting numbers
    int cr; 
    /// the K(x0, x_t) factor
    ConditionalBeliefPtr B_t_0;
    /// the single ton partition function like function 
    BeliefPtr B_0, B_t;
    /// the partition function 
    BeliefPtr Z; 
    
    Region(RegionGraphPtr _fg, Cluster _cl, int _cr,  ConditionalBelief* _cb=0); 
    LD getStateBlf(const StatePtr sptr) const ;
    vector<string> getVarIds() const { return cl.first; }  
  
    friend ostream& operator<<(ostream& os , Region& r) {
        Node n = (Node) r; 
        os<<n;
        os<<"r.cr: "<<r.cr<<"\n";
        os<<"r.cl: "<<r.cl<<"\n"; 
        return os; 
    }
    /// -sum_b b*log(b) 
    LD getEntropy();
    /// sum_b b*(-sum_a log f_a) where a is in r
    LD getEnthalpy(); 
    LD getF() { return (getEnthalpy() - getEntropy()); } 
    
};

typedef Region* RegionPtr; 
/**
 * Graph structure for region graph
 * implements the parent-child algorithm
 */
class RegionGraph: public MyGraph {
    FactorGraphPtr graphPtr; 
public:
  //  map<Region, BidirEdgeList, NodeCmp> adj_map;   
    RegionGraph(FactorGraphPtr _fg) { graphPtr = _fg; } 
    FactorGraphPtr getFactorGraphPtr() const { return graphPtr; } 
    
    StatePtr getMargState(StatePtr sptr,  RegionPtr rgn, RegionPtr chld);
    ~RegionGraph() { delete(graphPtr); } 
    /*
    pair<EdgePtr, EdgePtr> addEdge(Region nps, Region npt, bool checkExisting = false) {
        return MyGraph::addEdge(nps, npt, checkExisting); 
    }
    
    RegionPtr addNode(Region ns) {
        RegionPtr res =(RegionPtr) MyGraph::addNode(ns); 
        return res; 
    }
    */  
    
}; 

#endif /*MYGRAPH_H_*/

/* REVISION LOG:
 * 
 * 08 July 08 - Added getParents(), getChildren, getMsg() to MyGraph class 
 */
