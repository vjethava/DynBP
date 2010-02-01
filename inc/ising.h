#ifndef ISING_H_
#define ISING_H_
#include <fstream>
#include "mygraph.h"

/**
 * Graph having the linear Ising Model Structure ( nearest neighbour interactions only )
 */
class IsingGraph: public FactorGraph {
    int N1;
    int N2;
    map<string, LD> hMap; // single variable factors
    map<string, LD> jMap; // pair interactions
    map<string, NodePtr> npMap;
    gsl_rng *myRNG;
    LD H_sigma;
    LD J_sigma;
public:
    //IsingGraph() { ; }
    IsingGraph(int N1, int N2, int seed=0, LD H_sigma=1.0, LD J_sigma=1.0);
    bool isValid(int r, int c);
    void addSpinNodes();
    void addPairs();

    pair<int, int> getGraphSize() const {
        pair<int, int> res= make_pair<int, int>(N1, N2);
        return res;
    }
    void addEdge(NodePtr nps, NodePtr npt) {
        pair<EdgePtr, EdgePtr> pr = ((MyGraph*) (this))->addEdge(nps, npt);
        //((MyGraph*) (this))->addEdge(nps, npt);
        vector<INT> p0, p1;
        p0.push_back(0);
        p1.push_back(1);
        StatePtr s0 = getStatePtr(p0);
        StatePtr s1 = getStatePtr(p1);
        NodePtr pptr = getNode(nps);
        NodePtr cptr = getNode(npt);
        //     cout<<"parent: "<<(*pptr)<<" child: "<<(*cptr)<<"\n";
        updateMsg(pptr, cptr, s0, 1.0, false);
        updateMsg(pptr, cptr, s0, 1.0, true);
        updateMsg(pptr, cptr, s1, 1.0, true);
        updateMsg(pptr, cptr, s1, 1.0, false);
        ASSERT(pr.first != 0);
        ASSERT(pr.second != 0);
        //    cout<<"pr.first: "<<(*(pr.first))<<"\npr.second: "<<(*(pr.second))<<"\n";
    }
    LD getNodeEntropy(NodePtr np);
    LD getNodeEnthalpy(NodePtr np);
    LD getF();
    void bp_synch(bool ch); 
    ~IsingGraph() {  }
}
;

typedef IsingGraph* IsingGraphPtr;

class IsingRegionGraph: public RegionGraph {
    /// original ising graph
    IsingGraphPtr graphPtr;
    /// biggest region width
    int rx;
    /// biggest region size
    int ry;
    /// how much x-overlap between current and next region
    int cx;
    /// how much y-overlap between current and next region
    int cy;
    /// how many rows
    int N1;
    /// how many columns
    int N2;
    /// the clustering given by genMaximalRegions along with counting numbers
    map<Cluster, int, ClusterCmp> clusters;
    /// spin flip rate
    LD theta;
    /// the time of simulation
    LD delta_t;
    /// iteration count
    int iterNum;
public:
    IsingRegionGraph(IsingGraphPtr graph, int _rx, int _ry, int _cx, int _cy, LD theta=SPIN_FLIP_RATE, LD delta_t = DELTA_T);
    inline bool isValid(int r, int c) {
        return graphPtr->isValid(r, c);
    }
    Cluster genMaximalRegion(int x, int y); // collect set of maximal regions starting at (x, y)
    void genSubRegions(); // generate the smaller regions while still possible
    void initEdgeMsgs(pair<EdgePtr, EdgePtr> ep);
    LD getPrntM(const RegionPtr& rp, const StatePtr& sp);
    LD getChldM(const RegionPtr& rp, const StatePtr& sp0, const StatePtr& spT);
    LD getPrntM(const RegionPtr& rp, const StatePtr& sp0, const StatePtr& spT);
    ConditionalBeliefPtr getCondBlfForRegn(Region* np);
    void initBeliefs(RegionPtr& rp);
    void updateZ(RegionPtr& rp);
    LD updateB_t_0(RegionPtr& rp);
    LD updateB_t(RegionPtr& rp) {
        LD res = 0.0;
        map<StatePtr, LD, StatePtrCmp> bMp = rp->getBeliefPtr()->getMp();
        bool nanCheck=false; 
        foreach(iter, bMp) {
            LD b0 = rp->B_0->getPr(iter->first);
            LD bt = rp->B_t->getPr(iter->first);
            res = max(res, abs(b0 - bt));
            if(isnan(bt) || isnan(res)) nanCheck = true; 
        }
        if(!nanCheck) {
            foreach(iter, bMp) {
            LD b0 = rp->B_0->getPr(iter->first);
            LD bt = rp->B_t->getPr(iter->first);
            rp->B_0->setPr(iter->first, bt); 
            }
        } else   
            res = 0.0; 
        return res;
    }
    LD getMsgCtoP(EdgePtr e, StatePtr s0, StatePtr sT);
    bool bp_ppm(LD _th=SPIN_FLIP_RATE, LD _delta_t=DELTA_T);
    pair<LD, LD> update_edges(LD cFactor=1.0);
    // find the probability of variable = 1 from given region
    LD getMargPr(const RegionPtr& rp, string var);
    LD getF();
    bool renormalizeMessages(pair<LD, LD> maxMinMsgPair, LD lmt=1.0); 
    bool ppm_iter(LD _th, LD _delta_t, LD* minMsg, LD* maxMsg, LD* maxDev, LD* avgDev, map<string, vector<LD> >* nextBlfMpPtr=0, LD lmt=1.0);
    //  LD getApproxFreeEnergy();

};
#endif /*ISING_H_*/

