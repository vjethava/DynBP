#include "mygraph.h"
#include <cassert>
#include <vector>
#include <gsl/gsl_sf_exp.h>

using namespace std;

/// returns the pair <currentState, nextState> based on the jointState
pair<StatePtr, StatePtr> MyGraph::getChildFromMsgState(StatePtr sptr) {
    int num = sptr->num;
    ASSERT(num%2 == 0);
    vector<INT> cval, nval;
    fi(0, num/2) {
        cval.push_back(sptr->getData(i));
    }
    fi(num/2, num) {
        nval.push_back(sptr->getData(i));
    }
    StatePtr cst = getStatePtr(cval);
    StatePtr nst = getStatePtr(nval);
    pair<StatePtr, StatePtr> result = make_pair<StatePtr, StatePtr>(cst, nst);
    // cout<<"getChildFromMsg() msg: "<<*sptr<<"\n\t\tcst: "<<*cst<<"\n\t\t nst: "<<*nst<<endl;
    return result;
}

/// returns the list of parents of a node  
vector<NodePtr> MyGraph::getParents(NodePtr np) {
	map<NodePtr, BidirEdgeList, NodeCmp>::iterator myIter = adj_map.find(np); 
	vector<NodePtr> result; 
	if(myIter != adj_map.end()) { // iterator over the edges
		BidirEdgeList& bev = myIter->second;
		foreach(edge, bev.second) {
			result.push_back(edge->ends[0]);
		}
	}
	return result; 
}




/// returns the list of children of a node  
vector<NodePtr> MyGraph::getChildren(NodePtr np) {
	map<NodePtr, BidirEdgeList, NodeCmp>::iterator myIter = adj_map.find(np); 
	vector<NodePtr> result; 
	if(myIter != adj_map.end()) { // iterator over the edges
		BidirEdgeList& bev = myIter->second;
		foreach(edge, bev.first) {
			result.push_back(edge->ends[1]);
		}
	}
	return result; 
}


/**
 * Adds node ( checks against duplicate node insertion )
 */
NodePtr MyGraph::addNode(NodePtr n) {
    map<NodePtr, BidirEdgeList, NodeCmp>::iterator myIter = adj_map.find(n);
    NodePtr result = 0;
    if(myIter != adj_map.end()) {
        result = (NodePtr) (myIter->first);
    } else {
        VE ve;
        BidirEdgeList bev = make_pair<VE, VE >(ve, ve);
        adj_map.insert(make_pair<NodePtr, BidirEdgeList >(n, bev));
        myIter = adj_map.find(n);
        result = (NodePtr) (myIter->first);
        string nid = result->getNodeId();
        //     npMap.insert(make_pair<string, NodePtr>(nid, result));
        numNodes++;
    }
    return result;
}

StatePtr MyGraph::getStatePtr(vector<INT> vals) {
    int num = vals.size();
    INT data[num];
    for(int i=0; i < num; i++) {
        data[i] = vals[i];
        //cout<<"val: "<<vals[i]<<" data: "<<data[i]<<"\n";  
    }    
    StatePtr sp = new State(num, data);
    typedef map<StatePtr, int, StatePtrCmp>::iterator ITER;
    ITER sIt = sMp.find(sp);
    if(sIt == sMp.end()) {
        pair<ITER, bool> ires = sMp.insert(make_pair<StatePtr, int>(sp, 0));
        /*
        if(ires.second) { // insertion took place
            //  FPRINTF(stderr, "INSERTED INTO sMp: "); cerr<<(*(ires.first->first))<<"\n"; 
          } else { 
              FPRINTF(stderr, "ALREADY PRESENT IN sMp: "); cerr<<(*(ires.first->first))<<"\n"; 
            sleep(2); 
        }
        */ 
        sIt = ires.first;
    } else {
        delete(sp); 
    }
    sIt->second += 1; 
    return sIt->first;    
}
/**
 * Adds an edge to the graph
 * does not check for edge already present condition
 */
pair<EdgePtr, EdgePtr> MyGraph::addEdge(const NodePtr& ns, const NodePtr& nt, bool checkExisting) {
    //void MyGraph::addEdge(Node ns, Node nt) {
    NodePtr nps = addNode(ns);
    NodePtr npt = addNode(nt);
    bool found = false;
    pair<EdgePtr, EdgePtr> result;
    map<NodePtr, BidirEdgeList,NodeCmp >::iterator itS, itT;
    itS = adj_map.find(ns);
    itT = adj_map.find(nt);
    if(checkExisting) {
        EdgePtr e1 = 0, e2 = 0;
        int n1 = itS->second.first.size();
        int n2 = itT->second.second.size();
        int i1=0,i2=0;
        while( (e1 == 0) && (i1 < n1) ) {
            NodePtr cs = itS->second.first[i1].ends[0];
            NodePtr ct = itS->second.first[i1].ends[1];
            ASSERT(cs == nps);
            if(ct == npt) {
                e1 = &(itS->second.first[i1]);
            } else
                i1++;
        }

        while( (e2 == 0) && (i2 < n2) ) {
            NodePtr cs = itT->second.second[i2].ends[0];
            NodePtr ct = itT->second.second[i2].ends[1];
            ASSERT(ct == npt);
            if(cs == nps) {
                e2 = &(itT->second.second[i2]);
            } else
                i2++;
        }

        ASSERT( ((e1==0)&&(e2==0)) ||((e1!=0)&&(e2!=0)));
        if(e1 != 0) {
            //   FPRINTF(stderr, " Found edge existing: ("); cerr<<ns<<", "<<nt<<")\n";
            result = make_pair<EdgePtr, EdgePtr>(e1, e2);
            found = true;
        }

    }

    if(!found) {
        MyEdge edge(nps, npt);
        itS->second.first.push_back(edge);
        itT->second.second.push_back(edge);

        numEdges++;
        int sz1 = itS->second.first.size();
        int sz2 = itT->second.second.size();

        EdgePtr e1 = (EdgePtr) &(itS->second.first[sz1-1]);
        EdgePtr e2 = (EdgePtr) &(itT->second.second[sz2-1]);
        result = make_pair<EdgePtr, EdgePtr>(e1, e2);
    }
    return result;
}

pair<EdgePtr, EdgePtr> MyGraph::findEdge(NodePtr nps, NodePtr npt) {

    ASSERT( (nps != 0) && (npt != 0));
    NodePtr ns = nps;
    NodePtr nt = npt;
    pair<EdgePtr, EdgePtr> result;
    map<NodePtr, BidirEdgeList,NodeCmp >::iterator itS, itT;
    itS = adj_map.find(ns);
    itT = adj_map.find(nt);
    EdgePtr e1 = 0, e2 = 0;
    int n1 = itS->second.first.size();
    int n2 = itT->second.second.size();
    int i1=0,i2=0;
    while( (e1 == 0) && (i1 < n1) ) {
        NodePtr cs = itS->second.first[i1].ends[0];
        NodePtr ct = itS->second.first[i1].ends[1];
        ASSERT(cs == nps);
        if(ct == npt) {
            e1 = &(itS->second.first[i1]);
        } else
            i1++;
    }

    while( (e2 == 0) && (i2 < n2) ) {
        NodePtr cs = itT->second.second[i2].ends[0];
        NodePtr ct = itT->second.second[i2].ends[1];
        ASSERT(ct == npt);
        if(cs == nps) {
            e2 = &(itT->second.second[i2]);
        } else
            i2++;
    }

    ASSERT( ((e1==0)&&(e2==0)) ||((e1!=0)&&(e2!=0)));

    result = make_pair<EdgePtr, EdgePtr>(e1, e2);
    return result;

}


/**
 * Returns ptr to the node having same id ( or 0, if not found)
 */
NodePtr MyGraph::getNode(NodePtr np) {
    map<NodePtr, BidirEdgeList, NodeCmp>::iterator cIter = adj_map.find(np);
    if(cIter != adj_map.end()) {
        return ((NodePtr) (cIter->first));
    } else
        return 0;
}

/** 
 * Update the message between parent to child
 * mode = true : parent->child
 * mode = false : child -> parent
 */
void MyGraph::updateMsg(NodePtr parent, NodePtr child, StatePtr num, LD value, bool mode) {
    NodePtr pNode = parent; //(Node) *(parent);
    NodePtr cNode = child; //(Node) *(child);
    map<NodePtr, BidirEdgeList, NodeCmp>::iterator pIter, cIter;
    pIter = adj_map.find(pNode);
    cIter = adj_map.find(cNode);
    assert(pIter != adj_map.end());
    assert(cIter != adj_map.end());
    vector<MyEdge>::iterator peIter = pIter->second.first.begin();
    bool notFound = true;
    while(notFound && (peIter != pIter->second.first.end())) {
        // cout<<"updtMsg(): child"<<(*child)<<" peIter->ends[1]: "<< (*(peIter->ends[1]))<<"\n";
        if(child == peIter->ends[1]) {
            notFound = false;
        } else
            peIter++;
    }

    ASSERT(peIter != pIter->second.first.end());
    notFound = true;
    vector<MyEdge>::iterator ceIter = cIter->second.second.begin();
    while(notFound && (ceIter != cIter->second.second.end())) {
        if(parent == ceIter->ends[0]) {
            notFound = false;
        } else
            ceIter++;
    }
    ASSERT(ceIter != cIter->second.second.end());
    if(mode) { // update parent
        peIter->setMsg(num, value);
    } else { // update child - here num = StatePtr
        ceIter->setMsg(num, value);
    }

};

/**
 * This returns the message from parent->child (parent->bidirEdgeList.first) 
 * edge if forwardMsg==true else returns the reverse edge parent->child 
 * (child->bidirEdge.second) given the child state is jointSt 
 */ 
LD MyGraph::getMsg(NodePtr parent, NodePtr child, StatePtr jointSt, bool forwardMsg) {
	// assumes childSt is given as jointSt
	ASSERT(jointSt->num == 2*child->getNodeNumVar()); 
	if(forwardMsg) { 
		map<NodePtr, BidirEdgeList, NodeCmp>::iterator pIter = this->adj_map.find(parent); 
		ASSERT(pIter != adj_map.end()); 
		vector<MyEdge>& forwardEdges = pIter->second.first; 
		vector<MyEdge>::iterator edgeIter = forwardEdges.begin(); 
		while(edgeIter != forwardEdges.end()) {
			if(edgeIter->ends[1] == child) {
				MyEdge& edge = *edgeIter;
				LD message = edge.getMsg(jointSt);
				return message;  	
			}
			edgeIter++; 
		}	
	} else {
		map<NodePtr, BidirEdgeList, NodeCmp>::iterator cIter = this->adj_map.find(child); 
		ASSERT(cIter != adj_map.end()); 
		vector<MyEdge>& reverseEdges = cIter->second.second; 
		vector<MyEdge>::iterator edgeIter = reverseEdges.begin(); 
		while(edgeIter != reverseEdges.end()) {
			if(edgeIter->ends[0] == parent) {
				MyEdge& edge = *edgeIter;
				LD message = edge.getMsg(jointSt);
				return message;  	
			}
			edgeIter++; 
		}
	}
	cerr<<"\nCould not find the message: parent: "<<parent->getNodeId();
	cerr<<" child: "<<child->getNodeId()<<endl;
	cerr<<" pNV: "<<parent->getNodeNumVar()<<" cNV: "<<child->getNodeNumVar();
	cerr<<" jStNV: "<<jointSt->num<<endl;
	exit(-1); 	 
}

/**
 * Wrapper to provide node search using node id 
 */
NodePtr MyGraph::getNode(string nid) {
    Node node(nid, 1, 0);
    return getNode(&node);
}

Region::Region(RegionGraphPtr _rg, Cluster _cl, int _cr,  ConditionalBelief* _cb) {
    //FPRINTF(stderr, "Region(): "); cerr<<" cl: "<<_cl<<" cr: "<<_cr<<"\n"; fflush(stderr);
    rg = _rg;
    fg = rg->getFactorGraphPtr();
    cl = _cl; // the cluster
    vector<INT> varVals = fg->getVarVals();
    int Q = fg->getVarVals().size(); // number of values each variable takes
    int N = cl.first.size(); // the number of variables in the factor graph

    //FPRINTF(stderr, "Region: N=%d Q=%d vals: [", N, Q); foreach(iter, varVals) FPRINTF(stderr, "%d " , (*iter)); FPRINTF(stderr, "]\n");
    cr = _cr; // counting number
    setConditionalBelief(_cb);
    setNodeId(cl.getRegionId());
 //   BeliefPtr bptr = new Belief();
    vector<INT> vals;
    LD Z=0.0;
    int mm = 1;
    BeliefPtr newB = new Belief();
    for(int i=0; i < N; i++)
        mm = mm*Q;
    for(int i=0; i < mm; i++) { // set belief for each state
        int k = i;
        vals.clear();
        for(int j=0; j < N; j++) {
            int id = k%Q;
            int vl = varVals[id];
            vals.push_back(vl);
            k = k/Q;
        }
        StatePtr sptr = rg->getStatePtr(vals);
        LD blf = getStateBlf(sptr);
        newB->setPr(sptr, blf);
        Z += blf;
        // cerr<<i<<" mm: "<<mm<<" " <<*sptr<<"  blf: "<<blf<<" Z: "<<Z<<" bptr->getPr(sptr): "<<bptr->getPr(sptr)<<" allocPr: "<<bptr->getAllocPr()<<"\n"; // sleep(2);
    }
    map<StatePtr, LD, StatePtrCmp> mp = newB->getMp();
  /*  foreach(iter, mp) {
        bptr->setPr(iter->first, ((iter->second)/Z));
    }*/ 
    //FPRINTF(stderr, " belief: "); cerr<<(*bptr);
    setBelief(newB);

}

LD Region::getStateBlf(const StatePtr sptr) const {
    LD p = 1.0;
    //FPRINTF(stderr, "\ngetStateBlf: %s state: ", getNodeId().c_str());  cerr<<*sptr<<"\n"; fflush(stderr);
    foreach(iter, cl.second) { // iterate over the factors
        //  FPRINTF(stderr, "\n\t FACTOR: %s ", iter->c_str());
        vector<string> fvars = split(*iter, '_');

        vector<INT> cvals;
        for(int i=1; i < (int) fvars.size(); i++) {
            //    FPRINTF(stderr, " var: %s", fvars[i].c_str());
            int idx=0;

            while(cl.first[idx] != fvars[i]) {
                idx++;
                ASSERT(idx < (int) cl.first.size());
            }
            int valOfCurrVar = sptr->getData(idx);
            cvals.push_back(valOfCurrVar);
        }
        StatePtr currFctrStPtr = fg->getStatePtr(cvals);
        NodePtr currFctrNode = fg->getNode(*iter);
        ASSERT(currFctrNode != 0);
        LD currFctrP = currFctrNode->getBeliefPtr()->getPr(currFctrStPtr);
        ASSERT(currFctrP != BELIEF_NOT_FOUND);
        //     FPRINTF(stderr, " fctrSt: "); cerr<<(*currFctrStPtr)<<" pCurr: "<<currFctrP<<" p: "<<p;
        p *= currFctrP;
    }
    //FPRINTF(stderr, "\np=%lf\n", p);   sleep(10);
    return p;
}

StatePtr FactorGraph::getFactorStatePtr(StatePtr sComplete, string fctrNodeId) {
    vector<string> allVars = getVarIds();
    ASSERT(sComplete->num == (int) allVars.size());
    vector<string> fvars = split(fctrNodeId, '_');
    fvars.erase(fvars.begin()); // erase the F at the beginning
    vector<INT> stVal;
    foreach(iter, fvars) {
        int idx=0;
        while( (allVars[idx] != (*iter)) && (idx < (int)allVars.size()) ) {

            idx++;

        }
        ASSERT(idx < allVars.size());
        int cval = sComplete->getData(idx);
        stVal.push_back(cval);
    }

    StatePtr result = getStatePtr(stVal);
    /* FPRINTF(stderr, "sComplete: ");
     fi(0, allVars.size()) {
         FPRINTF(stderr, "%s=%d ", allVars[i].c_str(), sComplete->getData(i));
     }
     FPRINTF(stderr, " %s sCurr: ", fctrNodeId.c_str());
     fi(0, fvars.size()) {
         FPRINTF(stderr, "%s=%d ", fvars[i].c_str(), result->getData(i));
     }
     FPRINTF(stderr, "\n");
     */
    return result;
}

LD FactorGraph::getTrueZ() {

    LD Z =0.0;
    if(!dirtyZ)
        Z = this->trueZ;
    else {
        printf("computing Z: ");
        fflush(stdout);
        int N = getVarIds().size();
        int Q = getVarVals().size();
        int mm = 1;
        for(int i=0; i < N; i++) mm = mm*Q;
        printf("mm: %d", mm);  
        fi(0, mm) { // enumerate the states
            if(i%1000 == 0) {
                printf(".");
                fflush(stdout);
            }
            StatePtr si = getStatePtr(getIthStVec(i, N, getVarVals()));
            Z += getTrueP(si);
        }
        dirtyZ = false;
        this->trueZ = Z;
        printf(" finished\n");
    }

    return Z;
}
LD FactorGraph::getTrueP(StatePtr state) {
    int idx = getStPtrIdx(state);
    map<long long, LD>::iterator tpIter = trueProbMp.find(idx);
    LD p=1.0;
    bool calc = false;
    if(tpIter == trueProbMp.end()) {
        calc = true;
        int N = state->num;
        ASSERT(N == (int)getVarIds().size());

        foreach(iter, adj_map) {
            NodePtr nptr = (NodePtr) (iter->first);
            string id = nptr->getNodeId();
            LD cf = 1.0;
            if(isFactor(nptr)) {
                StatePtr sf = getFactorStatePtr(((StatePtr)state), ((string)id));
                cf = nptr->getBeliefPtr()->getPr(sf);
                p = p*cf;
                FPRINTF(stderr, "getTrueP() p: %f n: %s cf: %f\n", p, nptr->getNodeId().c_str(), cf);
            }
        }
   //     StatePtr ns=new State(*state);
        trueProbMp.insert(make_pair<long long, LD>(idx, p));
    } else {
        p = tpIter->second;
    }
    // FPRINTF(stderr, "getTrueP(");
    // cerr<<(*state)<<") calc:  "<<calc<<" p = "<<p<<"\n";
    return p;
}

vector<INT> FactorGraph::getIthStVec(int idx, int N, vector<INT> vals) {
    int Q = vals.size();
    vector<INT> cval;
    for(int j=0; j < N; j++) {
        int id = idx%Q;
        int vl = vals[id];
        cval.push_back(vl);
        idx = idx/Q;
    }
    return cval;
}


/// updates the child to parent message
LD FactorGraph::updateN(NodePtr np, NodePtr nc) {
    FPRINTF(stderr, "updateN(%s, %s):\n", np->getNodeId().c_str(), nc->getNodeId().c_str());
    vector<INT> varVals = getVarVals();
    map<StatePtr, LD, StatePtrCmp> newN;

    foreach(valIter, varVals) {
        StatePtr csptr = getStatePtr(*valIter); // state for this value
        newN.insert(make_pair<StatePtr, LD>(csptr, 1.0));
    }

    map<NodePtr, BidirEdgeList, NodeCmp>::iterator cIter = adj_map.find(nc);
    ASSERT(cIter != adj_map.end());
    bool match = false;
    foreach(edgeIter, cIter->second.second) { // iterate over back edges
        if(edgeIter->ends[0] == np) {
            match = true;
        } else { // other edges
            pair<EdgePtr, EdgePtr> currEdges = findEdge(edgeIter->ends[0], edgeIter->ends[1]);
            EdgePtr fwdEdge = currEdges.first;
            foreach(valIter, varVals) {
                StatePtr csptr = getStatePtr(*valIter); // state for this value
                LD currFwdM = fwdEdge->getMsg(csptr); // message for this state
                map<StatePtr, LD, StatePtrCmp>::iterator itCurr = newN.find(csptr);
                ASSERT(itCurr != newN.end());
                itCurr->second = itCurr->second*currFwdM;
                FPRINTF(stderr, "\t%s val: %d m: %g newN: %g\n", fwdEdge->ends[0]->getNodeId().c_str(), *valIter, currFwdM, itCurr->second);
            }

        }
    }
    ASSERT(match == true);
    pair<EdgePtr, EdgePtr> orig_edges = findEdge(np, nc);
    orig_edges.second->msgMp.clear();
    orig_edges.second->msgMp.insert(newN.begin(), newN.end());

}


/// updates the parent to child message
LD FactorGraph::updateM(NodePtr np, NodePtr nc) {
    FPRINTF(stderr, "updateM(%s, %s):\n", np->getNodeId().c_str(), nc->getNodeId().c_str());
    vector<INT> varVals = getVarVals();
    vector<string> fvars = split(np->getNodeId(), '_');
    fvars.erase(fvars.begin()); // erase the beginning F
    string cvar = nc->getNodeId();
    map<StatePtr, LD, StatePtrCmp> newN;
    int Q = varVals.size();
    int N = fvars.size();
    foreach(valIter, varVals) {
        StatePtr csptr = getStatePtr(*valIter); // state for this value
        LD currM=0.0;
        int mm=1;
        for(int i=0; i < (N-1); i++)
            mm*=Q;
        for(int i=0; i < mm; i++) {
            LD locM = 1.0;
            vector<INT> cval = getIthStVec(i, N-1, varVals);
            vector<INT> oval = cval;
            // code that inserts x_i at its proper position
            int idx=0;
            vector<INT>::iterator cIt = cval.begin();
            while(fvars[idx] != cvar) {
                idx++;
                cIt++;
                ASSERT(idx < fvars.size());
            }
            cval.insert(cIt, *valIter);

            // cerr<<"\tvi: "<<*valIter<<" os: "<<getStr<int>(oval)<<" ns: "<<getStr<int>(cval)<<"\n"; sleep(10);

            StatePtr fSt=getStatePtr(cval);
            LD f_a = np->getBeliefPtr()->getPr(fSt);
            locM *= f_a;
            map<NodePtr, BidirEdgeList, NodeCmp>::iterator cIter = adj_map.find(np);
            ASSERT(cIter != adj_map.end());
            bool match = false;
            FPRINTF(stderr, "\t v: %d cval: %s f_a: %g\n", *valIter, getStr(cval).c_str(), f_a);
            foreach(edgeIter, cIter->second.first) { // iterate over forward edge
                if(edgeIter->ends[1] == nc) {
                    match = true;
                } else {
                    pair<EdgePtr, EdgePtr> curr_edges =  findEdge(np, edgeIter->ends[1]);
                    EdgePtr revEdge = curr_edges.second;
                    string child_id = edgeIter->ends[1]->getNodeId();
                    int idx=0;
                    while(fvars[idx] != child_id) {
                        idx++;
                        ASSERT(idx < fvars.size());
                    }
                    int child_val = cval[idx];
                    StatePtr child_ptr = getStatePtr(child_val);
                    LD n_child = revEdge->getMsg(child_ptr);
                    locM *= n_child;
                    FPRINTF(stderr, "\t\tcId: %s val: %d n_child: %g locM: %g\n", child_id.c_str(), child_val, n_child, locM);
                    // sleep(2);
                }
            }
            currM += locM;
            ASSERT(match==true);
        }
        FPRINTF(stderr, "\tval: %d currM: %g\n", *valIter, currM);
        newN.insert(make_pair<StatePtr, LD>(csptr, currM));
    }

    pair<EdgePtr, EdgePtr> edges=findEdge(np, nc);
    EdgePtr fwdEdge = edges.first;
    fwdEdge->msgMp.clear();
    fwdEdge->msgMp.insert(newN.begin(), newN.end());
}
/// updates the belief at variable
bool FactorGraph::updateB(NodePtr nv) {
    ASSERT(isFactor(nv) == false);
    LD totDiff=0.0;
    map<NodePtr, BidirEdgeList, NodeCmp>::iterator cIter = adj_map.find(nv);
    FPRINTF(stderr, "updateB(%s)\n", nv->getNodeId().c_str());
    ASSERT(cIter != adj_map.end());
    LD Z=0.0;
    vector<INT> varVals = getVarVals();
    Belief b_old = *(nv->getBeliefPtr());
    Belief b_new;
    bool metNan=false;
    // update the beliefs
    foreach(valIter, varVals) {
        LD p=1.0;
        StatePtr scp = getStatePtr(*valIter);
        foreach(edgeIter, cIter->second.second) {
            pair<EdgePtr, EdgePtr> curr_edges = findEdge(edgeIter->ends[0], edgeIter->ends[1]);
            LD ev = curr_edges.first->getMsg(getStatePtr(*valIter));
            p*=ev;
            FPRINTF(stderr, "\t\t%s m: %g p: %g\n", edgeIter->ends[0]->getNodeId().c_str(), ev, p);
        }
        FPRINTF(stderr, "\tval: %d p: %g\n", *valIter, p);
        // sleep(2);
        if(isnan(p))
            metNan = true;
        b_new.setPr(scp, p);

        Z += p;
    }
    // do not recompute should not reach nan
    if(Z == 0.0) {
        metNan = true;
    }
    if(!metNan) {
        foreach(valIter, varVals) {
            StatePtr scp = getStatePtr(*valIter);
            LD p_old = b_old.getPr(scp);
            LD p_new = b_new.getPr(scp);
            nv->getBeliefPtr()->setPr(scp, p_new/Z);
            totDiff = max(totDiff, abs(p_old-p_new));
            FPRINTF(stderr, "\tZ: %g val: %d p_old: %g p_new: %g totDiff: %g\n", Z, *valIter, p_old, p_new, totDiff);
        }
        // sleep(1);
        //getchar();
        if(totDiff > theta)
            return true;
        else
            return false;
    } else
        return false;
}


LD FactorGraph::getMarginalProb(vector<pair<string, int> > vvp) {
    vector<string> fvars = getVarIds();
    vector<INT> varVals = getVarVals();
    int N = fvars.size();
    int Q = varVals.size();
    int Nf = vvp.size();
    LD Z = getTrueZ();
    LD Psum=0.0;
    int mm=1;
    for(int i=0; i < (N-Nf) ; i++)
        mm*= Q;
    for(int i=0; i < mm; i++) {
        vector<INT> cval = getIthStVec(i, N-Nf, varVals);
        foreach(iter, vvp) {
            string cvar = iter->first;
            int val = iter->second;
            int idx=0;
            vector<INT>::iterator cIt = cval.begin();
            while(fvars[idx] != cvar) {
                idx++;
                cIt++;
                ASSERT(idx < fvars.size());
            }
            cval.insert(cIt, val);
        }
        StatePtr cPtr =  getStatePtr(cval);
        Psum += getTrueP(cPtr);
    }
    LD res = Psum/Z;
    return res;
}
/**
 * Adds an edge to the graph
 * does not check for edge already present condition
 * - Additionally initializes the msgMp of each edge with value 1
 */
pair<EdgePtr, EdgePtr> FactorGraph::addEdge(NodePtr ns, NodePtr nt, bool checkExisting) {
    //void MyGraph::addEdge(Node ns, Node nt) {
    NodePtr nps = ns; //addNode(ns);
    NodePtr npt = nt; //addNode(nt);
    bool found = false;
    pair<EdgePtr, EdgePtr> result;
    map<NodePtr, BidirEdgeList,NodeCmp >::iterator itS, itT;
    itS = adj_map.find(nps);
    itT = adj_map.find(npt);
    if(checkExisting) {
        EdgePtr e1 = 0, e2 = 0;
        int n1 = itS->second.first.size();
        int n2 = itT->second.second.size();
        int i1=0,i2=0;
        while( (e1 == 0) && (i1 < n1) ) {
            NodePtr cs = itS->second.first[i1].ends[0];
            NodePtr ct = itS->second.first[i1].ends[1];
            ASSERT(cs == nps);
            if(ct == npt) {
                e1 = &(itS->second.first[i1]);
            } else
                i1++;
        }

        while( (e2 == 0) && (i2 < n2) ) {
            NodePtr cs = itT->second.second[i2].ends[0];
            NodePtr ct = itT->second.second[i2].ends[1];
            ASSERT(ct == npt);
            if(cs == nps) {
                e2 = &(itT->second.second[i2]);
            } else
                i2++;
        }

        ASSERT( ((e1==0)&&(e2==0)) ||((e1!=0)&&(e2!=0)));
        if(e1 != 0) {
            FPRINTF(stderr, " Found edge existing: (");
            cerr<<ns<<", "<<nt<<")\n";
            result = make_pair<EdgePtr, EdgePtr>(e1, e2);
            found = true;
        }

    }

    if(!found) { // actually go and add a new edge
        MyEdge edge(nps, npt);
        itS->second.first.push_back(edge);
        itT->second.second.push_back(edge);

        numEdges++;
        int sz1 = itS->second.first.size();
        int sz2 = itT->second.second.size();

        EdgePtr e1 = (EdgePtr) &(itS->second.first[sz1-1]);
        EdgePtr e2 = (EdgePtr) &(itT->second.second[sz2-1]);
        vector<INT> vals = getVarVals();
        ASSERT(vals.size() > 0);
        for(int i=0; i < (int) vals.size(); i++) {
            StatePtr spCurr = getStatePtr(vals[i]);
            e1->setMsg(spCurr, 1.0);
            e2->setMsg(spCurr, 1.0);
        }
        result = make_pair<EdgePtr, EdgePtr>(e1, e2);
    }
    return result;
}


bool FactorGraph::bp_iter(bool ch) {
    dirtyF = true;
    bool redo=false;
    foreach(iter, adj_map) {
        if(isFactor((NodePtr) (iter->first))) {
            foreach(edgeIter, iter->second.first) {
                updateM(edgeIter->ends[0], edgeIter->ends[1]);
            }
        } else
            if(updateB((NodePtr) (iter->first))) {
                //              cout<< " belief changed at: "<<iter->first->getNodeId()<<"\n";
                redo = true;
            }
    }

    //  PRINTF("Press any key redo = %d\n", redo); getchar();
    //   cout<<(*this);

    // sleep(10);
    vector<string> varIds = getVarIds();
    vector<INT> varVals = getVarVals();


    if(redo) {
        foreach(iter, adj_map) {
            if(isFactor((NodePtr) (iter->first))) {
                foreach(edgeIter, iter->second.first) {
                    updateN(edgeIter->ends[0], edgeIter->ends[1]);
                }
            }
        }
    }
    return redo;
}

void FactorGraph::bp_synch(bool ch) {
    bool redo;
    int count=0;
LOOP:
    count++;
    //    cout<<"\n\nbp_synch(iteratn: "<<count<<")";
    redo = bp_iter(ch);
    if(redo)
        goto LOOP;
}


StatePtr RegionGraph::getMargState(StatePtr sptr,  RegionPtr rgn, RegionPtr chld) {

    /*
    vector<string> vs = split(subr, ' '); 
    if(vs[0] == "R") vs.erase(vs.begin());

    vector<string> cl = split(regn, ' '); 
    if(cl[0] == "R") cl.erase(cl.begin()); 
    cout<<"\n cl: "; foreach(it, cl) { cout<<" "<<(*it); }
    */
    vector<string> vs, cl;
    vs = rgn->getVarIds();
    cl = chld->getVarIds();
    // cout<<"\n getMargState() vs: "; foreach(it, vs) { cout<<" "<<(*it); }
    // cout<<"\n getMargState() cl: "; foreach(it, cl) { cout<<" "<<(*it); }
    // getchar();
    vector<INT> vals;
    int i=0, j=0;
    ASSERT(sptr->num == cl.size());
    while( (i < cl.size()) && (j < vs.size()) ) {
        //   cout<<"i: "<<i<<"j: "<<j<<" cl[i]: "<<cl[i]<<" vs[j]: "<<vs[j]<<"\n";
        if(cl[i] == vs[j]) {
            vals.push_back(sptr->getData(j));
            i++;
        }

        j++;

    }
    ASSERT(j == vs.size()); // found all
    StatePtr res = getStatePtr(vals);
    /*cerr<<" r_big: "<<regn<<"\tst: "<<*sptr<<"\n";
    cerr<<" r_sml: "<<subr<<"\tst: "<<*res;*/
    return res;
}

/// -sum_b b*log(b)
LD Region::getEnthalpy() {
    LD H = 0.0;
    map<StatePtr, LD, StatePtrCmp> bMp = B_0->getMp();
    BeliefPtr bfPtr = getBeliefPtr(); 
    //cout<<"rH: "<<getNodeId()<<"\n";
    foreach(it, bMp) {
        LD f = bfPtr->getPr(it->first); 
        H += -it->second*log(f);
      //  cout<<"\t"<<*it->first<<" b: "<<it->second<<" f: "<<f<<" H: "<<H<<"\n"; 
    }
    //getchar();
    return H;
    
}

/// sum_b b*(-sum_a log(f_a) ) where a is in r
LD Region::getEntropy() {
    LD U=0.0;
    map<StatePtr, LD, StatePtrCmp> bMp = B_0->getMp();
    foreach(it, bMp) {
        LD d = it->second;
        if( (d > 0.0) && (!isnan(d))) {
            U += -d*log(d);
        }
    }
    return U;
}


/* REVISION LOG
 * ============
 * July 08 - Moved getChildFromMsgState() from GraphPPM to MyGraph
 */

