#include "GraphBP.h" 
using namespace std; 

pair<EdgePtr, EdgePtr> GraphBP::addEdge(const NodePtr& ns, const NodePtr& nt, bool isParent, bool checkExisting) {
	NodePtr nps = addNode(ns);
    NodePtr npt = addNode(nt);
   // cerr<<"GraphBP::addEdge("<<*nps<<", "<<*npt<<")\n"; 
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
        if(e1 != 0) {
            //   FPRINTF(stderr, " Found edge existing: ("); cerr<<ns<<", "<<nt<<")\n";
            result = make_pair<EdgePtr, EdgePtr>(e1, e2);
            found = true;
        }

    }

    if(!found) {
        MyEdge edge(nps, npt);
        if(isParent)
        	itS->second.first.push_back(edge);
		else
			itT->second.second.push_back(edge); 
			
        numEdges++;
        int sz1 = itS->second.first.size();
        int sz2 = itT->second.second.size();

        EdgePtr e1 = (EdgePtr) &(itS->second.first[sz1-1]);
        result = findEdge(ns, nt);
        // result = make_pair<EdgePtr, EdgePtr>(e1, e2);
    }
    return result;
}


pair<EdgePtr, EdgePtr> GraphBP::findEdge(NodePtr nps, NodePtr npt) {
	// cerr<<"GraphBP::findEdge("<<(*nps)<<", "<<(*npt)<<")\n"; 
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
    result = make_pair<EdgePtr, EdgePtr>(e1, e2);
    return result;

}



/// updates the child to parent message
LD GraphBP::updateN(NodePtr np, NodePtr nc) {
    FPRINTF(stderr, "GraphBP::updateN(%s, %s):\n", np->getNodeId().c_str(), nc->getNodeId().c_str());
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
LD GraphBP::updateM(NodePtr np, NodePtr nc) {
    FPRINTF(stderr, "GraphBP::updateM(%s, %s):\n", np->getNodeId().c_str(), nc->getNodeId().c_str());
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
