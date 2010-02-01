#include "mygraph.h"
#include <cassert>
#include <gsl/gsl_sf_exp.h>
#include "ising.h"
using namespace std;

/**
 * Factor graph having gaussian i.i.d. weights with nearest neighbour interactions
 */
IsingGraph::IsingGraph(int N1, int N2, int seed, LD H_sigma, LD J_sigma) {
    this->N1 = N1;
    this->N2 = N2;
    this->H_sigma = H_sigma;
    this->J_sigma = J_sigma;
    const gsl_rng_type* T;
    T = gsl_rng_default;
    myRNG = gsl_rng_alloc(T);
    gsl_rng_set(myRNG, seed);
    addSpinNodes();
    // can take only +1, -1
    vector<INT> nvals;
    nvals.push_back(0);
    nvals.push_back(1);
    setVarVals(nvals);
    addPairs();
}
/**
 * Adds the spin nodes to the required graph
 */
void IsingGraph::addSpinNodes() {
    StatePtr sptr[2];
    INT s[2] = {0, 1};
    for(int i=0; i < 2; i++) {
        //LD val = ((LD) i);
        // int j = (i==1)?1:0;
        sptr[i] = new State(1, &(s[i]));
    }
    LD pEach = 0.5;
    vector<string> varIds;
    for(int i=0; i < N1; i++) {
        for(int j=0; j < N2; j++) {
            stringstream ss(""), ssf("");
            ss<<"X"<<i<<"."<<j;
            string id = ss.str();
            ssf<<"F_"<<id;
            BeliefPtr bptr = new Belief();
            BeliefPtr fbptr = new Belief();
            LD hCurr = gsl_ran_gaussian(myRNG, H_sigma);
            LD pCurr =  gsl_sf_exp(hCurr)/(gsl_sf_exp(hCurr) + gsl_sf_exp(-hCurr));
            //   cout<<"\npCurr: "<<pCurr<<" hCurr: "<<hCurr;
            for(int k=0; k < 2; k++) {
                LD pToPut;
                if(k==0)
                    pToPut = pCurr;
                else
                    pToPut = (1.0 - pCurr);
                //     cout<<" k="<<k<<" pToPut: "<<pToPut<<" ";
                bptr->setPr(sptr[k], pEach); // initial belief equiprobable
                fbptr->setPr(sptr[k], pToPut); // factor based on random 'h'
                //       cout<<(*fbptr);
            }
            string fid = ssf.str();
            DynNodePtr nptr = new DynNode(id, 1, bptr);
            varIds.push_back(id);
            DynNodePtr fnptr = new DynNode(fid, 1, fbptr);
            //   cout<<(*fnptr)<<"\n";
            // function invokation from parent class
            addNode(nptr);
            addNode(fnptr);
            addEdge((NodePtr)fnptr,(NodePtr) nptr);
        }
    }
    this->setVarIds(varIds);
}

bool IsingGraph::isValid(int r, int c) {
    return ( (r >=0) && (r < N1) && (c>= 0) && (c < N2) );
}

/**
 * Adds the nearest neighbour interactions
 */
void IsingGraph::addPairs() {
    StatePtr sptr[4];
    INT d[4][2];
    for(int i=0; i < 4; i++) {
        // unsigned char modification
        d[i][0] = i/2;
        d[i][1] = i%2; 
     //   d[i][0] = 2*(i/2)-1;
     //   d[i][1] = 2*(i%2)-1;
        sptr[i] = new State(2, d[i]);
    }
    for(int i=0; i < N1; i++) {
        for(int j=0; j < N2; j++) {
            for(int k=0; k < 2; k++) {
                int xn, yn; // neighbour coords
                if(k==0) {
                    xn=i;
                    yn = j+1;
                } else {
                    xn = i+1;
                    yn = j;
                }
                if(isValid(xn, yn) ) { // have a factor between (i,j)-(xn, yn)
                    BeliefPtr bptr = new Belief();
                    LD jCurr = gsl_ran_gaussian(myRNG, J_sigma);
                    LD pPlus = gsl_sf_exp(jCurr);
                    LD pMinus = gsl_sf_exp(-jCurr);
                    LD pTotal = 2*(pPlus + pMinus);
                    for(int l=0; l < 4; l++) {
                        int sgn = -1*(sptr[l]->getData(0))*(sptr[l]->getData(1));
                        LD pRel = ((sgn > 0)?pPlus:pMinus);
                        LD pCurr = pRel/pTotal;
                        bptr->setPr(sptr[l], pCurr);
                    }
                    stringstream ss(""), ss1(""), ss2("");
                    ss1<<"X"<<i<<"."<<j;
                    string nid1 = ss1.str();
                    ss2<<"X"<<xn<<"."<<yn;
                    string nid2 = ss2.str();
                    ss<<"F_"<<nid1<<"_"<<nid2;
                    string fid = ss.str();
                    NodePtr nptr1 = getNode(nid1);
                    NodePtr nptr2 = getNode(nid2);
                    assert(nptr1 != 0);
                    assert(nptr2 != 0);
                    DynNodePtr fptr = new DynNode(fid, 2, bptr);
                    addNode(fptr);
                    addEdge(fptr, nptr1);
                    addEdge(fptr, nptr2);
                }
            }
        }
    }
}

/**
 * Region graph based on the underlying ising graph 
 */
IsingRegionGraph::IsingRegionGraph(IsingGraphPtr graph, int _rx, int _ry, int _cx, int _cy, LD _theta, LD _delta_t) :
RegionGraph(graph) {
    // ofstream edgeUpdateFileStream;
    // edgeUpdateFileStream.open("edgeupdates.txt", ios::trunc);
    //char str[80];
    //strftime(str, 80, "%X", tc );
    //edgeUpdateFileStream<<"code execution: "<<time(NULL)<<"\n";
    //edgeUpdateFileStream.close();
    iterNum=0;
    graphPtr = graph;
    rx = _rx;
    ry = _ry;
    cx = _cx;
    cy = _cy;
    this->delta_t = 0.0; // _delta_t;
    this->theta = 0.0; // _theta;

    map<Cluster, int, ClusterCmp> origCl, newCl;
    pair<int, int> fsize = graph->getGraphSize();
    N1 = fsize.first;
    N2 = fsize.second;
    //    cout<<"REGION_ISING ORIG GRAPH: "<<(*graphPtr)<<"\n";
    for(int i=0; i <  N1-cx; i += cx) {
        for(int j=0; j < N2-cy; j += cy) {
            Cluster res = genMaximalRegion(i, j);
            origCl.insert(make_pair<Cluster, int>(res, 1));
            // add TopLevel Regions
            RegionPtr nrp =  new Region(this, res, 1, 0);
            //     ConditionalBeliefPtr cptr = getCondBlfForRegn(newR);
            //newR->setConditionalBelief(cptr);
            initBeliefs(nrp);
            addNode(nrp);
        }
    }

    bool foundSubCluster = false;
    do {

        clusters.insert(origCl.begin(), origCl.end());
        foundSubCluster = false;
        foreach(iter, origCl) {
            map<Cluster, int, ClusterCmp>::iterator oIt = iter;
            oIt++;
            while(oIt != origCl.end()) {
                vector<string> nv;
                Cluster *clusterptr = new Cluster();
                bool validSubCluster = Cluster::intersection(clusterptr, oIt->first, iter->first);
                int cr = 1;
                if(validSubCluster) {
                    if(clusters.find(*clusterptr) == clusters.end()) { // have found new subcluster
                        foundSubCluster = true;
                    }
                    //    cout<<(*clusterptr);
                    //    cout<<"subset of\n\t"<<(oIt->first);
                    //    cout<<"\n\t"<<(iter->first);
                    //    cout<<"\nOVERALL SUBSET:\n";
                    if(newCl.find(*clusterptr) == newCl.end()) {
                        foreach(iterMp, clusters) {
                            if(Cluster::subset(*clusterptr, iterMp->first) ) {
                                //              cout<<"\t"<<iterMp->second<<"\t"<<iterMp->first;
                                if( (*clusterptr) != (iterMp->first) )
                                    cr = cr - iterMp->second;
                                //              cout<<"\n";
                            }

                        }
                        //      cout<<"cr Final: "<<cr<<"\n";
                        newCl.insert(make_pair<Cluster, int>(*clusterptr, cr));
                    }
                    // add the node to region graph

                    RegionPtr newR = new Region(this, *clusterptr, cr, 0);
                    //RegionPtr newR;
                    if(getNode(newR->getNodeId()) == 0) {
                        //     ConditionalBeliefPtr cptr = getCondBlfForRegn(newR);
                        //     newR->setConditionalBelief(cptr);
                        initBeliefs(newR);
                        addNode(newR);
                    } else {
                        newR = (RegionPtr) getNode(newR->getNodeId()); //   delete(tmpR);
                    }

                    // }
                    // addEdges from two parents
                    RegionPtr parentA = (RegionPtr) getNode(oIt->first.getRegionId());
                    RegionPtr parentB = (RegionPtr) getNode(iter->first.getRegionId());
                    ASSERT(parentA != 0);
                    ASSERT(parentB != 0);
                    // prevent self loops
                    if(parentA->getNodeId() != newR->getNodeId())  {
                        pair<EdgePtr, EdgePtr> ep = addEdge(parentA, newR, true);
                        initEdgeMsgs(ep);
                    }
                    if(parentB->getNodeId() != newR->getNodeId()) {
                        pair<EdgePtr, EdgePtr> ep = addEdge(parentB, newR, true);
                        initEdgeMsgs(ep);


                    }
                    

                } /*else {
                                                                                                                                    cout<<"INVALID: "; 
                                                                                                                                    cout<<" zero intersection\n"<<(oIt->first);
                                                                                                                                   cout<<"\n"<<(iter->first); 
                                                                                                                                }*/

                delete(clusterptr); 
                // all code above
                oIt++;
            }
        }
        /*
        cout<<" NEWCL: " ;
        foreach(iter, newCl) {
            cout<<"\n\t"<<iter->first<<" = "<<iter->second;
    }
        cout<<"\n foundSubCluster = "<<foundSubCluster<<"\n";
        */
        origCl = newCl;
        newCl.clear();

    } while(foundSubCluster);
}

/**
 * Generate the maximal regions 
 */
Cluster IsingRegionGraph::genMaximalRegion(int x, int y) {
    int i=0, j = 0;
    vector<string> vars, fctrs;
    while(isValid(x+i, y+j) && (i < rx)) {
        while(isValid(x+i, y+j) &&  (j < ry)) {
            stringstream ss("");
            ss<<"X"<<(x+i)<<"."<<(y+j);

            if(isValid(x+i+1, y+j) && ((i+1) < rx) ) { // horizontal bond
                stringstream ssf("");
                ssf<<"F_"<<"X"<<(x+i)<<"."<<(y+j)<<"_X"<<(x+i+1)<<"."<<(y+j);
                fctrs.push_back(ssf.str());
            }
            if(isValid(x+i, y+j+1) && ((j+1) < ry) ) { // horizontal bond
                stringstream ssf("");
                ssf<<"F_"<<"X"<<(x+i)<<"."<<(y+j)<<"_X"<<(x+i)<<"."<<(y+j+1);
                fctrs.push_back(ssf.str());
            }
            stringstream ssf1("");
            ssf1<<"F_"<<ss.str();
            fctrs.push_back(ssf1.str());
            vars.push_back(ss.str());

            // the j increment should be the last line
            j++;
        }

        // all code above
        i++;
        j=0;
    }
    /*
    PRINTF("vars: {");
    foreach(iter, vars) {
        PRINTF("%s ", iter->c_str());
}
    PRINTF("}\nfctrs: {");
    foreach(iter, fctrs) {
        PRINTF("%s ", iter->c_str());
}
    PRINTF("}\n");
    */
    sort(fctrs.begin(), fctrs.end());
    sort(vars.begin(), vars.end());
    Cluster res(vars, fctrs);
    return res;
}
/// this gives the conditional P(x_0 | x_dt) for the region
ConditionalBeliefPtr IsingRegionGraph::getCondBlfForRegn(Region* np) {
    BeliefPtr bptr = np->getBeliefPtr();
    map<StatePtr, LD, StatePtrCmp> bMp = bptr->getMp();
    ConditionalBeliefPtr cptr = new ConditionalBelief();
    foreach(iter1, bMp) {
        foreach(iter2, bMp) {
            // the energy term
            LD f_energy = iter2->second/iter1->second;
            int indc = 0;
            for(int i=0; i < iter1->first->num; i++) {
                if(iter1->first->getData(i) != iter2->first->getData(i))
                    indc++;
            }

            LD fctr = (theta*delta_t)/(1.0 - theta*delta_t);
            LD f_spins = pow(fctr, indc);
            LD overallFctr = f_spins*f_energy;
            cptr->setPr(iter1->first, iter2->first, overallFctr, iter1->second);
            //    cerr<<" s0: "<<*(iter1->first)<<" sT: "<<*(iter2->first)<<" fE: "<<f_energy<<" fS: "<<f_spins<<" fTot: "<<overallFctr<<"\n";
            //    getchar();
        }
    }
    return cptr;
}

/// initialize the messages as given
void IsingRegionGraph::initEdgeMsgs(pair<EdgePtr, EdgePtr> ep) {
    EdgePtr re = ep.second;
    RegionPtr rP = (RegionPtr) re->ends[0]; // parent region
    RegionPtr rC = (RegionPtr) re->ends[1]; // child region
    map<StatePtr, LD, StatePtrCmp> bMp = rC->getBeliefPtr()->getMp();
    vector<INT> vals;
    foreach(iter1, bMp) {


        foreach(iter2, bMp) {
            fi(0, iter1->first->num) {
                vals.push_back(iter1->first->getData(i));
            }
            fi(0, iter2->first->num) {
                vals.push_back(iter2->first->getData(i));
            }
            StatePtr sD = getStatePtr(vals);
            re->setMsg(sD, 1.0);
            vals.clear();
        }
    }
}

LD IsingRegionGraph::getPrntM(const RegionPtr& rp, const StatePtr& sp0, const StatePtr& spT) {
    State sChng = *sp0 + *spT;
    // cout<<" getPrntM("<<*rp<<", "<<(*sp0)<<", "<<(*spT)<<")\n"; fflush(stdout);
    return getPrntM(rp, &sChng);
}

LD IsingRegionGraph::getPrntM(const RegionPtr& rp, const StatePtr& sp) {
    LD result = 1.0;
    //  cout<<"getParentM("<<(*rp)<<"\nsp: "<<(*sp)<<")\n";
    map<NodePtr, BidirEdgeList, NodeCmp>::iterator iter =  adj_map.find(rp);
    ASSERT(iter != adj_map.end());
    // reverse edges
    vector<MyEdge>  ep = iter->second.second;
    // cout<<" ep.size(): "<<ep.size();
    foreach(iter, ep) {

        LD cMsg = iter->getMsg(sp);
        //     cout<<(*iter)<<"getPrntM() cMsg: "<< cMsg<<" result: "<<result<<"\n";

        ASSERT(cMsg != MSG_NOT_FOUND);
        result *= cMsg;
    }
    return result;
}

LD IsingRegionGraph::getChldM(const RegionPtr& rp, const StatePtr& sp0, const StatePtr& spT) {
    LD result = 1.0;
    map<NodePtr, BidirEdgeList, NodeCmp>::iterator iter =  adj_map.find(rp);
    ASSERT(iter != adj_map.end());
    vector<MyEdge> ep = iter->second.first;
    foreach(iter, ep) {
        pair<EdgePtr, EdgePtr> cpr = findEdge(iter->ends[0], iter->ends[1]);
        EdgePtr cedge = cpr.second; // found the reverse edge corresponding to forward edge
        // cout<<"getChldM() sp: "<<*sp0<<" child: "<<iter->ends[1]->getNodeId()<<"\n" ;        // getchar();

        RegionPtr rpT = (RegionPtr) iter->ends[1];
        vector<string> cl = rp->getVarIds();
        // cout<<"\ngetChildM_vars: "; foreach(it, cl) cout<<" "<<*it; cout<<"\n";
        StatePtr sc0 = getMargState(sp0, rp, rpT);
        StatePtr scT = getMargState(spT, rp, rpT);
        State sChng = *sc0 + *scT;
        double cMsg = cedge->getMsg(&sChng);
        //    cout<<" child: "<<(*rpT)<<" sp0: "<<*sp0<<" spT: "<<(*spT)<<" sChng: "<<(sChng)<<" cMsg: "<<cMsg<<"\n";
        result *= cMsg;
    }
    return result;
}

// initialize b_t, b_t_o
void IsingRegionGraph::initBeliefs(RegionPtr& rp) {
    BeliefPtr bfac = rp->getBeliefPtr();
    map<StatePtr, LD, StatePtrCmp> bMp = bfac->getMp();
    int N = bMp.size();
    LD bf0 = 1.0/( (double) N);
    LD bf1 = 1.0/( (double) (N*N));
    BeliefPtr b0 = new Belief();
    ConditionalBeliefPtr cb0 = new ConditionalBelief();
    BeliefPtr Z = new Belief();
    foreach(it1, bMp) {
        b0->setPr(it1->first, bf0);
        Z->setPr(it1->first, 1.0);
        foreach(it2, bMp) {
            cb0->setPr(it2->first, it1->first, bf1);
        }
    }
    rp->B_0 = b0;
    rp->B_t_0 = cb0;
    rp->Z = new Belief();
    rp->B_t = new Belief();

    // cout<<"\n"<<rp->getNodeId()<<" initB_0: "<<*(rp->Z);
}

void IsingRegionGraph::updateZ(RegionPtr& rp) {
    map<StatePtr, LD, StatePtrCmp> bMp = rp->getBeliefPtr()->getMp();

    foreach(it0, bMp) {
        LD Zs=0.0;
        foreach(itT, bMp) {
            StatePtr s0 = it0->first;
            StatePtr sT = itT->first;
            LD kCurr = rp->getConditionalBeliefPtr()->getPr(s0, sT);
            LD pContrib = getPrntM(rp, s0, sT);
            LD cContrib = getChldM(rp, s0, sT);
            LD powFac = 1.0/((double)rp->cr);
            LD term = kCurr*pow((pContrib/cContrib), powFac);
            Zs += term;
            //      cout<<*s0<<"->"<<*sT<<" ";
            //       printf("kCurr: %g pC: %g cC: %g powF: %g term: %g Zs: %g\n", kCurr, pContrib, cContrib, powFac, term, Zs);

        }
        //   cout<<"\n\t s0: "<<*(it0->first)<<" Z: "<<Zs<<"\n";
        rp->Z->setPr(it0->first, Zs);
    }

}
/// returns the maximum change relative
LD IsingRegionGraph::updateB_t_0(RegionPtr& rp) {
    LD Z=0.0;
    map<StatePtr, LD, StatePtrCmp> bMp = rp->getBeliefPtr()->getMp();
    ConditionalBeliefPtr b22 = new ConditionalBelief();
    Belief* b12 = new Belief();
    foreach(itT, bMp) {
        LD Zs=0.0;
        foreach(it0, bMp) {
            StatePtr sT = itT->first;
            StatePtr s0 = it0->first;
            LD bCurr = rp->B_0->getPr(s0);
            LD zCurr = rp->Z->getPr(s0);
            LD kCurr = rp->getConditionalBeliefPtr()->getPr(s0, sT);
            LD pContrib = getPrntM(rp, s0, sT);
            LD cContrib = getChldM(rp, s0, sT);
            LD powFac = 1.0/((double)rp->cr);
            LD term = kCurr*pow((cContrib/pContrib), powFac)*bCurr/zCurr;
            Zs += term;
            b22->setPr(s0, sT, term);
            // rp->B->setPr(s0, sT, term);
            // cout<<"upd_B_t_0 "<<*s0<<"->"<<*sT<<" ";
            // printf("kCurr: %g pC: %g cC: %g powF: %g bCurr: %g zCurr: %g term: %g Zs: %g\n", kCurr, pContrib, cContrib, powFac, bCurr, zCurr,  term, Zs);
            // getchar();
        }
        Z += Zs;
        b12->setPr(itT->first, Zs);
        // rp->B_t->setPr(itT->first, Zs);
        // cout<<"\n\t sT: "<<*(itT->first)<<" Z: "<<Zs<<"\n";
    }
    // normalized setting
    LD maxChange=0.0; 
    LD maxR=0.0; 
    LD oBlf=0.0; 
    LD nBlf=0.0; 
    StatePtr oSt; 
    if((!isnan(Z)) &&(Z  > 0.0)) {
        foreach(itT, bMp) {
            LD Zs = b12->getPr(itT->first);
            LD zOld = rp->B_t->getPr(itT->first); 
            foreach(it0, bMp) {
                StatePtr sT = itT->first;
                StatePtr s0 = it0->first;
                LD tcurr = b22->getPr(s0, sT);
                rp->B_t_0->setPr(s0, sT, tcurr/Z);
            }
            LD zDiff = abs(zOld-(Zs/Z));
           // if(Zs > 0.0) zDiff = zDiff/(Zs/Z); // ratio of change betwn old and new blf
            if(zDiff > maxChange) { 
                maxChange = zDiff;
                oBlf = zOld;
                oSt = itT->first;
                nBlf = Zs/Z; 
            }      
            rp->B_t->setPr(itT->first, Zs/Z);
            
         //   getchar(); 
        }
    }
    printf("updateB_t_0(%s)\nmChange: %g oBlf: %g nBlf: %g oSt: ", rp->getNodeId().c_str(), maxChange, oBlf, nBlf);
    cout<<*oSt<<"\n";  
    delete(b22);
    delete(b12);
    return maxChange; 
}
// assume the reverse edge
LD IsingRegionGraph::getMsgCtoP(EdgePtr e, StatePtr s0, StatePtr sT) {
    //ofstream edgeUpdateFileStream;
    //edgeUpdateFileStream.open("edgeupdates.txt", ios::app);
    State sJoint = *s0 + *sT;
    LD mOld = e->getMsg(&sJoint);
    LD mNew = 1.0;
    ASSERT(mOld != MSG_NOT_FOUND);
    RegionPtr rp = (RegionPtr) e->ends[0];
    RegionPtr rcp = (RegionPtr) e->ends[1];
    ASSERT( (s0->num == sT->num) && (s0->num == rcp->getVarIds().size()) );
    map<StatePtr, LD, StatePtrCmp> bMp = rp->getBeliefPtr()->getMp();
    LD Zs = 0.0;
    foreach(itT, bMp) {

        foreach(it0, bMp) {
            StatePtr sr0 = getMargState(it0->first, rp, rcp);
            StatePtr srT = getMargState(itT->first, rp, rcp);
            StatePtr spT = itT->first;
            StatePtr sp0 = it0->first;

            if( (*sr0 == *s0) && (*srT == *sT)) {
                // marginalized sum
                //edgeUpdateFileStream<<"\n\tmsgCtoP() s0: "<<*s0<<" sT: "<<(*sT)<<" sp0: "<<*sp0<<" spT: "<<*spT<<" sr0: "<<*sr0<<" srT: "<<*srT;

                LD bCurr = rp->B_0->getPr(sp0);
                LD zCurr = rp->Z->getPr(sp0);
                LD kCurr = rp->getConditionalBeliefPtr()->getPr(sp0, spT);
                LD pContrib = getPrntM(rp, sp0, spT);
                LD cContrib = getChldM(rp, sp0, spT);
                LD powFac = 1.0/((double)rp->cr);
                LD term = kCurr*pow((pContrib/cContrib), powFac)*bCurr/zCurr;
                Zs += term;
                //edgeUpdateFileStream<<"\n\tmsgCtoP() kCurr: "<<kCurr<<" pC: "<<pContrib<<" cC: "<<cContrib<<" powF: "<<powFac;
                //edgeUpdateFileStream<<" bCurr: "<<bCurr<<" zCurr: "<<zCurr<<" term: "<<term<<" Zs: "<<Zs;

            }
        }

    }


    double bChild = rcp->B_t_0->getPr(s0, sT);
    mNew = pow((1/bChild*Zs), rp->cr)*mOld;
    //edgeUpdateFileStream<<"\ngetMsgCtoP() mNew: "<<mNew<<" mOld: "<<mOld<<" bChild: "<<bChild<<" Zs: "<<Zs<<"\n";
    //edgeUpdateFileStream.close();
    return mNew;
}

pair<LD, LD> IsingRegionGraph::update_edges(LD cFactor) {
    iterNum++;
    LD mMax = 0.0;
    LD mMin = 1.0;
    foreach(nIter, adj_map) {
        RegionPtr rp = (RegionPtr) nIter->first;
        map<StatePtr, LD, StatePtrCmp> bMp = rp->getBeliefPtr()->getMp();
        foreach(edgeIter, nIter->second.second) {
            ASSERT(edgeIter->ends[0] == rp);
            pair<EdgePtr, EdgePtr> ep = findEdge(edgeIter->ends[0], edgeIter->ends[1]);
            ASSERT((ep.first != 0) && (ep.second != 0));
            EdgePtr cEdge = ep.second;
            //  cout<<"cEdge: "<<*cEdge<<"\n";
            foreach(stIter0, bMp) {

                StatePtr s0 = stIter0->first; // old state
                foreach(stIterT, bMp) {
                    StatePtr sT = stIterT->first;
                    LD mNew = getMsgCtoP(cEdge, s0, sT);
                    vector<INT> vals;
                    for(int i=0; i < (s0->num + sT->num); i++) {
                        int n;
                        if(i < s0->num) {
                            n = s0->getData(i);
                        } else
                            n = sT->getData(i-s0->num);
                        vals.push_back(n);
                    }
                    mNew = mNew*cFactor;
                    mMax = max(mNew, mMax);
                    mMin = min(mNew, mMin);
                    StatePtr sJoint = getStatePtr(vals);
                    LD mOld = cEdge->getMsg(sJoint);
                    cEdge->setMsg(sJoint, mNew);
                    /*
                        pair<EdgePtr, EdgePtr> ep2 = findEdge(edgeIter->ends[0], edgeIter->ends[1]);
                        LD mOld = ep2.second->getMsg(sJoint); 
                        ep2.second->setMsg(sJoint, mNew); 
                        LD nm = ep2.second->getMsg(sJoint);  
                     */
                    /*
                    cout<<"edge_p "<<cEdge->ends[0]->getNodeId()<<"\nedge_t"<<cEdge->ends[1]->getNodeId();
                    cout<<"\nsJoint: "<<*sJoint<<" mNew: "<<mNew<<" mOld: "<<mOld<<"\n";
                    getchar();
                    */
                }
            }
        }
    }
    return make_pair<LD, LD>(mMax, mMin);
}


bool IsingRegionGraph::bp_ppm(LD _th, LD _delta_t) {
    // initialize the conditional beliefs
    delta_t = _delta_t;
    theta = _th;
    foreach(nIter, adj_map) {

        RegionPtr rp = (RegionPtr) nIter->first;
        ConditionalBeliefPtr cptr = getCondBlfForRegn(rp);
        rp->setConditionalBelief(cptr);

    }

    double mres=0.0;
    foreach(nIter, adj_map) {
        RegionPtr rp = (RegionPtr) nIter->first;
        updateZ(rp);
        updateB_t_0(rp);
        cout<<"\n\tB_0: "<<*(rp->B_0)<<"\n\tB_t: "<<*(rp->B_t)<<"F: "<<(rp->getEnthalpy()-rp->getEntropy())<<"\n";
        LD res = updateB_t(rp);
        cout<<" maxChange: "<<res<<"\n";
        mres = max((double) mres, (double) res);
    }
    pair<LD, LD> mp = update_edges();
    cout<<"mMax: "<<mp.first<<" mMin: "<<mp.second<<" mres: "<<mres<<"\n";
    if(mres < 0.01)
        return true; // converged
    else
        return false;
}

/*
 * ch: whether to compare with true probability or bp probability
 * maxDev: the maximum deviation between two region belief and the corresponding b req (bp or true)
 * avgDiff: the maximum deviation between region averaged belief for node snf P_act
 */
bool IsingRegionGraph::ppm_iter(LD _th, LD _delta_t, LD* minMsg, LD* maxMsg, LD* maxDev, LD* avgDev, map<string, vector<LD> >* nextBlfMpPtr, LD lmt) {
    // initialize the conditional beliefs if already uninitialized
    if((delta_t != _delta_t) || (theta != _th)) {
        delta_t = _delta_t;
        theta = _th;
        foreach(nIter, adj_map) {

            RegionPtr rp = (RegionPtr) nIter->first;
            ConditionalBeliefPtr cptr = getCondBlfForRegn(rp);
            rp->setConditionalBelief(cptr);

        }
    }
    int fxdTimeCount=0;
FIXED_TIME_LOOP:      
    fxdTimeCount++; 
    LD maxChange = 0.0; 
    foreach(nIter, adj_map) {
        RegionPtr rp = (RegionPtr) nIter->first;
        updateZ(rp);
        LD currChange = updateB_t_0(rp);
        //   cout<<"\n\tB_0: "<<*(rp->B_0)<<"\n\tB_t: "<<*(rp->B_t)<<"F: "<<(rp->getEnthalpy()-rp->getEntropy())<<"\n";
        maxChange = max(maxChange, currChange); 
       
        
    }
    
  
    *maxDev = 0.0;
    *avgDev = 0.0;
    int numVar=0;
    IsingGraphPtr ig = graphPtr;
    foreach(iter, ig->adj_map) {
        if(!ig->isFactor(iter->first)) {
            numVar++;
            string vid = iter->first->getNodeId();
            LD sX=0.0;
            LD sXsqr = 0.0;
            LD count=0.0;
            vector<LD> cblfs;
            foreach(riter, this->adj_map) {
                string rid = riter->first->getNodeId();
                if(rid.find(vid) != string::npos) {
                    RegionPtr rpMy = (RegionPtr) riter->first;
                    LD rCurrPr = getMargPr(rpMy, vid);
                    sX += rCurrPr;
                    sXsqr += (rCurrPr)*(rCurrPr);
                    count += 1.0;
                    cblfs.push_back(rCurrPr);
                }

            }

            ASSERT(count != 0.0);
            LD muX = sX/count;
            foreach(ccIter, cblfs) { // compute max diff
                *maxDev = max(*maxDev, abs(muX-*ccIter));
            }
            LD varX = muX*muX - sXsqr/count;
            *avgDev += abs(varX);

            if(nextBlfMpPtr != 0) { // have to save the beliefs;
                map<string, vector<LD> >::iterator cccIter = nextBlfMpPtr->find(vid);
                if(cccIter == nextBlfMpPtr->end()) {
                    vector<LD> curr;
                    cccIter = nextBlfMpPtr->insert(make_pair<string, vector<LD> >(vid, curr)).first;
                }
                cccIter->second.push_back(muX);
            }
        }
    }
    *avgDev = (*avgDev)/((LD) numVar);

    pair<LD, LD> mp = update_edges();
    //pair<LD, LD> mp(1.0, 1.0);
    *maxMsg = mp.first;
    *minMsg = mp.second;
    bool res = renormalizeMessages(mp, lmt);
    
     cout<<"ppm_iter() maxChange: "<<maxChange<<" ";        
    if(maxChange < 0.001) {
        LD mres= 0.0; 
        
        foreach(nIter, adj_map) {
            RegionPtr rp = (RegionPtr) nIter->first;
            LD res = updateB_t(rp);
            mres = max(mres, res);
        }
        cout<<" converged_mres: "<<mres<<"\n"; 
    } else {
        cout<<" time_fixed_reiter "<<fxdTimeCount<<"\n";
        getchar(); 
        goto FIXED_TIME_LOOP;  
        
    }
    return res; // whether normalized
}

/*
 * If (node == factor) H = -\sum_\alpha b_\alpha * \log(b_\alpha)
 * else H = (d-1)*\sum_i b_i * \log(b_i) where d is the degree
 */
LD IsingGraph::getNodeEntropy(NodePtr np) {
    LD res = 0.0;
    if(isFactor(np) ) {

        vector<INT> vals = getVarVals();
        int Q = vals.size();
        vector<string> vars = split(np->getNodeId(), '_');
        vars.erase(vars.begin());
        int N = vars.size();
        map<StatePtr, LD, StatePtrCmp> bMp = np->getBeliefPtr()->getMp();
        BeliefPtr b_n = new Belief();
        LD Z=0.0;
        foreach(iter, bMp) {
            LD term = iter->second;
            for(int i=0; i < N; i++) {
                NodePtr npc = getNode(vars[i]); // find the child
                pair<EdgePtr, EdgePtr> ep = findEdge(np, npc);
                int cval = iter->first->getData(i); // the variable value
                StatePtr sp = getStatePtr(cval);
                LD cf = ep.second->getMsg(sp);
                term *= cf;
            }
            b_n->setPr(iter->first, term);
            Z += term;
        }
        BeliefPtr b_n2=new Belief();
        bool corrZ = false;
        if(Z > 0.0) {
            corrZ = true;
            foreach(iter, bMp) {
                StatePtr sp = iter->first;
                LD tm = b_n->getPr(sp)/Z;
                b_n2->setPr(sp, tm);
            }
        }

        if(corrZ) {
            foreach(iter, bMp) {

                StatePtr sp = iter->first;
                LD bc = b_n2->getPr(sp);
                if(bc > 0.0) {
                    res += -bc*log(bc);
                }
                //       cout<<"H: "<<np->getNodeId()<<" st: "<<*iter->first<<" bc: "<<bc<<" entropy(): "<<res<<"\n";
            }
        }
        delete(b_n);
        delete(b_n2);
    } else {
        map<StatePtr, LD, StatePtrCmp> bMp = np->getBeliefPtr()->getMp();

        foreach(iter, bMp) {
            if(iter->second > 0.0) {
                res += iter->second * log(iter->second);
            }
            //  cout<<"H: "<<np->getNodeId()<<" st: "<<*iter->first<<" bc: "<<iter->second<<" entropy(): "<<res<<"\n";
        }
        int d = getDegree(np);
        ASSERT(d > 0);
        res = (d-1)*res;

    }
    //getchar();
    return res;
}

/*
 * if(node == factor) U = -\sum_\alpha b_\alpha \log(f_\alpha)
 * else U = 0.0 // variable 
 */
LD IsingGraph::getNodeEnthalpy(NodePtr np) {
    if(isFactor(np)) {

        vector<INT> vals = getVarVals();
        int Q = vals.size();
        vector<string> vars = split(np->getNodeId(), '_');
        vars.erase(vars.begin());
        int N = vars.size();
        map<StatePtr, LD, StatePtrCmp> bMp = np->getBeliefPtr()->getMp();
        BeliefPtr b_n = new Belief();
        LD Z=0.0;
        foreach(iter, bMp) {
            LD term = iter->second;
            for(int i=0; i < N; i++) {
                NodePtr npc = getNode(vars[i]); // find the child
                pair<EdgePtr, EdgePtr> ep = findEdge(np, npc);
                int cval = iter->first->getData(i); // the variable value
                StatePtr sp = getStatePtr(cval);
                LD cf = ep.second->getMsg(sp);
                term *= cf;
            }
            b_n->setPr(iter->first, term);
            Z += term;
        }
        BeliefPtr b_n2 = new Belief();
        bool corrZ = false;
        if(Z > 0.0) {
            corrZ = true;
            foreach(iter, bMp) {
                StatePtr sp = iter->first;
                LD tm = b_n->getPr(sp)/Z;
                b_n2->setPr(sp, tm);
                //          cout<<"U: calcBa "<<np->getNodeId()<<" st: "<<*iter->first<<" tm: "<<tm<<" Z: "<<Z<<"\n";
            }
        }
        LD res = 0.0;
        if(corrZ) {
            foreach(iter, bMp) {
                StatePtr sp = iter->first;
                LD bc = b_n2->getPr(sp);
                LD lfc = log(iter->second);
                res += -bc*lfc;
                //        cout<<"U: getU "<<np->getNodeId()<<" st: "<<*iter->first<<" bc: "<<bc<<" lfc: "<<lfc<<"res: "<<res<<"\n";
            }
        }
        delete(b_n);
        delete(b_n2);
        //   getchar();
        return res;
    } else
        return 0.0;
}

LD IsingGraph::getF() {
    // if(!dirtyF) return Fapprox;
    //map<NodePtr, BidirEdgeList, NodeCmp>::iterator mIt =  adj_map.find(np);
    LD F=0.0;
    foreach(iter, adj_map) {
        F += getNodeEnthalpy(iter->first) - getNodeEntropy(iter->first);
    }
    Fapprox = F;
    dirtyF = false;
    return F;
}


LD IsingRegionGraph::getF() {
    //map<NodePtr, BidirEdgeList, NodeCmp>::iterator mIt =  adj_map.find(np);
    LD F=0.0;
    foreach(iter, adj_map) {
        RegionPtr rp = (RegionPtr) iter->first;
        F += rp->cr*(rp->getEnthalpy() - rp->getEntropy());

    }
    return F;
}


// find the prbability of node states
LD IsingRegionGraph::getMargPr(const RegionPtr& rp, string var) {
    // variable is present
    ASSERT(rp->getNodeId().find(var) != string::npos);
    int idx=0;
    vector<string> varids = rp->getVarIds();
    while(varids[idx] != var)
        idx++;
    LD Z=0.0, Zt = 0.0;
    map<StatePtr, LD, StatePtrCmp> bMp = rp->B_0->getMp();
    // cout<<"getMargPr() varids: "; foreach(iter, varids) cout<<" "<<*iter; cout<<" var: "<<var<<"\n";
    foreach(iter, bMp) {
        if(iter->first->getData(idx) == 1) {

            Z += iter->second;
            //       cout<<"\tst: "<<*(iter->first)<<" val: "<<iter->second<<" Z: "<<Z<<"\n";
        }
        Zt += iter->second;
    }

    return Z/Zt;
}

void IsingGraph::bp_synch(bool ch) {
    bool redo;
    int count=0;
    LD F_old = 0.0;
    LD F_new;
LOOP:
    if(ch)
        getTrueZ();
    count++;
    cout<<"bp_synch(iteratn: "<<count<<")";
    F_new = getF();
    cout<<" F_approx: "<<F_new;

    if(ch) {
        cout<<" F_true: "<<-log(getTrueZ());
        LD maxDiff=0.0;
        vector<string> varIds = getVarIds();
        vector<INT> varVals = getVarVals();
        foreach(iter, varIds) {
            string id = *iter;
            //   cout<<id;
            NodePtr np = getNode(id);
            ASSERT(np != 0);
            //   cout<<"\n"<<np->getNodeId()<<"\t";
            foreach(viter, varVals) {
                StatePtr sp = getStatePtr(*viter);

                LD p_bp = np->getBeliefPtr()->getPr(sp);
                //     printf("v: %d p: %g ", *viter, p_bp);
                pair<string, int> cp = make_pair<string, int>(*iter, *viter);
                vector<pair<string,int> > vp;
                vp.push_back(cp);
                LD p_act = getMarginalProb(vp);
                //       printf("p_t: %g ", p_act);
                maxDiff = max(maxDiff, abs(p_bp - p_act));
            }
        }
        cout<<" maxDiff: "<<maxDiff;
    }
    cout<<"\n";
    redo = bp_iter(ch);
    if(abs(F_old-F_new) > 0.0) {
        F_old = F_new;
        goto LOOP;
    }
}

bool IsingRegionGraph::renormalizeMessages(pair<LD, LD> mp, LD lmt) {
    LD maxM = mp.first;
    LD minM = mp.second;
    bool res = false;
    /*
    if(isnan(maxM) || (minM == 0.0) || isinf(maxM) ) {
        res = true;
        foreach(iter, adj_map) {
            foreach(eIter, iter->second.second) {
                foreach(mIter, eIter->msgMp) {
                    mIter->second = 1.0;
                }
            }
        }
        maxM = minM = 1.0;
    }
    */ 
    LD lmm = log(minM);

    if(maxM/minM > lmt) {
        res = true;
        cout<<" maxM: "<<maxM<<" minM: "<<minM <<"maxM/minM: "<<maxM/minM<<" lmt: "<<lmt<<"\n"; 
        foreach(iter, adj_map) {
            foreach(eIter, iter->second.second) {
                foreach(mIter, eIter->msgMp) {
                    mIter->second = log(mIter->second)-lmm + 1.0;
                    // mIter->second = 1.0;
                }
            }
        }
    }

    return res;
}



