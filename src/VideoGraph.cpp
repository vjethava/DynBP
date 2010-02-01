#include "VideoGraph.h"
#include "VideoReader.h"
#include <cassert>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <set>
#include <algorithm>
using namespace std;

bool VideoGraph::isValid(int x, int y, int xi, int yi, int xl, int yl) {
    if(xl == -1)
        xl = N1;
    if(yl == -1)
        yl = N2;
    return ((x>=xi) && (y>=yi) && (x < (xi+xl)) && (y < (yi+yl)));
}

VideoGraph::VideoGraph(VideoReader* _reader, bool _isGray, int numL, int _rx, int _ry, int _cx, int _cy) : GraphPPM() {
    FILE* nodeFile = fopen("node_list", "w");
    FILE* gFile = fopen("GNFP.log","w");
    fprintf(gFile, "#########################################################################################################################\n");
    fprintf(gFile, "# frame iter trueFrameDiff [minM maxM <ParDiff> <Cblf> blfC_max <blfC_max>]\n");
    fprintf(gFile, "#########################################################################################################################\n");

    fclose(gFile);
    fclose(nodeFile);
    reader = _reader;
    output_frame=NULL;
    local_mask = NULL;
    inpainted_frame = NULL;
    isGray = _isGray;
    numLvls =  numL;
    reader->setGrayScale(isGray);
    reader->setNumLvls(numLvls);
    reader->setGraph(this);
    fi(0, numLvls) { // push the level
        levels.push_back(i*256/numLvls);
    }
    hits = misses=0;
    rx= _rx;
    ry = _ry;
    cx= _cx;
    cy=_cy;
    N1 = reader->getN1();
    N2 = reader->getN2();
    // generates the appropriate graph structure
    genSubRegions();
}

VideoGraph::~VideoGraph() {
    // delete CPT
    //fclose(nodeFile);
    foreach(iter, blfMp) {
        delete(iter->second);
    }
}


void VideoGraph::updateCondBlfMp(IplImage* last_frame, IplImage* curr_frame, IplImage* mask_frame) {
    int mmsize=0;
    foreach(iter, adj_map) {
        VideoNode* vn = (VideoNode*) iter->first;
        string nid = vn->getNodeId();
        if(vn->isMasked == false) {
            vector<int> vals = VideoNode::getIntFromId(nid);
            vector<INT> cblf = reader->getDataFromFrame(last_frame, vals[0], vals[1], vals[2], vals[3]);
            vector<INT> nblf = reader->getDataFromFrame(curr_frame, vals[0], vals[1], vals[2], vals[3]);
            int msize = addCondBlf(cblf, nblf);
            mmsize = max(msize, mmsize);
        }
    }
    printf("updateCondBlfMp() hits: %d misses: %d state_mp_size: %d mmsize: %d\n", hits, misses, sMp.size(), mmsize);
}

int VideoGraph::addCondBlf(vector<INT> pastVals, vector<INT> nextVals) {
    int N = pastVals.size();
    assert(nextVals.size() == N);
    ConditionalBeliefPtr cptr = blfMp.find(N)->second;
    StatePtr cst = getStatePtr(pastVals);
    StatePtr nst = getStatePtr(nextVals);
    LD blf = cptr->getPr(cst, nst);
    BeliefPtr bptr =  cptr->getBeliefForState(cst);
    map<StatePtr, LD, StatePtrCmp> mp;
    int oldSize;
    bool spCase = false;
    if(bptr != 0) {
        mp = bptr->getMp();
        /* if(!((*cst) == (*nst))) {
             cout<<"cst: "<<(*cst)<<"\nnst: "<<(*nst)<<"\n";
             cout<<"old state map:\n"; 
             foreach(iter, mp) {
                 cout<<"\t"<<*(iter->first)<<" p: "<<iter->second<<"\n";
             }
         }*/
        spCase = true;
        oldSize = mp.size();
    } else {

        oldSize = 0;
    }
    bool wasHit = true;
    if(blf == BELIEF_NOT_FOUND) {
        cptr->setPr(cst, nst, 1.0);
        misses++;
        wasHit = false;
    } else {
        cptr->setPr(cst, nst, blf+1.0);
        hits++;
    }
    bptr =  cptr->getBeliefForState(cst);
    mp = bptr->getMp();
    int newSize = mp.size();

    LD nblf = cptr->getPr(cst, nst);
    /*
    if((spCase)) {
        
        cout<<" osize: "<<oldSize<<" newSize: "<<newSize<<" wasHit: "<<wasHit<<"\n";
    //    sleep(2);          
} */
    return oldSize;
}

void VideoGraph::setNodeCondBlf(VideoNode* vn) {
    if(vn->getCPTR() == 0) {

        int N = vn->getNodeNumVar();
        if(!isGray) {
            N *= 3; // rgb components viewed in a single row
        }
        if(blfMp.find(N) == blfMp.end()) { // insert a conditional belief
            ConditionalBeliefPtr cptr = new ConditionalBelief();
            blfMp.insert(make_pair<int, ConditionalBeliefPtr>(N, cptr));
        }
        ConditionalBeliefPtr cptr = blfMp.find(N)->second;
        vn->setCPTR(cptr);
    }
}

/// this function adds the requisite graph structure to the VideoGraph
void VideoGraph::genSubRegions() {
    // add the top-level regions
    //FILE* rgnFile = fopen("regionCreation.txt", "w");

    for(int i=0; i < N1; i+= cx) {
        for(int j=0; j < N2; j+=cy) {
            int wx = min(rx , N1-i);
            int wy = min(ry , N2-j);
            if((wx > 0) && (wy > 0)) {
                stringstream ss("");
                ss<<i<<"_"<<j<<"_"<<wx<<"_"<<wy;
                VideoNode* vn =  new VideoNode(ss.str());
                vn->setCR(1);
                setNodeCondBlf(vn);

                // top level nodes
                // cout<<" top level node: "<<*vn<<"\n"; // sleep(1);
                addNode(vn);
                //  fprintf(rgnFile, "%s\n", vn->getNodeId().c_str());

            }

        }
    }

    //fprintf(rgnFile, "L2_HORIZ\n"); //, vn->getNodeId.c_str());

    // add the next horizontal overlap region as top-level
    for(int i=0; i < N1; i+= cx) {
        for(int j=0; j < N2; j+= cy) {
            int wx = min(rx , N1-i);
            int wy = min(ry , N2-j);
            int xn = i+cx;
            int wxn = min(rx, N1-xn);
            if( (wxn > 0) && (wy > 0) && (wx > 0)) { // overlap exists
                stringstream ss("");
                int wc = i+rx-xn;
                ss<<xn<<"_"<<j<<"_"<<wc<<"_"<<wy;
                VideoNode* vn2 = new VideoNode(ss.str());
                vn2->setCR(-1);
                stringstream ssA(""), ssB("");
                ssA<<i<<"_"<<j<<"_"<<wx<<"_"<<wy;
                VideoNode* vnA = (VideoNode*) getNode(ssA.str());
                //       int nxw = min(rx, N1 - xn);
                ssB<<xn<<"_"<<j<<"_"<<wxn<<"_"<<wy;
                VideoNode* vnB = (VideoNode*) getNode(ssB.str());
                if( (vnA != 0) && (vnB != 0)
                        && VideoNode::isSubset(vn2->getNodeId(), vnA->getNodeId())
                        && VideoNode::isSubset(vn2->getNodeId(), vnA->getNodeId())
                  ) {
                    addNode(vn2);
                    addEdge(vnA, vn2, true);
                    addEdge(vnB, vn2, true);
                } else {
                    delete(vn2);
                }
            }
        }
    }

    //fprintf(rgnFile, "L2_VERT\n"); //, vn->getNodeId.c_str());

    // add the next vertical overlap region as top-level
    for(int i=0; i < N1; i+= cx) {
        for(int j=0; j < N2; j+= cy) {
            int wx = min(rx , N1-i);
            int wy = min(ry , N2-j);
            int yn = j+cy;
            int wyn = min(ry, N2-yn);
            if( (wyn > 0) && (wx > 0) && (wy > 0)) { // overlap exists
                stringstream ss("");
                int wc = j+ry-yn;
                ss<<i<<"_"<<yn<<"_"<<wx<<"_"<<wc;
                VideoNode* vn2 = new VideoNode(ss.str());
                vn2->setCR(-1);
                stringstream ssA(""), ssB("");
                ssA<<i<<"_"<<j<<"_"<<wx<<"_"<<wy;
                VideoNode* vnA = (VideoNode*) getNode(ssA.str());
                // int nyw = min(ry, N2 - yn);
                ssB<<i<<"_"<<yn<<"_"<<wx<<"_"<<wyn;
                VideoNode* vnB = (VideoNode*) getNode(ssB.str());
                if( (vnA != 0) && (vnB != 0)
                        && VideoNode::isSubset(vn2->getNodeId(), vnA->getNodeId())
                        && VideoNode::isSubset(vn2->getNodeId(), vnA->getNodeId())
                  ) {
                    addNode(vn2);
                    addEdge(vnA, vn2, true);
                    addEdge(vnB, vn2, true);
                } else {
                    delete(vn2);
                }

            }
        }
    }


    //    //fprintf(rgnFile, "\nL3\n");
    //    // add the third level regions
    //    for(int i=0; i < N1; i += cx) {
    //        for(int j=0; j < N2; j += cy) {
    //            if( isValid(i-cx, j-cy) && isValid(i+rx, j+ry) ) { // central blocks
    //                VideoNode* vn = new  VideoNode(VideoNode::getIdFromInt(i, j, rx-cx, ry-cy));
    //                vn->setCR(1);
    //                string sA = VideoNode::getIdFromInt(i-cx, j, rx, ry-cy);
    //                VideoNode* vnA = (VideoNode*) getNode(sA);
    //                VideoNode* vnB = (VideoNode*) getNode(VideoNode::getIdFromInt(i, j, rx-cx, ry));
    //                VideoNode* vnC = (VideoNode*) getNode(VideoNode::getIdFromInt(i, j, rx, ry-cy));
    //                VideoNode* vnD = (VideoNode*) getNode(VideoNode::getIdFromInt(i, j-cy, rx-cx, ry));
    //                setNodeCondBlf(vn);
    //                BeliefPtr bptr = new Belief(true);
    //                vn->setBelief(bptr);
    //
    //                addNode(vn);
    //                //fprintf(rgnFile, "%s\n", vn->getNodeId().c_str());
    //                bool bA=false, bB=false, bC=false, bD=false;
    //                if(vnA != 0) {
    //                    addEdge(vnA, vn, true);
    //                }
    //                if(vnB != 0) {
    //                    addEdge(vnB, vn,true );
    //                }
    //                if(vnC != 0) {
    //                    addEdge(vnC, vn, true);
    //                }
    //                if(vnD != 0) {
    //                    addEdge(vnD, vn, true);
    //                }
    //            }
    //        }
    //    }
    //    //fclose(rgnFile);
    //
}

int VideoGraph::updateNodeBlf(VideoNode* vn, IplImage* frame) {
    string nid = vn->getNodeId();
    vector<int> vals = VideoNode::getIntFromId(nid);
    vector<INT> cblf = reader->getDataFromFrame(frame, vals[0], vals[1], vals[2], vals[3]);
    StatePtr sptr = getStatePtr(cblf);
    BeliefPtr blfPtr = vn->getGlobalBlfPtr();
    if(blfPtr == 0) {
        int key = vn->getNodeNumVar();
        BeliefPtr res = getOrigBlfFromMp(key);
        vn->setGlobalBlfPtr(res);
    }

    if(isGray)
        ASSERT(sptr->num == vn->getNodeNumVar());
    else
        ASSERT(sptr->num == (vn->getNodeNumVar()*3));
    vn->getGlobalBlfPtr()->setPr(sptr);
    map<StatePtr, LD, StatePtrCmp> mp = vn->getGlobalBlfPtr()->getMp();
    return mp.size();
}

LD VideoGraph::getNodeMaskPct(VideoNode* vn, IplImage* mask_frame) {
    LD result = 0.0;
    LD count = 0.0;
    vector<int> L = VideoNode::getIntFromId(vn->getNodeId());
    BwImage maskImage(mask_frame);
    for(int i=L[0]; i < (L[0]+L[2]); i++) {
        for(int j=L[1]; j < (L[1]+L[3]); j++) {
            if(maskImage[i][j] != 0) {
                count = count + 1.0;
            }
        }
    }
    result = count/vn->getNodeNumVar();
    /* if(count != 0) {
         cout<<" node: "<<vn->getNodeId()<<" count: "<<count<<" res: "<<result<<"\n"; getchar();
     }*/
    return result;
}


void VideoGraph::updateMaskStatus(IplImage* mask_frame) {
    FPRINTF(stderr,"updateMaskStatus()\n");
    foreach(iter, adj_map) {
        VideoNode* vn = (VideoNode*) iter->first;
        if(getNodeMaskPct(vn, mask_frame) < 1.0e-9) {
            vn->isMasked = false;
        } else {
            vn->isMasked = true;
        }
    }
}
int VideoGraph::updateAllNodes(IplImage* frame, IplImage* mask_frame) {
    int res=0;
    Timer timer;
    FPRINTF(stderr, "updateAllNodes()\n");
    foreach(iter, adj_map) {
        VideoNode* vn = (VideoNode*) iter->first;
        // if(VideoNode::intersection(vn->getNodeId(), "25_125_65_50") == "") {
        timer.tick();
        if(vn->isMasked == false) {
            int n = updateNodeBlf(vn, frame);
            res += n;

            //      cout<<" node blf: "<<vn->getNodeId()<<" size: "<<n; getchar();
        }

    }
    cout<<" done\n";
    return res;
}

BeliefPtr VideoGraph::getApproxBlf(VideoNode* vn, StatePtr sOld) {
    Timer timer;

    cout<<"getApproxBlf() node: "<<vn->getNodeId(); // <<"state: "<<*sOld<<"\n";
    map<StatePtr, LD, StatePtrCmp> bMp = pOrigMp.find(vn->getNodeNumVar())->second->getMp();
    vector<pair<StatePtr, LD> > states;
    LD avgDist=0.0, numDist=0.0;
    vector<LD> dists;
    foreach(iter, bMp) {
        //  cout<<"sOld: "<<*sOld<<" size: "<<sOld->num<<endl;
        //  cout<<"iter: "<<*iter->first<<" size: "<<iter->first->num<<endl;
        //   getchar();
        ASSERT(iter->first->num == sOld->num);
        states.push_back(make_pair<StatePtr, LD>(iter->first, iter->second));
        LD currD = State::SSD(*sOld, *iter->first);
        avgDist += currD;
        numDist += 1.0;
        dists.push_back(currD);
    }
    sort(dists.begin(), dists.end());
    int idx = min(10, (int)dists.size()-1);
    LD criteriaDist = dists[idx];
    ASSERT(numDist > 0.0);

    avgDist = avgDist/numDist;
    //ANN_Cmp cmp(this, vn, NULL, sOld);
    //sort(states.begin(), states.end(), cmp);
    LD totFreq = vn->getGlobalBlfPtr()->getTotalFreq();
    // cout<<" node: "<<vn->getNodeId()<<" state: "<<*sOld<<" totalFreq: "<<totFreq<<"\n";
    BeliefPtr bptr = new Belief();
    BeliefPtr tempPtr = new Belief();
    LD sumT = 0.0;
    foreach(iter, states) {
        ASSERT(sOld->num == iter->first->num);
        LD diff = State::SSD(*sOld, *(iter->first));

        LD distFac = 1.0/(1.0+diff);
        LD pSt = (iter->second);
        BeliefPtr currBlfPtr = vn->getCPTR()->getBeliefForState(iter->first);

        if((currBlfPtr != 0) && (diff <= criteriaDist)) {
            map<StatePtr, LD, StatePtrCmp> currMp = currBlfPtr->getMp();
            foreach(mpIter, currMp) {
                timer.tick();
                StatePtr sNew = mpIter->first;
                LD pAdd = mpIter->second;
                LD pOld = max(0.0, tempPtr->getPr(sNew));
                pAdd = pAdd*distFac*pSt/totFreq;
                sumT = sumT + pAdd;
                LD pNew = pOld + pAdd;
                tempPtr->setPr(sNew, pNew);
                //            printf("\ndFac: %f pSt: %f pOld: %f pNew: %f\n", distFac, pSt, pOld, pNew);
                //            cout<<"\t"<<*(iter->first)<<" p: "<<iter->second<<" diff: "<<diff<<"\n";
                //            getchar();
            }
        }

    }
    // cout<<"sumT: "<<sumT<<" allocPr: "<<tempPtr->getAllocPr()<<"\n"; getchar();
    map<StatePtr, LD, StatePtrCmp> sMp = tempPtr->getMp();
    foreach(iter, sMp) {
        bptr->setPr(iter->first, ((iter->second)/sumT));
    }
    cout<<" done\n";
    delete(tempPtr);
    return bptr;
}

bool VideoGraph::isClose(StatePtr s1, StatePtr s2, int NL, LD distS) {
    ASSERT(s1->num == s2->num);
    LD div = ((LD) 256)/((LD) NL);
    LD ssd = State::SSD(*s1, *s2);
    LD normD = ssd/div*sqrt(s1->num);
    if(normD < distS)
        return true;
    else
        return false;
}



/*
vector<pair<StatePtr, LD> > VideoGraph::getNextStateCands(VideoNode* vn, IplImage* frame, LD candT) {
    Timer timer; 
    vector<pair<StatePtr, LD> > result, tempVec; 
    StatePtr sptr = readNodeFromFrame(vn, frame);
    int snum = sptr->num;
    ConditionalBeliefPtr cptr = vn->getCPTR(); 
    map<StatePtr, Belief, StatePtrCmp> *cbMp = cptr->getCbMpPtr();
    BeliefPtr cblf = cptr->getBeliefForState(sptr); 
    ASSERT(cblf != 0);
    map<StatePtr, LD, StatePtrCmp> cMp = cblf->getMp(); 
    cout<<"vn: "<<vn->getNodeId()<<" cbMp.size(): "<<cbMp->size()<<" cMp.size(): "<<cMp.size()<<"\n";
    getchar();
    StatePtr sptr =   
    foreach(iter, cMp) {
        
    }
    foreach(iter, cMp) { 
        tempVec.push_back(make_pair<StatePtr, LD>(iter->first, iter->second)); 
    }
    foreach(iter, *cbMp) {
        
    
    
    }
    return result;   
}
*/

/// ultimately has to look at gaussian mixture model of the K-NN states
vector<pair<StatePtr, LD> > VideoGraph::getNextStateCands(VideoNode* vn, IplImage* frame, LD candT) {
    Timer timer;
    // simplification 1 - look at just next hits
    // map<StatePtr, LD, StatePtrCmp > candMp;
    vector<pair<StatePtr, LD> > candMp;
    StatePtr sptr = readNodeFromFrame(vn, frame);
    int N = vn->getNodeNumVar();
    if(!isGray)
        N*=3;
    ASSERT(N == sptr->num);

    ASSERT(blfMp.find(N) != blfMp.end());
    ConditionalBeliefPtr cptr = vn->getCPTR();
    BeliefPtr cblf = cptr->getBeliefForState(sptr);
    if(cblf == 0) {
        cblf = getApproxBlf(vn, sptr);

        map<StatePtr, LD, StatePtrCmp> cMp = cblf->getMp();
        cout<<"vn: "<<vn->getNodeId()<<" cMp.size(): "<<cMp.size()<<endl;
        foreach(mpIter, cMp) {
            cptr->setPr(sptr, mpIter->first, mpIter->second);
        }
        delete(cblf);
        cblf =  cptr->getBeliefForState(sptr);
        ASSERT(cblf != 0);
    }
    int i=0 ;
    map<StatePtr, LD, StatePtrCmp> mp = cblf->getMp();
    LD ap = cblf->getAllocPr();
    int mpSize = mp.size();
    //if((mpSize > 1) ) {
    //    printf("\ngetNextStateCands(%s): ", vn->getNodeId().c_str());
    //    cout<<"\norig: "<<*sptr;
    //    cout<<" mpSize: "<<mp.size()<<" ap: "<<ap<<"\n"<<flush;
    //}
    LD pT = 0.1, pCov = 0.0;
    int numa = 0;
    if(mpSize > 1) {
        LD pT = 1.0/((LD) mpSize);
    }
    //  vector<StatePtr> candStates;
    foreach(iter, mp) {

        ASSERT(cblf->getAllocPr() > 0);
        LD p = ((LD) iter->second)/((LD) ap);

        //if(mpSize > 1) {
        //  cout<<"\ti: "<<i++<<" "<<flush;
        //  cout<<*(iter->first)<<" p: "<<p<<"\n";
        // }


        //cout<<"op.size(): "<<op.size()<<" "; getchar();
        //StatePtr spn = getStatePtr(op);
        StatePtr spn = iter->first;
        ASSERT(spn != 0);
        int maxNum = ((isGray)?(rx*ry):(3*rx*ry));
        ASSERT((spn->num > 0) && (spn->num <= (maxNum) ));
        // foreach(iter3, op) { cout<<" "<<*iter3<<flush; }

        //candMp.insert(make_pair<StatePtr, LD>(spn , p));
        timer.tick();
        if(p > pT*candT) {
            //       candStates.push_back(spn);
            candMp.push_back(make_pair<StatePtr, LD>(spn , p));
            numa++;
            pCov += p;
        }
    }
    int div = 256/numLvls;
    if(candMp.size() <= 1) {
        StatePtr cand0 = candMp[0].first;
        map<StatePtr, LD, StatePtrCmp> myMp2 = vn->getBeliefPtr()->getMp();
        vector<pair<StatePtr, LD> > prVec;
        foreach(iter, myMp2) {

            int N = iter->first->num;
            LD diff = State::SSD(*cand0, *(iter->first));

            LD val = diff/((LD) div);
            //cout<<" addState: N: "<<N<<"\ns0:"<<*cand0<<"*\niter:"<<*iter->first<<" NL: "<<numLvls<<" val: "<<val<<"\n";
            //getchar();
            if(diff/div < 1.0) {
                candMp.push_back(make_pair<StatePtr, LD>(iter->first, iter->second));
            }
        }
    }

    if(pCov < 1.0) {
        //FPRINTF(stderr, "getNextStateCand(node: %s) numCand: %d numAcc: %d pCov: %f\n",
        //            vn->getNodeId().c_str(), mpSize, numa, pCov);
    }
    return candMp;
    //   return candMp;
}


/// read up the node and get equivalent state pointer
StatePtr VideoGraph::readNodeFromFrame(VideoNode* vn, IplImage* frame) {
    string nid = vn->getNodeId();
    vector<int> vals = VideoNode::getIntFromId(nid);
    vector<INT> cblf = reader->getDataFromFrame(frame, vals[0], vals[1], vals[2], vals[3]);
    StatePtr sptr = getStatePtr(cblf);
    return sptr;
}

vector<pair<VideoNode*, LD> > VideoGraph::getIncompleteList(IplImage* mask_frame) {
    vector<pair<VideoNode*, LD> > incompleteList;
    cout<<"getIncompleteList() \n";
    foreach(iter, adj_map) {
        VideoNode* vn = (VideoNode*) iter->first;
        LD pct = getNodeMaskPct(vn, mask_frame);
        if(pct > 0.0) {
            incompleteList.push_back(make_pair<VideoNode*, LD>(vn, pct));
        }
    }
    sort(incompleteList.begin(), incompleteList.end(), PRCMP2<VideoNode*, LD>());
    // foreach(iter, incompleteList) {
    //   cout<<" node: "<<iter->first->getNodeId()<<" pct: "<<iter->second<<endl;
    //    getchar();
    //}
    return incompleteList;
}

void VideoGraph::PPM_step2(IplImage* frame, LD candT) {
    Timer timer;
    parAugmentFound = parAugmentNotFound = 0;
    foreach(iter,  adj_map) { // for all the nodes - initialize the candidate maps
        VideoNode* vn = (VideoNode*) iter->first;
        vector<pair<StatePtr, LD> > candMp = getNextStateCands(vn, frame, candT);
        vn->setNextStateCandMp(candMp);
        vector<pair<StatePtr, LD> > pMp = vn->getNextStateCandMp();
        /*if(pMp.size() > 1) {
            cout<<"node: "<<vn->getNodeId()<<"\n";
            foreach(iter, pMp) {
                cout<<"\t"<<*(iter->first)<<" "<<iter->second<<"\n";
            }
          
    }*/
    }
    int totAugmentSize=0;
    int childAugments=0;
    if(candT > 0.0) {
        // parent augmentation

        foreach(iter, adj_map) {
            VideoNode* vn = (VideoNode*) iter->first;
            StatePtr sOld = readNodeFromFrame(vn, frame);
            updateParentAugments(vn, sOld);
            totAugmentSize += vn->augmentsMp.size();
            timer.tick();
        }

        // child augmentation

        foreach(iter, adj_map) {
            VideoNode* vn = (VideoNode*) iter->first;
            StatePtr sOld = readNodeFromFrame(vn, frame);
            int aSizeOld = vn->augmentsMp.size();
            updateChildAugments(vn, sOld);
            childAugments += vn->augmentsMp.size()-aSizeOld;
            timer.tick();
        }
    }
    foreach(iter, adj_map) {
        VideoNode* vn = (VideoNode*) iter->first;

        StatePtr sOld = readNodeFromFrame(vn, frame);

        ASSERT(sOld != 0);
        vn->finalizeCands(sOld);
        vn->setOldState(sOld);
        timer.tick();
        //cout<<" mode: "<<vn->getNodeId()<<" cands: "<<vn->nextStateMp.size()<<"\n"; getchar();
    }

    FPRINTF(stderr, "PPM_step2() parAugmentNotFound: %d parAugments: %d childAugments: %d\n", parAugmentNotFound, totAugmentSize, childAugments);
}


bool VideoGraph::getAbsolutePixel(string nid, int i, int j, StatePtr sptr, vector<INT>* result) {
    vector<int> L = VideoNode::getIntFromId(nid);
    if(!isValid(i, j, L[0], L[1], L[2], L[3]))
        return false;
    result->clear();
    int idx;
    if(isGray) {
        idx = (i-L[0])*L[3] + (j-L[1]);
        INT val = sptr->getData(idx);
        result->push_back(val);
    } else {
        idx = 3*((i-L[0])*L[3] + (j-L[1]));
        for(int i=0; i < 3; i++) {
            INT val = sptr->getData(idx + i);
            result->push_back(val);
        }
    }
    return true;
}



StatePtr VideoGraph::getChildStatePtr(string parId, string childId, StatePtr parState) {
    vector<int> LC = VideoNode::getIntFromId(childId);
    vector<INT> result, current;
    //cout<<" par: "<<parId<<" child: "<<childId<<" parSt: "<<*parState<<"\n";
    for(int i=LC[0]; i < LC[0]+LC[2]; i++) {
        for(int j=LC[1]; j < LC[1]+LC[3]; j++) {
            bool res = getAbsolutePixel(parId, i, j, parState, &current);
            if(res != true) {
                cout<<"getChildStatePtr() par: "<<parId<<" child: "<<childId<<" parSt: "<<*parState<<"\n";
                cout<<"i: "<<i<<" j: "<<j;
                getchar();
            }
            ASSERT(res == true);
            int idx=0;
            foreach(iter, current) {
                //      cout<<"\ti: "<<i<<" j: "<<j<<" color: "<<idx++<<" val: "<<((int)*iter)<<"\n";
                result.push_back(*iter);
            }
        }
    }
    StatePtr cptr = getStatePtr(result);
    return cptr;
}

void VideoGraph::updateChildAugments(VideoNode* vn, StatePtr oldState) {

    map<NodePtr, BidirEdgeList, NodeCmp>::iterator mIter = adj_map.find(vn);
    ASSERT(mIter != adj_map.end());
    vector<MyEdge> outEdges = mIter->second.first;
    string nodeId = vn->getNodeId();
    LD pFac =  vn->getCPTR()->getAllocPr(oldState);
    map<StatePtr, LD, StatePtrCmp> localMp;
    vector<pair<StatePtr, LD> > selfCands = vn->getNextStateCandMp();
    vector<pair<string, StatePtr> > unaccountedMinors;
    foreach(edgeIter, outEdges) {
        VideoNode* vc = (VideoNode*) edgeIter->ends[1];
        string childId = vc->getNodeId();

        localMp.clear();
        foreach(iter, selfCands) {
            StatePtr currChildState = getChildStatePtr(nodeId, childId, iter->first);
            LD pCurr = iter->second;
            localMp.insert(make_pair<StatePtr, LD>(currChildState, pCurr));
        }

        vector<pair<StatePtr, LD> > childCands = vc->getNextStateCandMp();
        foreach(stMpIter, childCands) {
            StatePtr currChildSt = stMpIter->first;
            LD pCurr = stMpIter->second;
            if(localMp.find(currChildSt) == localMp.end()) {
                unaccountedMinors.push_back(pair<string, StatePtr>(childId, currChildSt));
            }
        }
    }
    int N0 = unaccountedMinors.size();
    // need child augmentations
    if(unaccountedMinors.size() > 0) {
        BeliefPtr cblf = vn->getCPTR()->getBeliefForState(oldState);
        ASSERT(cblf != 0);
        int i=0 ;
        map<StatePtr, LD, StatePtrCmp> mp = cblf->getMp();
        vector<pair<StatePtr, LD> > vec1;
        foreach(iter, selfCands) {
            mp.erase(iter->first);
        }
        foreach(iter, mp) {
            vec1.push_back(make_pair<StatePtr, LD>(iter->first, iter->second));
        }
        //printf("vec size: %d umSize: %d\n", vec1.size(), N0);
        int N = N0;
        sort(vec1.begin(), vec1.end(), SPLD_Cmp());
        vector<pair<StatePtr, LD> >::iterator iter = vec1.begin();
        bool oper = false;
        while( (N > 0) && (iter != vec1.end()) ) {
            //   cout<<*iter->first<<" "<<iter->second<<"\n";
            bool add=false;
            vector<vector<pair<string, StatePtr> >::iterator > toDeleteIters;
            foreach(umIter, unaccountedMinors) {

                string cId = umIter->first;

                StatePtr childState = umIter->second;
                StatePtr currState = getChildStatePtr(nodeId, cId, iter->first);
                // cout<<" child: "<<cId<<" st: "<<*childState<<" margSt: "<<*currState<<"\n";
                if(*currState == *childState) {
                    oper = true;
                    add = true;
                    toDeleteIters.push_back(umIter);
                }
            }
            foreach(tdIter, toDeleteIters) {
                unaccountedMinors.erase(*tdIter);
                N--;
            }
            if(add) {
                vn->addAugmentCand(iter->first, iter->second);
                i++;
            }
            iter++;
        }
        //cout<<"node: "<<nodeId<<" N: "<<N0<<" added: "<<i<<"\n";
    }
}

void VideoGraph::updateParentAugments(VideoNode* vn, StatePtr oldState) {
    map<NodePtr, BidirEdgeList, NodeCmp>::iterator mIter = adj_map.find(vn);
    ASSERT(mIter != adj_map.end());
    vector<MyEdge> inEdges = mIter->second.second;
    string nodeId = vn->getNodeId();
    LD pFac =  vn->getCPTR()->getAllocPr(oldState);
    foreach(edgeIter, inEdges) {
        VideoNode* vp = (VideoNode*) edgeIter->ends[0];
        string parId = vp->getNodeId();
        vector<pair<StatePtr, LD> > parCands = vp->getNextStateCandMp();
        foreach(stMpIter, parCands) {
            StatePtr parState = stMpIter->first;
            StatePtr cStPtr = getChildStatePtr(parId, nodeId, parState);
            LD pChild = vn->getCPTR()->getPr(oldState, cStPtr);
            if(pChild == BELIEF_NOT_FOUND) {
                /*cout<<"par: "<<parId<<" st: "<<*(stMpIter->first)<<"\n";
                  cout<<"chd: "<<nodeId<<" old_st: "<<*oldState<<" pFac: "<<pFac<<"\n";
                  cout<<"chd: "<<nodeId<<" new_st: "<<*cStPtr<<"\n" ;
                  vector<int> LC = VideoNode::getIntFromId(nodeId);
                  vector<INT> cp, pp; 
                  for(int i=LC[0]; i < LC[0]+LC[2]; i++) {
                      for(int j=LC[1]; j < LC[1]+LC[3]; j++) {
                          bool res = getAbsolutePixel(parId, i, j, parState, &pp);
                          bool res2 = getAbsolutePixel(nodeId, i, j, cStPtr, &cp); 
                          cout<<" resp: "<<res<<" resc: "<<res2<<" pp[0]: "<<((int)pp[0])<<" cp[0]: "<<((int)cp[0])<<"\n";
                         
                          // getchar(); 
                      }
                  }
                  cout<<"parent next state map:\n";
                  ConditionalBeliefPtr cptr = vn->getCPTR();
                  ASSERT(cptr != 0); 
                  BeliefPtr cblf = cptr->getBeliefForState(oldState);
                  ASSERT(cblf != 0); 
                  map<StatePtr, LD, StatePtrCmp> mp = cblf->getMp();
                  foreach(localIter, mp) {
                      cout<<"\tst: "<<*(localIter->first)<<" p: "<<localIter->second<<endl; 
                  }*/
                pChild = stMpIter->second;
                // augmented child states
                vn->getCPTR()->setPr(oldState, cStPtr, pChild);
                parAugmentNotFound++;
                // cout<<" case 2 pChild: "<<pChild<<"\n";
                // getchar();

            } else {
                // cout<<" no problemo\n";
                pChild = pChild/pFac;

            }

            bool res = vn->addAugmentCand(cStPtr, pChild);

            if(res == true) {
                parAugmentFound++;
                /*
                cout<<"par: "<<parId<<" st: "<<*(stMpIter->first)<<"\n"; 
                cout<<"chd: "<<nodeId<<" st: "<<*cStPtr<<" p: "<<pChild<<" pFac: "<<pFac<<"\n";
                //getchar();*/
            }
        }
    }
}


/// holes(grayscale) convention = 255 (hole) ... 0 (no hole)
void VideoGraph::getNextFramePPM(IplImage* frame, IplImage* mask_frame) {
    static int calls=0;
    calls++;
    int count = 0;
    LD candT = 0.10;
    //   FILE* file = fopen("GNFP.log", "a+");
STATE_INIT:
    PPM_step2(frame, candT);
MSG_REINIT:
    PPM_step3();
INNER_LOOP:
    PPM_inner_loop();
    //    PPM_LN(1.0e20);
    //    updateNodesInfo();
    //fillFrame1(false);
    fillFrame0();
    count++;
    string report = get_status();
    string str = get_status();
    //  double diff = reader->getFrameDifference(true_frame, output_frame);
    //  cout<<"PPM_iteration: "<<count<<" diff: "<<diff<<"\n";

    cout<<"stats: ["<<str<<"]\n";
    // FPRINTF(file, "%d %d %f %s\n", calls, count, diff, str.c_str());
    // fflush(file);
    //  if( (avgChngInBlfs > 0.01) && (count <= 2) ) goto INNER_LOOP;
    if( (count <= 2) )
        goto INNER_LOOP;
    //fclose(file);
    //cout<<"another iteration(y/n): "; // if(getchar() != 'n') //    goto INNER_LOOP;
}

void VideoGraph::fillFrame1(bool lookMax) {
    Timer timer;
    FPRINTF(stderr, "fillFrame1() fills by looking at all big regions\n");
    int div = 256/numLvls;
    if(isGray) {
        BwImage image(output_frame);
        LD pixels[N1][N2][numLvls];
        for(int i=0; i < N1; i++)
            for(int j=0; j < N2; j++)
                for(int k=0; k < numLvls; k++)
                    pixels[i][j][k] = 0.0;
        foreach(iter, adj_map) {
            VideoNode* vn = (VideoNode*) iter->first;
            if(vn->getNodeNumVar() != (rx*ry))
                continue;
            vector<int> L = VideoNode::getIntFromId(vn->getNodeId());
            map<StatePtr, LD, StatePtrCmp> cMp = vn->bNew->getMp();
            LD localPixels[L[2] ][L[3] ][numLvls];
            fi(0, L[2]) fj(0, L[3]) fk(0, numLvls) localPixels[i][j][k] = 0;
            foreach(mpIter, cMp) {
                StatePtr stPtr = mpIter->first;
                LD currPr = mpIter->second;
                fi(0, L[2]) {
                    fj(0, L[3]) {
                        int idx = i*L[3] + j;
                        INT pVal = stPtr->getData(idx);
                        int array_idx = (( (int) pVal)-0)/div;
                        if(!lookMax) {
                            pixels[i+L[0] ][j+L[1] ][array_idx] += currPr;
                        } else {
                            localPixels[i][j][array_idx] += currPr;
                        }
                    }
                }
            }
            if(lookMax) {
                fi(0, L[2]) {
                    fj(0, L[3]) {
                        fk(0, numLvls) {
                            LD opr = pixels[i+L[0] ][j+L[1] ][k];
                            LD npr = localPixels[i][j][k];
                            pixels[i+L[0] ][j+L[1] ][k] = max(opr, npr);
                        }
                    }
                }
            }
        }
        fi(0, N1) {
            fj(0, N2) {
                int maxIdx = 0;
                LD maxVal = 0.0;
                fk(0, numLvls) {
                    if(pixels[i][j][k] > maxVal) {
                        maxVal = pixels[i][j][k];
                        maxIdx = k;
                    }
                }
                INT pVal = (INT) (maxIdx*div);
                image[i][j] = pVal;
            }
        }
        /*
        cvNamedWindow("fillFrame1", CV_WINDOW_AUTOSIZE);
        cvShowImage("fillFrame1", image.getImg());
        cvWaitKey(1000);
        cvDestroyWindow("fillFrame1");
        */
    } else {
        RgbImage image(output_frame);
        LD pixels[N1][N2][numLvls][3];
        fi(0, N1) fj(0, N2) fk(0, numLvls) for(int c=0; c < 3; c++)
            pixels[i][j][k][c] = 0.0;
        foreach(iter, adj_map) {
            VideoNode* vn = (VideoNode*) iter->first;
            if(vn->getNodeNumVar() != (rx*ry))
                continue;
            vector<int> L = VideoNode::getIntFromId(vn->getNodeId());
            map<StatePtr, LD, StatePtrCmp> cMp = vn->bNew->getMp();
            LD localPixels[L[2] ][L[3] ][numLvls][3];
            fi(0, L[2]) fj(0, L[3]) fk(0, numLvls) for(int c=0; c< 3; c++)
                localPixels[i][j][k][c] = 0;
            foreach(mpIter, cMp) {
                StatePtr stPtr = mpIter->first;
                LD currPr = mpIter->second;
                fi(0, L[2]) {
                    fj(0, L[3]) {
                        for(int c=0; c< 3; c++) {
                            int idx = i*L[3] + j;
                            INT pVal = stPtr->getData(3*idx+c);
                            int array_idx = (( (int) pVal)-0)/div;
                            if(!lookMax) {
                                pixels[i+L[0] ][j+L[1] ][array_idx][c] += currPr;
                            } else {
                                localPixels[i][j][array_idx][c] += currPr;
                            }
                        }
                    }
                }
            }
            if(lookMax) {
                fi(0, L[2]) {
                    fj(0, L[3]) {
                        fk(0, numLvls) {
                            for(int c=0; c< 3; c++) {
                                LD opr = pixels[i+L[0] ][j+L[1] ][k][c];
                                LD npr = localPixels[i][j][k][c];
                                pixels[i+L[0] ][j+L[1] ][k][c] = max(opr, npr);
                            }
                        }
                    }
                }
            }
        }
        fi(0, N1) {
            fj(0, N2) {
                int maxIdx[3];
                LD maxVal[3];
                for(int c=0; c < 3; c++) {
                    maxIdx[c] = 0;
                    maxVal[c] = 0.0;
                }
                fk(0, numLvls) {

                    for(int c=0; c < 3; c++) {

                        if(pixels[i][j][k][c] > maxVal[c]) {
                            maxVal[c] = pixels[i][j][k][c];
                            maxIdx[c] = k;
                        }
                    }
                }

                image[i][j].r = ((INT) (maxIdx[0]*div));
                image[i][j].g = ((INT) (maxIdx[1]*div));
                image[i][j].b = ((INT) (maxIdx[2]*div));
            }
        }
        /*
        cvNamedWindow("fillFrame1", CV_WINDOW_AUTOSIZE);
        cvShowImage("fillFrame1", image.getImg());
        cvWaitKey(1000);
        cvDestroyWindow("fillFrame1");
        */
    }


}

void VideoGraph::fillFrame0() {
    FPRINTF(stderr, "fillFrame0()\n");
    if(isGray) {
        //    IplImage* local_frame=cvCreateImage(cvSize(N2,N1), IPL_DEPTH_8U,1);
        BwImage img(output_frame);
        bool filled[N1][N2];
        fi(0, N1) fj(0, N2) filled[i][j]=false;
        foreach(iter, adj_map) {
            VideoNode* vn = (VideoNode*) iter->first;
            if(vn->getNodeNumVar() != (rx*ry))
                continue;
            else {
                vector<int> L = VideoNode::getIntFromId(vn->getNodeId());
                pair<StatePtr, LD> bestPair = getBestState0(vn);

                StatePtr bestState = bestPair.first;
                LD bestPr = bestPair.second;

                for(int i=0; i < rx; i++) {
                    for(int j=0; j < ry; j++) {
                        int x = i+ L[0];
                        int y = j+ L[1];

                        if(filled[x][y])
                            continue;
                        else {
                            int idx = i*L[3]+j;
                            img[x][y] = bestState->getData(idx);
                            filled[x][y] = true;
                            //   printf("x: %d y: %d pixel: %d prob: %f\n", x, y, ((int) img[x][y]), bestPr);
                        }
                    }
                }
            }
        }
        //  cvNamedWindow("fillFrame0", CV_WINDOW_AUTOSIZE);
        //  cvShowImage("fillFrame0", img.getImg());
        //  cvWaitKey(0);
        //  cvDestroyWindow("fillFrame0");
        //   cvCopy(local_frame, output_frame, 0);
        //   cvReleaseImage(&local_frame);

    } else {
        RgbImage img(output_frame);
        bool filled[N1][N2];
        fi(0, N1) fj(0, N2) filled[i][j]=false;
        foreach(iter, adj_map) {
            VideoNode* vn = (VideoNode*) iter->first;
            if(vn->getNodeNumVar() != (rx*ry))
                continue;
            else {
                vector<int> L = VideoNode::getIntFromId(vn->getNodeId());
                pair<StatePtr, LD> bestPair = getBestState0(vn);
                StatePtr bestState = bestPair.first;
                ASSERT(bestState->num == (3*rx*ry));
                LD bestPr = bestPair.second;
                for(int i=0; i < rx; i++) {
                    for(int j=0; j < ry; j++) {
                        int x = i+ L[0];
                        int y = j+ L[1];

                        if(filled[x][y])
                            continue;
                        else {
                            int idx = i*L[3]+j;
                            img[x][y].r = bestState->getData(3*idx);
                            img[x][y].g = bestState->getData(3*idx+1);
                            img[x][y].b = bestState->getData(3*idx+2);
                            filled[x][y] = true;
                            //              printf("x: %d y: %d pixel.r: %d prob: %f\n", x, y, ((int) img[x][y].r), bestPr);
                            //              printf("x: %d y: %d pixel.g: %d prob: %f\n", x, y, ((int) img[x][y].g), bestPr);
                            //              printf("x: %d y: %d pixel.b: %d prob: %f\n", x, y, ((int) img[x][y].b), bestPr);
                            //              getchar();
                        }
                    }
                }
            }
        }
        //   cvNamedWindow("fillFrame0", CV_WINDOW_AUTOSIZE);
        //    cvShowImage("fillFrame0", img.getImg());

    }

}

pair<StatePtr, LD> VideoGraph::getBestState0(VideoNode* vn) {
    map<StatePtr, LD, StatePtrCmp> bMp = vn->bNew->getMp();
    StatePtr bestRes = 0;
    LD bestP = 0.0;
    foreach(iter, bMp) {
        if(iter->second > bestP) {
            bestP = iter->second;
            bestRes = iter->first;
        }
    }
    pair<StatePtr, LD> result = make_pair<StatePtr, LD>(bestRes, bestP);
    return result;
}

bool VideoGraph::isMaskedPixel(int x, int y, IplImage* mask_frame) {
    BwImage image(mask_frame);
    return (image[x][y] != 0);
}

vector<pair<int, int> > VideoGraph::getNodeUnfilledPixels(VideoNode* vn, IplImage* mask_frame) {
    vector<int> L = VideoNode::getIntFromId(vn->getNodeId());
    vector<pair<int, int> > result;
    for(int i=L[0]; i < (L[0]+L[2]); i++) {
        for(int j=L[1]; j < (L[1] + L[3]); j++) {
            if(isMaskedPixel(i, j, mask_frame) ) {
                pair<int, int> cp = make_pair<int, int>(i, j);
                result.push_back(cp);
            }
        }
    }
    return result;

}

vector<pair<StatePtr, LD> > VideoGraph::getMatchingStates(VideoNode* vn, IplImage* frame, IplImage* mask_frame, bool approx) {
    Timer timer;
    cerr<<"getMatchingStates: "<<vn->getNodeId()<<" ";
    ASSERT(pOrigMp.find(vn->getNodeNumVar()) != pOrigMp.end());
    BeliefPtr bptr = pOrigMp.find(vn->getNodeNumVar())->second;
    map<StatePtr, LD, StatePtrCmp> mp = bptr->getMp();
    vector<pair<StatePtr, LD> > result;
    vector<int> L = VideoNode::getIntFromId(vn->getNodeId());
    BwImage maskImg(mask_frame);
    string nid = vn->getNodeId();
    int passed =0;
    int failed = 0;
    cout<<" mp.size(): "<<mp.size()<< " ";
    foreach(mpIter, mp) {
        StatePtr cptr = mpIter->first;
        StatePtr fptr = readNodeFromFrame(vn , frame);
        int i = L[0], j = L[1];
        bool misMatch = false;
        while((!misMatch) && (i < (L[0]+L[2])) ) {
            while((!misMatch) && (j < (L[1]+L[3])) ) {
                vector<INT> r1, r2;
                timer.tick();
                if(maskImg[i][j] == 0) { // check only when masking true
                    getAbsolutePixel(nid, i, j, cptr, &r1);
                    getAbsolutePixel(nid, i, j, fptr, &r2);

                    //cout<<"r1: "; foreach(r1it, r1) cout<<*r1it<< " ";
                    //cout<<"\nr2: "; foreach(r2it, r2) cout<<*r2it<< " ";
                    //getchar();
                    ASSERT(r1.size() == r2.size());
                    fk(0, r1.size()) {
                        if(r1[k] != r2[k])
                            misMatch = true;
                    }
                }
                j++;
            }
            i++;
        }
        if(misMatch == false) {
            failed++;
            result.push_back(make_pair<StatePtr, LD>(cptr, mpIter->second));
        } else {
            passed++;
        }
    }
    cout<<" failed: "<<failed<<" passed: "<<passed<<" res.size(): "<<result.size();
    //getchar();
    if((approx == false) || (failed > 0))  {
        sort(result.begin(), result.end(), PRCMP2<StatePtr, LD>(false));
        cerr<<"\n";
        return result;
    }

    LD maxDiff = 10000000.0;
    StatePtr best=NULL;
    foreach(mpIter, mp) {
        StatePtr cptr = mpIter->first;
        StatePtr fptr = readNodeFromFrame(vn , frame);
        LD currDiff = 0.0;
        for(int i=L[0]; i < (L[0]+L[2]); i++) {
            for(int j=L[1]; j < (L[1]+L[3]); j++) {
                vector<INT> r1, r2;
                timer.tick();
                if(maskImg[i][j] == 0) { // check only when masking true
                    getAbsolutePixel(nid, i, j, cptr, &r1);
                    getAbsolutePixel(nid, i, j, fptr, &r2);
                    fk(0, r1.size()) {
                        currDiff += abs(((int) r1[k])-((int) r2[k]));
                    }
                }
            }
        }
        if(currDiff < maxDiff) {
            maxDiff = currDiff;
            best = cptr;
        }
    }

    LD p = mp.find(best)->second;
    result.push_back(make_pair<StatePtr, LD>(best, p));
    return result;
}

vector<pair<StatePtr, LD> > VideoGraph::getMaximalPixelPr(VideoNode* vn, int i, int j) {
    // cout<<" node: "<<vn->getNodeId()<<" i: "<<i<<" j: "<<j<<endl; getchar();
    vector<pair<StatePtr, LD> > result;
    map<int, Belief*>::iterator miter = pOrigMp.find(vn->getNodeNumVar());
    ASSERT(miter != pOrigMp.end());
    BeliefPtr bptr = miter->second;
    //vn->getGlobalBlfPtr();
    ASSERT(bptr != 0);
    map<StatePtr, LD, StatePtrCmp> prMp;
    map<StatePtr, LD, StatePtrCmp> bMp = bptr->getMp();
    vector<int> L = VideoNode::getIntFromId(vn->getNodeId());
    int idx = (i-L[0])*L[3] + (j-L[1]);
    ASSERT(bMp.size() > 0);
    foreach(iter, bMp) {
        StatePtr sptr = iter->first;
        StatePtr resPtr = NULL;
        if(isGray) {
            ASSERT(sptr->num == vn->getNodeNumVar());
            resPtr = getStatePtr(sptr->getData(idx));
        } else {
            ASSERT(sptr->num == 3*vn->getNodeNumVar());
            vector<INT> pvals;
            pvals.push_back(sptr->getData(idx*3));
            pvals.push_back(sptr->getData(idx*3+1));
            pvals.push_back(sptr->getData(idx*3+2));
            resPtr = getStatePtr(pvals);
        }
        if(prMp.find(resPtr) == prMp.end())
            prMp.insert(make_pair<StatePtr, LD>(resPtr, 0.0));
        map<StatePtr, LD,StatePtrCmp>::iterator miter= prMp.find(resPtr);
        miter->second = miter->second + 1.0;
    }
    foreach(iter, prMp) {
        result.push_back(make_pair<StatePtr, LD>(iter->first, iter->second));
    }
    ASSERT(result.size() > 0);
    sort(result.begin(), result.end(), PRCMP2<StatePtr, LD>(false));
    return result;
}

void VideoGraph::clearAllNodes() {
    foreach(iter, pOrigMp) {
        iter->second->clear();
    }
}

void VideoGraph::clearCondBlfMp() {
    foreach(iter, blfMp) {
        iter->second->clear();
    }
}

IplImage* VideoGraph::inPaint(IplImage* input_frame, IplImage* mask_frame) {
    Timer timer;
    cout<<"inPaint(): \n";
    vector<pair<VideoNode*, LD>  > listP;
    cvNamedWindow("inPaint_pred", CV_WINDOW_AUTOSIZE);
    cvNamedWindow("maskFrame", CV_WINDOW_AUTOSIZE);
    if(local_mask == NULL)
        local_mask=cvCreateImage(cvSize(N2,N1), IPL_DEPTH_8U,1);
    if(inpainted_frame == NULL)
        inpainted_frame = cvCreateImage(cvSize(N2, N1), input_frame->depth, input_frame->nChannels);
    updateMaskStatus(mask_frame);
    cvCopy(input_frame, inpainted_frame, 0);
    cvCopy(mask_frame, local_mask, 0);
    BwImage local_image(local_mask);
    listP = getIncompleteList(local_mask);
    int filledNum = 0;
    set<string> old_failed, new_failed;
    bool approx = false;
LOOP:
    /*
        listP.clear(); 
        listP = getIncompleteList(local_mask); 
        if(listP.empty()) { ; }
        else {*/

    VideoNode* vn = listP[0].first;
    vector<pair<int, int> > unfilledPixels = getNodeUnfilledPixels(vn, local_mask);
    if(unfilledPixels.size() > 0) {
        vector<pair<StatePtr, LD> > validStates = getMatchingStates(vn, inpainted_frame, local_mask, approx);
        if(validStates.size() > 0) {
            vector<int> L = VideoNode::getIntFromId(vn->getNodeId());
            cout<<"vn: "<<vn->getNodeId()<<" validStates.size(): "<<validStates.size()<<"\n";

            StatePtr optState = validStates[0].first;
            vector<INT> pval;
            foreach(iter, unfilledPixels) {

                int x = iter->first;
                int y = iter->second;
                cout<<" updated pixel: ("<<x<<", "<<y<<")\n";
                getAbsolutePixel(vn->getNodeId(), x, y, optState, &pval);
                if(isGray) {
                    BwImage image(inpainted_frame);
                    image[x][y] = pval[0];
                } else {
                    RgbImage image(inpainted_frame);


                    image[x][y].r = pval[0];
                    image[x][y].g = pval[1];
                    image[x][y].b = pval[2];
                }
                local_image[x][y] = 0;
            }
            filledNum +=  unfilledPixels.size();
            cout<<" total pixels filled: "<<     filledNum<<endl;
        } else {
            new_failed.insert(vn->getNodeId());

        }

    }
    cvSaveImage("temp.png", inpainted_frame);
    listP.erase(listP.begin());
    if(listP.empty()) {
        listP = getIncompleteList(local_mask);
        if(new_failed == old_failed)
            approx = true;
        else {
            old_failed.clear();
            old_failed = new_failed;
            new_failed.clear();
        }
    }
    if(listP.size() > 0) {
        goto LOOP;
    }

    /* // single fill
    vector<pair<StatePtr, LD> > pixelPr = getMaximalPixelPr(vn, cpixel.first, cpixel.second);
    // set the pixel value
    int x = cpixel.first; 
    int y = cpixel.second;
    if(isGray) {
        BwImage image(inpainted_frame); 
        image[cpixel.first][cpixel.second] = pixelPr[0].first->getData(0); 
} else {
        RgbImage image(inpainted_frame); 
        StatePtr sptr = pixelPr[0].first; 
        ASSERT(sptr->num == 3); 
        image[cpixel.first][cpixel.second].r = pixelPr[0].first->getData(0); 
        image[cpixel.first][cpixel.second].g = pixelPr[0].first->getData(1); 
        image[cpixel.first][cpixel.second].b = pixelPr[0].first->getData(2); 
}
    cout<<"filled pixel: "<<cpixel.first<<", "<<cpixel.second<<"\n";
    local_image[x][y] = 0; 
    cvShowImage("inPaint_pred", inpainted_frame); 
    cvShowImage("maskFrame", local_mask); 
    // cvWaitKey(1000);
    // getchar();  
    goto LOOP;
}
    */
ENDLOOP:
    return inpainted_frame;
}

bool VideoGraph::isPartialMatch(StatePtr s1, StatePtr s2, StatePtr mask) {
    int N = mask->num;
    ASSERT( (s1->num == s2->num));
    if(mask->num == s1->num) {
        fi(0, N) {
            if( (mask->getData(i) == 0) && (s1->getData(i) != s2->getData(i)) )
                return false;
        }
    } else {
        fi(0, N) {
            if( (mask->getData(i) == 0) &&
                    ( (s1->getData(3*i) != s2->getData(3*i))
                      || (s1->getData(3*i+1) != s2->getData(3*i+1))
                      || (s1->getData(3*i+2) != s2->getData(3*i+2))  ))
                return false;
        }


    }
    return true;
}

void VideoGraph::getCandidates(IplImage* last_frame, IplImage* input_frame,  IplImage* mask_frame) {
    updateMaskStatus(mask_frame);
    foreach(iter, adj_map) {
        VideoNode* vn= (VideoNode*) iter->first;
        vn->getBeliefPtr()->clear();
        vn->pHat->clear();
        vn->bNew->clear();
        int numCands=0;
        string mode="";
        StatePtr lastState = readNodeFromFrame(vn, last_frame);
        StatePtr nextState = readNodeFromFrame(vn, input_frame);
        int N = lastState->num;
        vn->getBeliefPtr()->setPr(lastState, 1.0);
        if(! (vn->isMasked) ) {
            numCands++;
            vn->pHat->setPr(lastState, nextState, 1.0);
            vn->bNew->setPr(nextState, 1.0);
            mode = "UNMASKED";
        } else {


            StatePtr maskState = readNodeFromFrame(vn, mask_frame);
            LD sumP = 0.0;
            vector<INT> L = VideoNode::getIntFromId(vn->getNodeId());
            
            if(numCands == 0) {
                ConditionalBelief* cbptr = blfMp.find(N)->second;
                map<StatePtr, Belief, StatePtrCmp> *cbMp =cbptr->getCbMpPtr();
                if(cbMp->find(lastState) != cbMp->end()) {
                    map<StatePtr, LD, StatePtrCmp> stMp = cbMp->find(lastState)->second.getMp();
                    foreach(stIter, stMp) {
                        if((vn->isMasked == false) || (isPartialMatch(stIter->first, nextState, maskState) )) {
                            vn->bNew->setPr(stIter->first, 0.1);
                            vn->pHat->setPr(lastState, stIter->first, stIter->second);
                            numCands++;
                            sumP +=  stIter->second;
                        }
                        mode = "CBMP_GLOBAL";
                    }
                    map<StatePtr, LD, StatePtrCmp> nMp = vn->bNew->getMp();
                    int K = nMp.size();
                    LD pEq = 1.0/((LD) K);
                    foreach(stIter, nMp) {
                        vn->bNew->setPr(stIter->first, pEq);
                        vn->pHat->setPr(lastState, stIter->first, stIter->second/sumP);
                    }
                }
            }
            if(numCands == 0) {
                ConditionalBelief* currCbPtr =  (vn->localCondBlfPtr);
                map<StatePtr, Belief, StatePtrCmp> *cbMp =currCbPtr->getCbMpPtr();
                if(cbMp->find(lastState) != cbMp->end()) {
                    map<StatePtr, LD, StatePtrCmp> stMp = cbMp->find(lastState)->second.getMp();

                    foreach(stIter, stMp) {
                        if((vn->isMasked == false) || (isPartialMatch(stIter->first, nextState, maskState) )) {
                            vn->bNew->setPr(stIter->first, 0.1);
                            vn->pHat->setPr(lastState, stIter->first, stIter->second);
                            numCands++;
                            sumP +=  stIter->second;
                        }
                        mode = "CBMP_LOCAL";
                    }

                    map<StatePtr, LD, StatePtrCmp> nMp = vn->bNew->getMp();
                    int K = nMp.size();
                    LD pEq = 1.0/((LD) K);
                    foreach(stIter, nMp) {
                        vn->bNew->setPr(stIter->first, pEq);
                        vn->pHat->setPr(lastState, stIter->first, stIter->second/sumP);
                    }
                }
            }

            if(numCands == 0) {
                map<StatePtr, LD, StatePtrCmp> bMp = vn->localBlfPtr->getMp();
                vector<pair<StatePtr, LD> > exact_matches;
                LD sumPD = 0.0;
                foreach(stIter, bMp) {
                    LD diff= 0.0;
                    fi(0, stIter->first->num) {
                        if(maskState->getData(i) == 0) {
                            diff += abs( (stIter->first->getData(i) - nextState->getData(i))/((LD) numLvls));
                        }
                    }
                    LD pCurr = (stIter->second)/(1.0+diff*diff);
                    sumPD += pCurr;
                    exact_matches.push_back(make_pair<StatePtr, LD>(stIter->first, pCurr));
                }
                LD pEq = 1.0/((LD) exact_matches.size());
                foreach(stIter, exact_matches) {
                    vn->pHat->setPr(lastState, stIter->first, stIter->second/sumPD);
                    vn->bNew->setPr(stIter->first, pEq);
                    numCands++;
                    mode="APPROX";
                }
            }
         
        }
        if(mode != "UNMASKED")
            cout<<"vn: "<<vn->getNodeId()<<" mode: "<<mode<<" numCands: "<<numCands<<"\n";
    }
}


void VideoGraph::initLocals() {
    foreach(iter, adj_map) {
        VideoNode* vn = (VideoNode*) iter->first;
        vn->localBlfPtr = new Belief(true);
        vn->localCondBlfPtr = new ConditionalBelief();
    }
}

/*   vector<pair<StatePtr, LD> > matches = getMatchingStates(vn, input_frame, mask_frame, false);
                   LD sumP = 0.0; 
                   if(matches.size() > 0.0) {
                       LD pEq = 1.0/((LD) matches.size()); 
                       foreach(stIter, matches) { sumP += stIter->second; }
                       for(int i=0; i < min((int)matches.size(), 10); i++) {
                           //foreach(stIter, matches) {
                           StatePtr currSt = matches[i].first; 
                           LD pSt = matches[i].second;
                     //      cout<<" i: "<<i<<" st: "<<*currSt<<" p: "<<pSt<<endl; getchar();  
                           vn->bNew->setPr(currSt, pEq);  
                           vn->pHat->setPr(lastState, currSt, pSt);  
                           numCands++;
                          mode="EXACT";  
                       }
                   } else { // what to do 
                        // cout<<" approximate match required: "<<vn->getNodeId()<<" state: "<<*lastState<<endl; getchar();  
                        vector<pair<StatePtr, LD> > matches = getMatchingStates(vn, input_frame, mask_frame, true);
                        vn->pHat->setPr(lastState, matches[0].first, matches[0].second);  
                        vn->bNew->setPr(matches[0].first, 0.5);
                        numCands++;  
                        mode="APPROX"; 
                   }*/
