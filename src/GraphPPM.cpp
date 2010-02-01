#include "GraphPPM.h"
#include "main.h"

/**
 * Initializes child->parent messages (child->BidirEdgeList.second)
 * to 1.0 for all nodes. 
 * Note: PPM/GBP child->parent algorithm does not use forward messages 
 * from parent->child. 
 */
void GraphPPM::PPM_step3() {
    Timer timer;
    FPRINTF(stderr, "PPM_step3() initializes messages\n");
    foreach(iter, adj_map) {
        RegionPPM* vn = (RegionPPM*) iter->first;
        vector<MyEdge> inEdges = iter->second.second;
        foreach(edgeIter, inEdges) {
            edgeIter->msgMp.clear();
        }
        vector<MyEdge> outEdges = iter->second.first;
        //foreach(edgeIter, outEdges) { edgeIter->msgMp.clear(); }
        // write the parent lambda
        //cout<<"\nnode: "<<vn->getNodeId()<<"\n";
        foreach(edgeIter, iter->second.second) {
            /*map<StatePtr, LD, StatePtrCmp> oldStateMp = vn->getBeliefPtr()->getMp();
            foreach(stIter, oldStateMp) {
                BeliefPtr nStBlf = vn->pHat->getBeliefForState()
        } */
            //  cout<<"edge: "<<*edgeIter<<"\n";
            map<StatePtr, Belief, StatePtrCmp> *cbMpPtr =   vn->pHat->getCbMpPtr();
            foreach(stIter, *cbMpPtr) {
                StatePtr currSt = stIter->first;
                BeliefPtr nStBlf = vn->pHat->getBeliefForState(currSt);
                map<StatePtr, LD, StatePtrCmp> newStateMp = nStBlf->getMp();
                foreach(newStIter, newStateMp) {
                    StatePtr nextSt = newStIter->first;
                    StatePtr jointSt = getJointState(currSt, nextSt);
                    //  iter->second.second.msgMp.insert(make_pair<StatePtr, LD>(jointSt, 1.0));
                    edgeIter->setMsg(jointSt, 1.0);
                    timer.tick();

                }
            }
        }

        //       getchar();
    }

}
// compute the parent lambda
void GraphPPM::PPM_step4a() {
    Timer timer;
    FPRINTF(stderr, "PPM_step4a() updates lambdaPar\n");
    int found=0, notFound=0;
    foreach(iter, adj_map) {
        RegionPPM* vn = (RegionPPM*) iter->first;
        vector<MyEdge> inEdges = iter->second.second;
        map<StatePtr, Belief, StatePtrCmp> *cbMpPtr =   vn->pHat->getCbMpPtr();
        foreach(stIter, *cbMpPtr) {
            StatePtr currSt = stIter->first;
            BeliefPtr nStBlf = vn->pHat->getBeliefForState(currSt);
            map<StatePtr, LD, StatePtrCmp> newStateMp = nStBlf->getMp();
            foreach(newStIter, newStateMp) {
                StatePtr nextSt = newStIter->first;
                StatePtr jointSt = getJointState(currSt, nextSt);
                LD value = 1.0;
                foreach(edgeIter, inEdges) {
                    timer.tick();
                    LD currVal = edgeIter->getMsg(jointSt);
                    if(currVal > 0.0 ) {
                        value *= currVal;
                        // cout<<" PPM_step4a() node: "<<vn->getNodeId()<<" currVal: "<<currVal<<endl; getchar();
                        found++;
                    } else {
                        notFound++;
                    }

                }

                vn->lambdaPar->setPr(currSt, nextSt, value);

            }
        }
    }

    FPRINTF(stderr, "PPM_step4a() found: %d notFound: %d\n", found, notFound);
}



// update child contribution
void GraphPPM::PPM_step4b() {
    Timer timer;
    FPRINTF(stderr, "PPM_step4b() updates lambdaChild\n");
    int found=0, notFound = 0;
    foreach(iter, adj_map) {
        RegionPPM* vn = (RegionPPM*) iter->first;
        string nodeId = vn->getNodeId();
        vector<MyEdge> outEdges = iter->second.first;
        map<StatePtr, Belief, StatePtrCmp> *cbMpPtr =   vn->pHat->getCbMpPtr();

        foreach(stIter, *cbMpPtr) {
            StatePtr currStPar = stIter->first;
            BeliefPtr nStBlf = vn->pHat->getBeliefForState(currStPar);
            map<StatePtr, LD, StatePtrCmp> newStateMp = nStBlf->getMp();
            foreach(newStIter, newStateMp) {
                StatePtr nextStPar = newStIter->first;
                LD value = 1.0;
                foreach(edgeIter, outEdges) {
                    NodePtr vc = edgeIter->ends[1];
                    pair<EdgePtr, EdgePtr> edgePr = findEdge(vn, vc);
                    EdgePtr backEdge = edgePr.second;
                    string childId = edgeIter->ends[1]->getNodeId();
                    StatePtr currSt = getChildStatePtr(vn, vc, currStPar);
                    StatePtr nextSt = getChildStatePtr(vn, vc, nextStPar);
                    StatePtr jointSt = getJointState(currSt, nextSt);
                    LD currVal = backEdge->getMsg(jointSt);
                    timer.tick();
                    if(currVal > 0.0) {
                        value *= currVal;
                        //cout<<" PPM_step4b() vn: "<<nodeId<<" vc: "<<childId<<" sJoint: "<<*jointSt<<endl;
                        found++;
                    } else {
                        notFound++;
                    }
                }

                vn->lambdaChild->setPr(currStPar, nextStPar, value);
            }
        }
    }

    FPRINTF(stderr, "PPM_step4b() found: %d notFound: %d\n", found, notFound);
}


void GraphPPM::PPM_step5a() {
    Timer timer;
    FPRINTF(stderr, "PPM_step5a() updates lambdaR\n");
    foreach(iter, adj_map) {
        RegionPPM* vn = (RegionPPM*) iter->first;
        string nodeId = vn->getNodeId();
        int CR = vn->getCR();
        map<StatePtr, Belief, StatePtrCmp> *cbMpPtr =   vn->pHat->getCbMpPtr();
        //   cout<<"node: "<<nodeId<<" cbMpPtr.size(): "<<cbMpPtr->size()<<endl;
        // GBP MODIFICATION BEGIN 
        LD sumC = 0.0; 
        // GBP MODIFICATION END 
        foreach(stIter, *cbMpPtr) {
            StatePtr currSt = stIter->first;
            BeliefPtr nStBlf = vn->pHat->getBeliefForState(currSt);
            map<StatePtr, LD, StatePtrCmp> newStateMp = nStBlf->getMp();
            
            // ORIGINAL PPM CODE BEGIN  
            // LD sumC = 0.0;
            // ORIGINAL PPM CODE END
            
            foreach(newStIter, newStateMp) {
                StatePtr nextSt = newStIter->first;
                LD pHatC = vn->pHat->getPr(currSt, nextSt);
                ASSERT(pHatC >= 0.0);
                LD parC = vn->lambdaPar->getPr(currSt, nextSt);
                ASSERT(parC >= 0.0);
                LD childC = vn->lambdaChild->getPr(currSt, nextSt);
              	// code changed here 
                ASSERT(childC >= 0.0);

                LD currC = (LD) pHatC*pow(((double)(parC*childC)), ((1.0)/((double) CR)));
                sumC += currC;
                timer.tick();
            }

            LD bCurr = vn->getBeliefPtr()->getPr(currSt);
            //      cout<<"node: "<<nodeId<<" state: "<<*currSt<<" bCurr: "<<bCurr<<endl;
            ASSERT(bCurr != BELIEF_NOT_FOUND);
            
            // GBP MODIFICATION 
            // LD newVal = 1.0/sumC;
            // ORIGINAL PPM CODE
            LD newVal = bCurr/sumC;
            // 
            vn->lambdaR->setPr(currSt, newVal);
        }

    }

}

void GraphPPM::PPM_step5b() {
    Timer timer;
    FPRINTF(stderr, "PPM_step5b() updates the joint belief\n");
    foreach(iter, adj_map) {

        RegionPPM* vn = (RegionPPM*) iter->first;
        vn->bJointRev->clear();
        string nodeId = vn->getNodeId();
        int CR = vn->getCR();
        ConditionalBelief cbTmp;
        map<StatePtr, Belief, StatePtrCmp> *cbMpPtr =   vn->pHat->getCbMpPtr();
        LD sumC = 0.0;
        
        
        foreach(stIter, *cbMpPtr) {
            StatePtr currSt = stIter->first;
            BeliefPtr nStBlf = vn->pHat->getBeliefForState(currSt);
            map<StatePtr, LD, StatePtrCmp> newStateMp = nStBlf->getMp();
            LD lambda_rC =  vn->lambdaR->getPr(currSt);
            ASSERT(lambda_rC != BELIEF_NOT_FOUND);
            foreach(newStIter, newStateMp) {
                StatePtr nextSt = newStIter->first;
                if(lambda_rC > 0.0) {
                       
                    LD pHatC = vn->pHat->getPr(currSt, nextSt);
                    ASSERT(pHatC >= 0.0);
                    LD parC = vn->lambdaPar->getPr(currSt, nextSt);
                    ASSERT(parC >= 0.0);
                    LD childC = vn->lambdaChild->getPr(currSt, nextSt);
                    ASSERT(childC >= 0.0);
                    LD ngbrC = (LD) pow(((double)(parC*childC)), ((1.0)/((double) CR)));
                    LD currC = pHatC*lambda_rC*ngbrC;
                    cbTmp.setPr(currSt, nextSt, currC);
                    sumC += currC;
                    timer.tick();
                } else {
                   cbTmp.setPr(currSt, nextSt, 0.0);
                }
           //     cout<<nodeId<<"\n\tcSt: "<<currSt<<"\n\tnSt:"<<nextSt<<"\npHat: "<<pHatC<<" parC: "<<parC<<" chdC: "<<childC<<" lR: "<<lambda_rC<<"\n"; 
           //     getchar();
            }
        }
        if(sumC == 0.0) {
            foreach(stIter, *cbMpPtr) {
                StatePtr currSt = stIter->first;
                BeliefPtr nStBlf = vn->pHat->getBeliefForState(currSt);
                map<StatePtr, LD, StatePtrCmp> newStateMp = nStBlf->getMp();
                cout<<" node: "<<vn->getNodeId(); 
                foreach(newStIter, newStateMp) {
                    StatePtr nextSt = newStIter->first;
                    LD pHatC = vn->pHat->getPr(currSt, nextSt);
                    cout<<" cST: "<<*currSt<<" nST: "<<*nextSt<<" p: "<<pHatC<<"\n"; 
                }
            }
       //     getchar(); 
        }
        // fill in the normalized values into bJoint
        foreach(stIter, *cbMpPtr) {
            StatePtr currSt = stIter->first;
            BeliefPtr nStBlf = vn->pHat->getBeliefForState(currSt);
            map<StatePtr, LD, StatePtrCmp> newStateMp = nStBlf->getMp();
            LD lambda_rC =  vn->lambdaR->getPr(currSt);
            ASSERT(lambda_rC != BELIEF_NOT_FOUND);
            foreach(newStIter, newStateMp) {
                StatePtr nextSt = newStIter->first;
                LD currC = cbTmp.getPr(currSt, nextSt);
                ASSERT(currC != BELIEF_NOT_FOUND);
              //  LD pCurr = max(((double)(currC/sumC)), ((double)1e-9));
                LD pCurr = currC/sumC; 
                vn->bJointRev->setPr(nextSt, currSt, pCurr);
  //              cout<<"\tcst: "<<*currSt<<"\n\tnst: "<<*nextSt<<" pCurr: "<<pCurr<<"\n";
            }
        }
  //      getchar(); 
        map<StatePtr, Belief, StatePtrCmp> *MpPtr =   vn->bJointRev->getCbMpPtr();
        // update b_last and b_new
        foreach(stIter, *MpPtr) {
            StatePtr currNewSt = stIter->first;
            LD pOld = vn->bNew->getPr(currNewSt);
            //             ASSERT(pOld != BELIEF_NOT_FOUND);
            vn->bLast->setPr(currNewSt, pOld);
            LD pNew = 0.0; 
            map<StatePtr, LD, StatePtrCmp> mpLocal = stIter->second.getMp(); 
            foreach(innerIter, mpLocal) {
                pNew += innerIter->second; 
            }
            //LD pNew = stIter->second.getAllocPr();
            //ASSERT( (pNew >= 0.0) && (pNew <= 1.0));
            vn->bNew->setPr(currNewSt, pNew);
        }
    }
}

LD GraphPPM::getMargPr_5c(RegionPPM* vp, RegionPPM* vc, StatePtr cStOld, StatePtr cStNew) {

    LD sumR=0.0;
    map<StatePtr, Belief, StatePtrCmp> *cbMpPtr =   vp->bJointRev->getCbMpPtr();
    foreach(mpIter, *cbMpPtr) {
        StatePtr parStNew = mpIter->first;

        StatePtr redStNew =  getChildStatePtr(vp, vc, parStNew);
        // cout<<" redN: "<<*redStNew<<" ";
        LD diff = State::SSD(*cStNew, *redStNew);
        if(diff == 0.0) {
            map<StatePtr, LD, StatePtrCmp> oldStateMp = mpIter->second.getMp();
            foreach(innerIter, oldStateMp) {

                StatePtr parStOld = innerIter->first;
                StatePtr redStOld = getChildStatePtr(vp, vc, parStOld);
                LD innerDiff = State::SSD(*cStOld, *redStOld);
                //       cout<<" redStOld: "<<*redStOld<<" p: "<<innerIter->second<<"\n";
                if(innerDiff == 0.0) {
                    sumR += innerIter->second;
                } // else {
                //  cout<<"node: "<<vc->getNodeId();
                //  cout<<"\ncStOld: "<<*cStOld;
                //  cout<<"\nredStOld: "<<*redStOld;

                // }
            }
        }

    }
    //  cout<<"par: "<<vp->getNodeId()<<" vc: "<<vc->getNodeId()<<" sumR: "<<sumR<<"\n\tstO: "<<*cStOld<<"\n\tstN: "<<*cStNew<<endl;
    return sumR;
}
void GraphPPM::PPM_step5c() {
    Timer timer;
    FPRINTF(stderr, "PPM_step5c() updates the messages ");
    minMsg = maxMsg = 1.0;
    foreach(iter, adj_map) {
        RegionPPM* vn = (RegionPPM*) iter->first;
        int CR_c = vn->getCR();
        // no inedges
        if(iter->second.second.size() == 0)
            continue;
        else {
            vector<MyEdge>* inEdges = &(iter->second.second);
            foreach(edgeIter, *inEdges) {
                RegionPPM* vp = (RegionPPM*) edgeIter->ends[0];
                int CR_p =  vp->getCR();
                double CR_factor = ((double) CR_p*CR_c)/((double) (CR_p-CR_c));
                map<StatePtr, Belief, StatePtrCmp>* cbMpPtr =   vn->bJointRev->getCbMpPtr();
                foreach(mpIter, *cbMpPtr) {
                    StatePtr cStNew = mpIter->first;

                    map<StatePtr, LD, StatePtrCmp> oldStateMp = mpIter->second.getMp();
                    foreach(innerIter, oldStateMp) {
                        StatePtr cStOld = innerIter->first;
                        LD sumR = getMargPr_5c(vp, vn, cStOld, cStNew);
                        StatePtr jointSt = getJointState(cStOld, cStNew);
                        //ASSERT(innerIter->second > 0.0);
                    
                        timer.tick();
                        if((sumR > 0.0) && (innerIter->second) && (CR_p != CR_c)) {
                            double lbRc = log((double) innerIter->second);
                            double lbR = log((double) sumR);
                            LD oldMsg = edgeIter->getMsg(jointSt);
                  			// code changed here 
                            //ASSERT(oldMsg > 0.0); 
                            double lm0=-10.0; 
                            if(oldMsg > 0.0) lm0 = log((double) oldMsg);
                            double lm1 = lm0 + CR_factor*(lbR-lbRc);
                            LD newMsg = (LD) exp(lm1);
                            edgeIter->setMsg(jointSt, newMsg);
                            //      printf("sumR: %f lbR: %f lbRc: %f CRF: %f m0: %f lm0: %f lm1: %f m1: %f\n", sumR, lbR, lbRc, CR_factor, oldMsg, lm0, lm1, newMsg);
                            minMsg = (LD) min(((double) newMsg), ((double) minMsg));
                            maxMsg = (LD) max(((double) newMsg), ((double) maxMsg));

                        } else {
                            StatePtr cStOld = innerIter->first;
                            LD sumR = getMargPr_5c(vp, vn, cStOld, cStNew);
                            StatePtr jointSt = getJointState(cStOld, cStNew);
                            edgeIter->setMsg(jointSt, 1.0);  
                        }
                    }

                }

            }
        }
    }
    FPRINTF(stderr, " done\n");
}

string GraphPPM::get_status() {
    stringstream ss("");
    ss<<" "<<SR; 
    ss<<" | "<<minMsg<<" "<<maxMsg;
    ss<<" | "<<avgDiffFromPars<<" "<<maxDevFromPars<<" " <<avgMaxCaseDevFromPars;
    ss<<" | "<<avgChngInBlfs<<" "<<maxChngInBlfs<<" "<<avgMaxCaseChngInBlfs;
    
    return ss.str();
}

void GraphPPM::PPM_inner_loop() {
    clock_t c1 = clock(); 
    Timer timer; 
    PPM_step4a();
    PPM_step4b();
    PPM_step5a();
    PPM_step5b();
    PPM_step5c();
    timer.tick(); 
    clock_t c2 = clock();
    // cout<<" last iteration time: "<<timer.getElapsed()<<" c2 - c1: "<<(c2-c1)<<"\n"; 
    // getchar();  
}
/*
LD GraphPPM::getChildPr(RegionPPM* vp, RegionPPM* vn, StatePtr cSt) {
    vector<pair<StatePtr, LD> > parStates;
    RegionPPM* vp = (RegionPPM*) edgeIter->ends[0];
      map<StatePtr, LD, StatePtrCmp> childMp = vn->bNew->getMp();
    string parId = vp->getNodeId();
    map<StatePtr, LD, StatePtrCmp> parMp = vp->bNew->getMp();
    foreach(mpIter, parMp) {
        parStates.push_back(make_pair<StatePtr, LD>(mpIter->first, mpIter->second));
    }
    LD currParDiff=0.0;
    foreach(childStMpIter, childMp) {
        StatePtr cSt = childStMpIter->first;
        ANN_Cmp cmp(this, vp, vn, cSt);
        sort(parStates.begin(), parStates.end(), cmp);
        //  cout<<"child: "<<nodeId<<" par: "<<parId<<" cSt: "<<*cSt<<endl;
        bool match=true;
        // int numM=0;
        LD pMarg=0.0;
        vector<pair<StatePtr, LD> >::iterator pIter = parStates.begin();
        while(match && (pIter != parStates.end())) {
            StatePtr currSt = getChildStatePtr(vp, vn, pIter->first);
            LD diff = State::SSD(*currSt, *cSt);
            if(diff > 0.0) {
                match = false;
            } else {
                pMarg += pIter->second;
            }
            pIter++;
            //     cout<<"\t pSt: "<<*cCurrSt<<" d: "<<diff<<" p: "<<parIter->second<<"\n";
        }
        LD pChild = childStMpIter->second;
        currParDiff += abs(pMarg - pChild);
 
 
    }
}*/

void GraphPPM::updateNodesInfo() {
    FPRINTF(stderr, "updateNodesInfo() ");
    Timer timer;
    // reinitialize the graph statistics
    avgDiffFromPars=0.0;
    avgChngInBlfs=0.0;
    maxChngInBlfs=0.0;
    maxDevFromPars=0.0; 
    avgMaxCaseChngInBlfs=0.0;
    SR = 0.0; 
    avgMaxCaseDevFromPars=0.0;
    int parDevCounter=0; 
    foreach(iter, adj_map) {
        RegionPPM* vn = (RegionPPM*) iter->first;
        string nodeId = vn->getNodeId();
        vector<MyEdge> inEdges = iter->second.second;
        vn->avgDevFromPar=0.0;
        vn->maxDevFromPar=0.0;
        
        // find disagreements level with parent
        map<StatePtr, LD, StatePtrCmp> childMp = vn->bNew->getMp();
        int numM = inEdges.size();
        int cNum = childMp.size();
        LD totalDevFromPar=0.0;
        foreach(edgeIter, inEdges) {
          

            //vector<pair<StatePtr, LD> > parStates;
            RegionPPM* vp = (RegionPPM*) edgeIter->ends[0];
            string parId = vp->getNodeId();
            map<StatePtr, LD, StatePtrCmp> parMp = vp->bNew->getMp();
            /*foreach(mpIter, parMp) {
                parStates.push_back(make_pair<StatePtr, LD>(mpIter->first, mpIter->second));
                
            }*/ 
            LD currParDiff=0.0;
            foreach(childStMpIter, childMp) {
                StatePtr cSt = childStMpIter->first;
                //ANN_Cmp cmp(this, vp, vn, cSt);
                //sort(parStates.begin(), parStates.end(), cmp);
             //   cout<<"child: "<<nodeId<<" par: "<<parId<<" cSt: "<<*cSt<<" ";
                bool match=true;
                // int numM=0;
                LD pMarg=0.0;
                
                foreach(pIter, parMp) {
                    
                    
                    StatePtr currSt = getChildStatePtr(vp, vn, pIter->first);
                    LD diff = State::SSD(*currSt, *cSt);
                    if(diff > 0.0) {
                        match = false;
                    } else {
                        pMarg += pIter->second;
                    }
                   
                    //     cout<<"\t pSt: "<<*cCurrSt<<" d: "<<diff<<" p: "<<parIter->second<<"\n";
               //     cout<<*pIter->first<<" => "<<pIter->second<<" red: "<<*currSt<<" ssd: "<<diff<<" pMarg: "<<pMarg<<"\n"; 
                    
                }
                /*
                vector<pair<StatePtr, LD> >::iterator pIter = parStates.begin();
                while(match && (pIter != parStates.end())) {
                 
                }*/ 
                LD pChild = childStMpIter->second;
                LD localParDiff = abs(pMarg - pChild);
                vn->maxDevFromPar = max(localParDiff, vn->maxDevFromPar);
                currParDiff += abs(pMarg - pChild);
            //    cout<<"pMarg: "<<pMarg<<" pChild: "<<pChild<<" diff: "<<localParDiff<<" max: "<<vn->maxDevFromPar<<"\n";
            //    getchar(); 
            }
            timer.tick();
            totalDevFromPar += currParDiff;
        }
        
        if( (numM > 0) && ( cNum > 0)) {
            vn->avgDevFromPar = totalDevFromPar/((LD) (cNum*numM));
            this->maxDevFromPars = max(maxDevFromPars, vn->maxDevFromPar);
            parDevCounter++; 
            this->avgMaxCaseDevFromPars += vn->maxDevFromPar; 
            cout<<"vn: "<<nodeId<<" avgDev: "<<vn->avgDevFromPar<< " maxDev: "<<vn->maxDevFromPar<<endl; getchar();
        }
        
        LD totalBlfDiff = 0.0;
        LD maxBlfDiff = 0.0;
        foreach(stMpIter, childMp) {
            // compute belief change statistics
            LD oBlf = (LD) max( ((double)vn->bLast->getPr(stMpIter->first)), 0.0);
            LD diff = (LD) abs(((double) (oBlf- stMpIter->second)));
            if(diff > maxBlfDiff) maxBlfDiff = diff;
            totalBlfDiff += diff;
        }
        vn->avgChngInBlfs = totalBlfDiff/cNum;
        vn->maxChngInBlfs = maxBlfDiff;
        
        avgDiffFromPars += vn->avgDevFromPar;
        avgChngInBlfs += vn->avgChngInBlfs;
        maxChngInBlfs = (LD) max( ((double)vn->maxChngInBlfs), ((double) maxChngInBlfs) );
        avgMaxCaseChngInBlfs += vn->maxChngInBlfs;
        
        // update the entropy terms
        vn->SR = 0.0; 
        foreach(stMpIter, childMp) {
            LD newBelief = stMpIter->second;
            if(newBelief > 0.0) {
                vn->SR += -newBelief*log(newBelief);
            }
        }
        this->SR += vn->getCR()*vn->SR; 
        cout<<"node: "<<nodeId<<" pD: "<<vn->avgDevFromPar<<" avgC: "<<vn->avgChngInBlfs<<" maxC: "<<vn->maxChngInBlfs<<"\n"; getchar();
    }
    LD NN = (LD) adj_map.size();
    avgChngInBlfs =  avgChngInBlfs/NN;
    maxChngInBlfs =  maxChngInBlfs;
    avgMaxCaseDevFromPars = avgMaxCaseDevFromPars/((LD) parDevCounter);
    avgMaxCaseChngInBlfs =  avgMaxCaseChngInBlfs/NN;
    avgDiffFromPars =  avgDiffFromPars/NN;
    FPRINTF(stderr, " done\n");
}


void GraphPPM::PPM_LN(LD thetaLN) {
    FPRINTF(stderr, "PPM_thetaLN() ");
    if(minMsg == 0.0)
        minMsg = 1.0;
    if(maxMsg/minMsg > thetaLN) {
        FPRINTF(stderr, " updating ");
        foreach(iter, adj_map) {
            RegionPPM* vn = (RegionPPM*) iter->first;
            vector<MyEdge> *inEdges = &(iter->second.second);
            foreach(edgeIter, *inEdges) {
                foreach(stIter, edgeIter->msgMp) {
                    stIter->second = log(stIter->second) - log(minMsg) + 1.0;
                }
            }
        }
    }
    FPRINTF(stderr, "\n");
}


/*
void GraphPPM::PPM_step5c() {
    Timer timer; 
    minMsg = maxMsg = 1.0; 
    FPRINTF(stderr, "PPM_step5c() updates the messages\n");    
    int found=0, notFound=0, sf=0, snf=0; 
    foreach(iter, adj_map) {
        RegionPPM* vn = (RegionPPM*) iter->first;
        int CR = vn->getCR(); 
        vector<MyEdge> outEdges = iter->second.first;
        map<StatePtr, Belief, StatePtrCmp> *cbMpPtr = vn->pHat->getCbMpPtr();
        foreach(edgeIter, outEdges) {
            RegionPPM* vc = (RegionPPM*) edgeIter->ends[1]; 
            pair<EdgePtr, EdgePtr> edgePr = findEdge(edgeIter->ends[0], edgeIter->ends[1]);
            EdgePtr revEdge = edgePr.second;
          //  cout<<" node: "<<vn->getNodeId()<<" child: "<<vc->getNodeId()<<" epf: "<<*edgePr.first<<" epR: "<<*revEdge<<endl;
          //  cout<<" ep size: "<<edgePr.first->msgMp.size()<<" "<<edgePr.second->msgMp.size()<<endl;
          //  getchar(); 
            foreach(childMsgIter, revEdge->msgMp) {
                StatePtr childStPtr = childMsgIter->first;
                LD childMsg = childMsgIter->second; 
                LD childFac =(LD) pow(((double)childMsg), (1.0/((double) CR))); 
                map<StatePtr, Belief, StatePtrCmp> *cbMpPtr =   vn->pHat->getCbMpPtr();
                LD sumC = 0.0;
                
                
                
            //    cout<<" testSt: "<<childStPtr<<endl; 
                foreach(stIter, *cbMpPtr) {
                    StatePtr currStPar = stIter->first;                    
                    BeliefPtr nStBlf = vn->pHat->getBeliefForState(currStPar);
                    map<StatePtr, LD, StatePtrCmp> newStateMp = nStBlf->getMp();
               
                    foreach(newStIter, newStateMp) {
                        StatePtr nextStPar = newStIter->first;
                        StatePtr currSt = getChildStatePtr(vn, vc, currStPar);        
                        StatePtr nextSt = getChildStatePtr(vn, vc, nextStPar);
                        StatePtr jointSt = getJointState(currSt, nextSt);
                        bool choice= false;
                          timer.tick(); 
                    //    cout<<"currSt: "<<*currSt<<"\nnextSt: "<<*nextSt<<"\njntSt: "<<*jointSt<<" choice: "<<choice<<" sumC: "<<sumC<<"\n";          
                        if(((*jointSt) == (*childStPtr))) { // do the summing
                            choice=true; 
                            LD bPar = vn->bJointRev->getPr(nextStPar, currStPar);
                            sumC += bPar; 
                          
                     //   cout<<"\tbPar: "<<bPar<<" sumC: "<<sumC<<endl; sleep(2); 
                        }
                //      cout<<"\ttrialSt: "<<jointSt<< " added: "<<choice<<endl; 
                    }
                    
               
                }
                
                 
                // sumC = max(((double)sumC), 1e-9); 
                pair<StatePtr, StatePtr> cStPr = getChildFromMsgState(childStPtr); 
                LD pChild = vc->bJointRev->getPr(cStPr.second, cStPr.first); // reversed notation
                bool lastFound = true; 
                if(pChild == BELIEF_NOT_FOUND) {
                    notFound++;  
                    pChild = 1e-9; 
                    lastFound=false; 
                } else {
                    found++; 
                }
                ASSERT(pChild > 0.0); 
                LD newMsg;
                if(sumC > 0.0) { 
                    newMsg = pow(((double)sumC*childFac/pChild), ((double) CR));
                    sf++;
                }
                else { 
                    newMsg = 1.0; // no update
                    snf++; 
                }
                // cout<<" edge: "<<*revEdge<<" nst: "<<(*cStPr.second);
                //cout<<" childFac: "<<childFac<<" pChild: "<<pChild<<" sumC: "<<sumC<<" newM: "<<newMsg<<" found: "<<lastFound<<"\n"; 
                minMsg = (LD) min(((double) newMsg), ((double) minMsg));
                maxMsg = (LD) max(((double) newMsg), ((double) maxMsg));
            //    updateMsg(vn, vc, childStPtr, newMsg);
                pair<EdgePtr, EdgePtr> edgePair = findEdge(vn, vc);  
            //    cout<<" old msg: "<<edgePair.second->getMsg(childStPtr)<<" newMsg: "<<newMsg<<endl; getchar();  
                edgePair.second->setMsg(childStPtr, newMsg); 
            } 
        }
    }
   
    FPRINTF(stderr, "PPM_step5c() childfound: %d notFound: %d sumCfound: %d sumCnotFound: %d\n", found, notFound, sf, snf);
}
*/


