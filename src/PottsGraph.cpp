#include "PottsCB.h"
#include "PottsGraph.h"
#include "PottsModel.h"
#include "VideoNode.h"
#include <cstdio>
#include <gsl/gsl_sf_exp.h>
using namespace std;

bool PottsGraphPPM::isValid(int x, int y) {
	int N1 = model->N1;
	int N2 = model->N2;
	return ((x>=0) && (y>=0) && (x < N1) && (y< N2));
}

void PottsGraphPPM::genBETHE() {
	Timer timer;
	FPRINTF(stderr, "PottsGraphPPM::genBETHE()\n");
	int N1 = model->N1;
	int N2 = model->N2;
	fi(0, N1) {
		fj(0, N2-1) {
			stringstream ss("");
			ss<<i<<"_"<<j<<"_"<<1<<"_"<<2;
			VideoNode* vn = new VideoNode(ss.str());
			vn->setCR(1);
			addNode(vn);
			timer.tick();
		}
	}

	fi(0, N1-1) {
		fj(0, N2) {
			stringstream ss("");
			ss<<i<<"_"<<j<<"_"<<2<<"_"<<1;
			VideoNode* vn = new VideoNode(ss.str());
			vn->setCR(1);
			addNode(vn);
			timer.tick();
		}
	}

	fi(0, N1) {
		fj(0, N2) {
			stringstream ss("");
			ss<<i<<"_"<<j<<"_"<<1<<"_"<<1;
			VideoNode* vn = new VideoNode(ss.str());
			int cr = 1;
			addNode(vn);
			VideoNode* vnc = (VideoNode*) getNode(vn);
			if(isValid(i-1, j)) {
				cr--;
				string id = VideoNode::getIdFromInt( (i-1), j, 2, 1);
				VideoNode* vnp = (VideoNode*) getNode(id);
				ASSERT(vnp != 0);
				addEdge(vnp, vnc);
			}

			if(isValid(i, j-1)) {
				cr--;
				string id = VideoNode::getIdFromInt( i, (j-1), 1, 2);
				VideoNode* vnp = (VideoNode*) getNode(id);
				ASSERT(vnp != 0);
				addEdge(vnp, vnc);
			}

			if(isValid(i+1, j)) {
				cr--;
				string id = VideoNode::getIdFromInt( i, j, 2, 1);
				VideoNode* vnp = (VideoNode*) getNode(id);
				ASSERT(vnp != 0);
				addEdge(vnp, vnc);
			}

			if(isValid(i, j+1)) {
				cr--;
				string id = VideoNode::getIdFromInt( i, j, 1, 2);
				VideoNode* vnp = (VideoNode*) getNode(id);
				ASSERT(vnp != 0);
				addEdge(vnp, vnc);
			}
			vnc->setCR(cr);
			timer.tick();
			//   cout<<" vnc: "<<vnc->getNodeId()<<" cr: "<<cr<<" getCR(): "<<vnc->getCR()<<endl; getchar(); 

		}
	}
	FPRINTF(stderr, "done\n");
}

///// this function adds the requisite graph structure to the VideoGraph
//void PottsGraphPPM::genSubRegions() {
//	for (int i=0; i < model->N1; i+= cx) {
//		for (int j=0; j < model->N2; j+=cy) {
//			int wx = min(rx, model->N1-i);
//			int wy = min(ry, model->N2-j);
//			if ((wx > 0) && (wy > 0)) {
//				stringstream ss("");
//				ss<<i<<"_"<<j<<"_"<<wx<<"_"<<wy;
//				VideoNode* vn = new VideoNode(ss.str());
//				vn->setCR(1);
//				addNode(vn);
//			}
//		}
//	}
//
//	// add the next horizontal overlap region as top-level
//	for (int i=0; i < model->N1; i+= cx) {
//		for (int j=0; j < model->N2; j+= cy) {
//			int wx = min(rx, model->N1-i);
//			int wy = min(ry, model->N2-j);
//			int xn = i+cx;
//			int wxn = min(rx, model->N1-xn);
//			if ( (wxn > 0) && (wy > 0) && (wx > 0)) { // overlap exists
//				stringstream ss("");
//				int wc = i+rx-xn;
//				ss<<xn<<"_"<<j<<"_"<<wc<<"_"<<wy;
//				VideoNode* vn2 = new VideoNode(ss.str());
//				vn2->setCR(-1);
//				stringstream ssA(""), ssB("");
//				ssA<<i<<"_"<<j<<"_"<<wx<<"_"<<wy;
//				VideoNode* vnA = (VideoNode*) getNode(ssA.str());
//				//       int nxw = min(rx, N1 - xn);
//				ssB<<xn<<"_"<<j<<"_"<<wxn<<"_"<<wy;
//				VideoNode* vnB = (VideoNode*) getNode(ssB.str());
//				if ( (vnA != 0) && (vnB != 0) && VideoNode::isSubset(
//						vn2->getNodeId(), vnA->getNodeId())
//						&& VideoNode::isSubset(vn2->getNodeId(),
//								vnA->getNodeId()) ) {
//					addNode(vn2);
//					addEdge(vnA, vn2, true);
//					addEdge(vnB, vn2, true);
//				} else {
//					delete(vn2);
//				}
//			}
//		}
//	}
//	//fprintf(rgnFile, "L2_VERT\n"); //, vn->getNodeId.c_str());
//
//	// add the next vertical overlap region as top-level
//	for (int i=0; i < model->N1; i+= cx) {
//		for (int j=0; j < model->N2; j+= cy) {
//			int wx = min(rx, model->N1-i);
//			int wy = min(ry, model->N2-j);
//			int yn = j+cy;
//			int wyn = min(ry, model->N2-yn);
//			if ( (wyn > 0) && (wx > 0) && (wy > 0)) { // overlap exists
//				stringstream ss("");
//				int wc = j+ry-yn;
//				ss<<i<<"_"<<yn<<"_"<<wx<<"_"<<wc;
//				VideoNode* vn2 = new VideoNode(ss.str());
//				vn2->setCR(-1);
//				stringstream ssA(""), ssB("");
//				ssA<<i<<"_"<<j<<"_"<<wx<<"_"<<wy;
//				VideoNode* vnA = (VideoNode*) getNode(ssA.str());
//				// int nyw = min(ry, N2 - yn);
//				ssB<<i<<"_"<<yn<<"_"<<wx<<"_"<<wyn;
//				VideoNode* vnB = (VideoNode*) getNode(ssB.str());
//				if ( (vnA != 0) && (vnB != 0) && VideoNode::isSubset(
//						vn2->getNodeId(), vnA->getNodeId())
//						&& VideoNode::isSubset(vn2->getNodeId(),
//								vnA->getNodeId()) ) {
//					addNode(vn2);
//					addEdge(vnA, vn2, true);
//					addEdge(vnB, vn2, true);
//				} else {
//					delete(vn2);
//				}
//
//			}
//		}
//	}
//	//fprintf(rgnFile, "\nL3\n");
//	// add the third level regions
//	for (int i=0; i < model->N1; i += cx) {
//		for (int j=0; j < model->N2; j += cy) {
//			if (model->isValid(i-cx, j-cy) && model->isValid(i+rx, j+ry) ) { // central blocks
//				VideoNode* vn = new VideoNode(VideoNode::getIdFromInt(i, j, rx-cx, ry-cy));
//				vn->setCR(1);
//				string sA = VideoNode::getIdFromInt(i-cx, j, rx, ry-cy);
//				VideoNode* vnA = (VideoNode*) getNode(sA);
//				VideoNode* vnB = (VideoNode*) getNode(VideoNode::getIdFromInt(
//						i, j, rx-cx, ry));
//				VideoNode* vnC = (VideoNode*) getNode(VideoNode::getIdFromInt(
//						i, j, rx, ry-cy));
//				VideoNode* vnD = (VideoNode*) getNode(VideoNode::getIdFromInt(
//						i, j-cy, rx-cx, ry));
//
//				//      BeliefPtr bptr = new Belief(true);
//				//      vn->setBelief(bptr);
//
//				addNode(vn);
//				//fprintf(rgnFile, "%s\n", vn->getNodeId().c_str());
//				bool bA=false, bB=false, bC=false, bD=false;
//				if (vnA != 0) {
//					addEdge(vnA, vn, true);
//				}
//				if (vnB != 0) {
//					addEdge(vnB, vn, true);
//				}
//				if (vnC != 0) {
//					addEdge(vnC, vn, true);
//				}
//				if (vnD != 0) {
//					addEdge(vnD, vn, true);
//				}
//			}
//		}
//	}
//	//fclose(rgnFile);
//}

// modification to genSubRegions using cluster summation approach since S incorrect
// in original approach

/**
 * Modified genSubRegions - uses cluster counting approach 
 */
void PottsGraphPPM::genSubRegions() {
	map<string, int> fullList;
	map<string, int> currList;
	map<string, int> newList;
	int N1 = model->N1;
	int N2 = model->N2;
	bool filled[N1][N2];
	fi(0, N1) fj(0, N2) filled[i][j] = false;
	for (int i=0; i < model->N1; i+= cx) {
		for (int j=0; j < model->N2; j+=cy) {
			int wx = min(rx, model->N1-i);
			int wy = min(ry, model->N2-j);
			bool necessary = false;
			for (int l=i; l < i+wx; l++) {
				for (int k=j; k < j+wy; k++) {
					if (filled[l][k] == false) {
						filled[l][k] = true;
						necessary = true;
					}
				}
			}
			ASSERT( (wx >0) && (wy > 0) );
			if (necessary) {
				stringstream ss("");
				ss<<i<<"_"<<j<<"_"<<wx<<"_"<<wy;
				VideoNode* vn = new VideoNode(ss.str());
				vn->setCR(1);
				addNode(vn);
				fullList.insert(make_pair<string, int>(ss.str(), 1));
				currList.insert(make_pair<string, int>(ss.str(), 1));
			}
		}
	}
	int count = 0;

SUBLEVEL: 
//	cout<<" count: "<<count++<<" curr list: "<< endl;
//	foreach(fiter, currList) {
//		cout<<" node: "<<fiter->first<<" cr: "<<fiter->second<<endl;
//	}
	
	foreach(iter, currList) {
		map<string, int>::iterator iter2 = iter; iter2++;
		while(iter2 != currList.end()) {
			string nid1 = iter->first;
			string nid2 = iter2->first;
			string nid3 = VideoNode::intersection(nid1, nid2);

			if( (nid3 != "") && (fullList.find(nid3) == fullList.end()) && (newList.find(nid3) == newList.end()) ) {
				bool subsubregion = false;
				foreach(niter, newList) {
					if(VideoNode::isSubset(nid3, niter->first)) subsubregion = true;
				}
				if(!subsubregion) {
					int newcr = 1;
					foreach(fiter, fullList) { // compute CR
						if(VideoNode::isSubset(nid3, fiter->first)) {
							newcr = newcr - fiter->second;
						}
					}
					// add the node
					VideoNode* vn = new VideoNode(nid3);
					vn->setCR(newcr);
					addNode(vn);
					newList.insert(make_pair<string, int>(nid3, newcr));
//					cout<<" added node: "<<nid3<<" cr: "<<newcr<<endl;  
					foreach(citer, currList) { // add edges
						if(VideoNode::isSubset(nid3, citer->first)) {
							VideoNode* vp = (VideoNode*) getNode(citer->first);
							VideoNode* vc = (VideoNode*) getNode(nid3);
							ASSERT( (vp != 0) && (vc != 0) );
							addEdge(vp, vc, true);
	//						cout<<"added edge "<<vp->getNodeId()<<" to "<<vc->getNodeId()<<endl; 
						}
					}
				}
			}
			iter2++;
		}
	}
	if (newList.size() > 0) {
		fullList.insert(newList.begin(), newList.end());
		currList.clear();
		currList.insert(newList.begin(), newList.end());
		newList.clear();
		goto SUBLEVEL;
	}
}
VideoNode* PottsGraphPPM::addNode(VideoNode* vn) {
	//   cout<<"addNode() "<<vn->getNodeId()<<" CR: "<<vn->getCR()<<"\n";
	VideoNode* nn = (VideoNode*) ((MyGraph*) this)->addNode(vn);
	//ConditionalBeliefPtr cptr =(ConditionalBelief*) new PottsCB(this, vn);
	ConditionalBeliefPtr cptr = new ConditionalBelief();
	nn->pHat = cptr;
	return nn;
}

vector<INT> PottsGraphPPM::getIthStVec(int idx, int N) {
	vector<INT> cval;
	for (int j=0; j < N; j++) {
		int id = idx%model->Q;
		int vl = id;
		cval.push_back(vl);
		idx = idx/model->Q;
	}
	return cval;
}

void PottsGraphPPM::initOrigP(VideoNode* vn, bool rndm) {
	Timer timer;
	// FPRINTF(stderr, "initOrigP(%s, rndm: %d) ", vn->getNodeId().c_str(), rndm);
	Belief* tmpBlf = new Belief();
	BeliefPtr vnPtr = vn->getBeliefPtr();
	ASSERT(vn != 0);
	ASSERT(vnPtr != 0);
	vnPtr->clear();
	int N = vn->getNodeNumVar();
	LD sumP = 0.0;
	int Q = model->Q;
	int L = 1;
	fi(0, N) L*=Q;
	for (int i=0; i < L; i++) {
		timer.tick();
		vector<INT> cval = getIthStVec(i, N);
		StatePtr sptr = getStatePtr(cval);
		ASSERT(sptr->num == vn->getNodeNumVar());
		LD p = 1.0/L;
		if (rndm) {
			p = (LD) ((LD)(rand()%1000))/1000.0;
		}
		sumP += p;
		if (rndm) {
			tmpBlf->setPr(sptr, p);
		} else {
			vnPtr->setPr(sptr, p);
		}
	}
	if (rndm) {
		map<StatePtr, LD, StatePtrCmp> mp = tmpBlf->getMp();
		foreach(iter, mp) {
			vnPtr->setPr(iter->first, iter->second/sumP);

		}
	}
	/*
	 cout<<" vn: "<<vn->getNodeId()<< " rndm: "<<rndm<<" origP: \n"; 
	 map<StatePtr, LD, StatePtrCmp> mp = vn->getBeliefPtr()->getMp();
	 foreach(iter, mp) {
	 LD pSt = vn->getBeliefPtr()->getPr(iter->first);
	 cout<<"\t"<<(*iter->first)<<" => "<<pSt<<"\n"; 
	 }
	 getchar();
	 */
	delete(tmpBlf);
	//   FPRINTF(stderr, "done\n");

}

StatePtr PottsGraphPPM::getChildStatePtr(NodePtr vp, NodePtr vc,
		StatePtr parState) {
	string pid = vp->getNodeId();
	string cid = vc->getNodeId();
	vector<int> LP = VideoNode::getIntFromId(pid);
	vector<int> LC = VideoNode::getIntFromId(cid);
	vector<INT> vs;
	for (int i=LP[0]; i < (LP[0]+LP[2]); i++) {
		for (int j=LP[1]; j < (LP[1]+LP[3]); j++) {
			int idx = (i-LP[0])*LP[3]+(j-LP[1]);
			INT val = parState->getData(idx);
			if ( (i >= LC[0]) && (i < (LC[0]+LC[2])) && (j >= LC[1]) && (j
					< (LC[1]+LC[3]))) {
				vs.push_back(val);
			}
		}
	}
	ASSERT(vs.size() == vc->getNodeNumVar());
	StatePtr cState = getStatePtr(vs);
	return cState;
}

LD PottsGraphPPM::getStatePr(VideoNode* vn, StatePtr sptr) {
	vector<int> L = VideoNode::getIntFromId(vn->getNodeId());
	LD result = 1.0;
	// cout<<"getStatePr() : "<<vn->getNodeId()<<" state: "<<*sptr<<"\n";
	for (int i=0; i < L[2]; i++) {
		for (int j=0; j < L[3]; j++) {
			int x = i+L[0];
			int y = j+L[1];
			int idx = i*L[3]+j;
			stringstream ssf("");
			ssf<<x<<"."<<y;
			string csid = ssf.str();
			INT cval = sptr->getData(idx);
			StatePtr svPtr = getStatePtr(cval);
			LD sfac = model->getP(csid, svPtr);
			result *= sfac;
			//    cout<<"\t"<<csid<<" => "<<sfac<<" res: "<<result<<"\n"; getchar();

			// next horizontal bond
			if (i < (L[2]-1)) {
				int xn = x+1;
				int yn = y;
				int idn = (i+1)*L[3]+j;
				stringstream ssp("");
				ssp<<x<<"."<<y<<"_"<<xn<<"."<<yn;
				string cpid = ssp.str();
				INT nval = sptr->getData(idn);
				vector<INT> vp;
				vp.push_back(cval);
				vp.push_back(nval);
				StatePtr spPtr = getStatePtr(vp);
				LD sPr = model->getP(cpid, spPtr);

				result *= sPr;
				//     cout<<"\t"<<cpid<<" => "<<sPr<<" res: "<<result<<"\n"; getchar();
			}
			// next vertical bond
			if (j < (L[3]-1)) {
				int xn = x;
				int yn = y+1;
				int idn = (i)*L[3]+(j+1);
				stringstream ssp("");
				ssp<<x<<"."<<y<<"_"<<xn<<"."<<yn;
				string cpid = ssp.str();
				INT nval = sptr->getData(idn);
				vector<INT> vp;
				vp.push_back(cval);
				vp.push_back(nval);
				StatePtr spPtr = getStatePtr(vp);
				LD sPr = model->getP(cpid, spPtr);
				result *= sPr;
				//     cout<<"\t"<<cpid<<" => "<<sPr<<" res: "<<result<<"\n"; getchar();
			}

		}
	}
	return result;
}

LD PottsGraphPPM::get_p_hat(VideoNode* vn, StatePtr si, StatePtr so, LD theta,
		LD dt) {
	int N = vn->getNodeNumVar();
	int Q = model->Q;
	int Ine= 0, Ieq=0;
	if (theta*dt != 0.5) {
		fi(0, N) { // compute the non-equal part
			if(si->getData(i) != so->getData(i))
			Ine++;
			else
			Ieq++;
		}
	}
	LD fac1 = 1.0, fac2 = 1.0;
	if ( (Ine != 0) || (Ieq != 0)) {
		fac1 = (LD) pow(((double) (theta*dt)), ((double) (Ine)));
		fac2 = (LD) pow(((double) (1.0 - theta*dt)), ((double) (Ieq)));
	}
	LD pSo = getStatePr(vn, so);
	LD pSi = getStatePr(vn, si);
	// cout<<"p_hat(): "<<"si: "<<*si<<" so: "<<*so<<" Ine: "<<Ine<<" Ieq: "<<Ieq<<" pSi: "<<pSi<<" pSo: "<<pSo<<" fac1:"<<fac1<< " fac2: "<<fac2<<"\n"; getchar();

	LD result = 0.0;
	if (pSi > 0.0) {
		result = (pSo/pSi)*(fac1)*(fac2);
	}
	return result;
}

void PottsGraphPPM::initCondP(VideoNode* vn, LD theta, LD dt) {
	Timer timer;
	int N = vn->getNodeNumVar();
	vn->pHat->clear();
	int L = 1;
	int Q = model->Q;
	fi(0, N) L*=Q;
	FPRINTF(stderr, "initCondP(%s, theta=%g, dt=%g) num_states: %d", vn->getNodeId().c_str(), theta, dt, (L*L));
	LD sumP = 0.0;
	fi(0, L) {
		fj(0, L) {
			timer.tick();
			vector<INT> vi = getIthStVec(i, N);
			vector<INT> vo = getIthStVec(j, N);
			StatePtr si = getStatePtr(vi);
			StatePtr so = getStatePtr(vo);
			LD pCurr = get_p_hat(vn, si, so, theta, dt);
			//       cout<<" si: "<<*si<<" so: "<<*so<<" p: "<<pCurr<<" sumP: "<<sumP<<endl;   
			vn->pHat->setPr(si, so, pCurr);
			sumP += pCurr;
		}
	}
	//  getchar(); 
	// normalize the messages
	fi(0, L) {
		fj(0, L) {
			timer.tick();
			vector<INT> vi = getIthStVec(i, N);
			vector<INT> vo = getIthStVec(j, N);
			StatePtr si = getStatePtr(vi);
			StatePtr so = getStatePtr(vo);
			LD pCurr=0.0;
			LD pTmp = vn->pHat->getPr(si, so);
			if(sumP > 0.0) {
				pCurr = pTmp/sumP;
			}
			// cout<<"node: "<<vn->getNodeId()<<" ["<<*si<<"] -> ["<<*so<<"] = "<<pCurr<<"\n"; getchar();
			vn->pHat->setPr(si, so, pCurr);

		}
	}
	/*   cout<<"vn: "<<vn->getNodeId()<<" condBlf: "; 
	 map<StatePtr, Belief, StatePtrCmp>* cbMpPtr = vn->pHat->getCbMpPtr();  
	 foreach(cbIter, *cbMpPtr) {
	 StatePtr sCurr =  cbIter->first; 
	 BeliefPtr bPtr = &(cbIter->second);
	 map<StatePtr, LD, StatePtrCmp> mp2 = bPtr->getMp();  
	 foreach(bpIter, mp2) {
	 cout<<" si: "<<*(bpIter->first)<<" so: "<<*sCurr<< " p: "<<bpIter->second<<endl; 
	 }
	 }
	 getchar(); */
	FPRINTF(stderr, " done\n");
}

void PottsGraphPPM::PPM_algo(LD theta, LD dt) {
	FPRINTF(stderr, "PPM_algo(theta = %g,dt =  %g) ", theta, dt);
	bool rndmInit=false;
	FILE* file = fopen("POTTS.log", "a+");
	fprintf(
			file,
			"##############################################################################################\n");
	fprintf(
			file,
			"# ITER  -lnZ_approx HR [SR  minM maxM <ParDiff> PD_max <PD_max> <Cblf> blfC_max <blfC_max>]\n");
	fprintf(
			file,
			"##############################################################################################\n");
	int count=0;
	STEP_2: foreach(iter, adj_map) {
		VideoNode* vn = (VideoNode*) iter->first;
		initOrigP(vn, rndmInit);
		initCondP(vn, theta, dt);
	}
	STEP_3: PPM_step3();
	INNER_LOOP: for (int i=0; i < 10; i++) {
		PPM_inner_loop();
		updateNodesInfo();
		string stats = get_status();
		LD HR = getHR();
		LD lnZR = HR-SR; // SR is a class variable
		fprintf(file, "%d: %g %g [%s]\n", count, lnZR, HR, stats.c_str() );
		fflush(file);
		count++;
		cout<<count<<": -lnZ: " <<lnZR<<" HR: "<<HR<<stats<<"\n";
		getchar();
	}
	fclose(file);
}

LD PottsGraphPPM::getHR() {
	LD result = 0.0;
	foreach(np, adj_map) {
		VideoNode* vn = (VideoNode*) np->first;
		map<StatePtr, LD, StatePtrCmp> mp = vn->bNew->getMp();
		LD HR = 0.0;
		foreach(iter, mp) {
			StatePtr sptr = iter->first;
			LD blfSt = iter->second;
			LD fSt = getStatePr(vn, sptr);
			LD hSt = - log(fSt);
			HR += blfSt*hSt;
		}
		result += vn->getCR()*HR;
		//cout<<" node: "<<vn->getNodeId()<<" HR: "<<HR<<"\n";
	}
	return result;

}

void PottsGraphPPM::PPM_algo2(LD theta, LD dt) {
	FPRINTF(stderr, "PPM_algo(theta = %g,dt =  %g) ", theta, dt);
	bool rndmInit=true;
	FILE* file = fopen("potts_algo2.log", "w");
	fprintf(file, "# model: %s\n", model->getStats().c_str());
	fprintf(file, "# theta: %f dt: %f\n", theta, dt);
	fprintf(file,"##############################################################################################\n");
	fprintf(file,"# ITER  -lnZ_approx HR [SR  minM maxM <ParDiff> PD_max <PD_max> <Cblf> blfC_max <blfC_max>]\n");
	fprintf(file,"##############################################################################################\n");
	int count=0;
	stringstream nss("");
	nss<<"0_0_"<<rx<<"_"<<ry;
	NodePtr np = getNode(nss.str());
	ASSERT(np != 0);
	VideoNode* n00 = (VideoNode*) np;
	FILE* file2 = fopen("b00.log", "w");
	fprintf(file2, "# theta: %f dt: %f\n", theta, dt);
	fprintf(file2, "# time iter b00\n");

	STEP_2:
	//  initFG();
	foreach(iter, adj_map) {
		VideoNode* vn = (VideoNode*) iter->first;
		initOrigP(vn, rndmInit);
		initCondP(vn, theta, dt);
	}

	//    fg->bp_synch(false);

	STEP_3: PPM_step3();
	INNER_LOOP: PPM_inner_loop();
	updateNodesInfo();

	LD b00 = getVarMargPr(n00, 0, 0, 1);
	LD HR = getHR();
	LD lnZR = HR-SR;

	fprintf(stderr, "%f %d b: %f -lnZ: %f\n", currTime, count, b00,lnZR);
	count++;

	fflush(file2);
	cout<<get_status()<<"\n";
	fprintf(file, "t: %f iter: %d %s\n", currTime, count, get_status().c_str());
	fflush(file);
	bool sc = stoppingCriterion(count);
	//   printf(" maxBlf: %f maxPD: %f count: %f sc: %d\n", maxChngInBlfs, maxDevFromPars, count, sc); getchar();
	if ( (maxDevFromPars < 0.05) || (count > 5) || (maxChngInBlfs < 0.02)) {

		LD ct = updateTime(dt);
		fprintf(file2, "%f %d %f %f\n", currTime, count, b00, lnZR);
		count = 0;
		if (ct > 5.0) {
			fclose(file);
			fclose(file2);
			return;
		} else {
			goto STEP_3;
		}
	} else {

		goto INNER_LOOP;
	}

}

void PottsGraphPPM::initFG() {
	int Q = model->Q;
	vector<INT> varVals;
	for (int i=0; i < Q; i++) {
		varVals.push_back((INT) i);
	}
	fg->setVarVals(varVals);
	LD pEach = 1.0/((LD) Q);
	for (int i=0; i < model->N1; i++) {
		for (int j=0; j < model->N2; j++) {
			// set singleton nodes
			stringstream ssp(""), ssf(""), ssn("");
			ssp<<i<<"."<<j;
			ssn<<"X"<<i<<"."<<j;
			Belief* bptr = new Belief();
			LD sumP = 0.0;
			foreach(valIter, varVals) {
				StatePtr spf = fg->getStatePtr(*valIter);
				//                cout<< ssp.str()<<" state: "<<*spf<<"\n"; getchar();
				bptr->setPr(spf, pEach);
			}
			//            DynNodePtr xn = new DynNode(ssn.str(), 1, bptr);
			//            fg->addNode(xn);
		}
	}
	for (int i=0; i < model->N1; i++) {
		for (int j=0; j < model->N2; j++) {
			// set singleton factors
			stringstream ssp(""), ssf(""), ssn("");
			ssp<<i<<"."<<j;
			ssn<<"X"<<i<<"."<<j;
			ssf<<"F_X"<<i<<"."<<j;
			string fid = ssf.str();
			string nid = ssn.str();
			string pid = ssp.str();
			Belief* bptr = new Belief();
			LD sumP = 0.0;
			foreach(valIter, varVals) {
				StatePtr spf = fg->getStatePtr(*valIter);
				//  cout<< ssp.str()<<" state: "<<*spf<<"\n"; getchar();
				LD pVal = model->getP(ssp.str(), spf);

				sumP += pVal;
			}

			foreach(valIter, varVals) {
				StatePtr spf = fg->getStatePtr(*valIter);
				//   cout<< ssp.str()<<" state: "<<*spf<<"\n"; getchar();
				LD pVal = model->getP(ssp.str(), spf);
				fg->getNode(ssf.str())->getBeliefPtr()->setPr(spf, pVal/sumP);
				//                bptr->setPr(spf, pVal/sumP);
			}
			//            NodePtr xn = fg->getNode(ssn.str());
			//            DynNodePtr fn = new DynNode(nid, 1, bptr);
			//            fg->addNode(fn);
			//            fg->addEdge(fn, xn);
		}
	}

	for (int i=0; i < model->N1; i++) {
		for (int j=0; j < model->N2; j++) {
			if (j < (model->N2 - 1)) {
				// set singleton factors
				stringstream ssp(""), ssf(""), ssn1(""), ssn2("");
				ssp<<i<<"."<<j<<"_"<<i<<"."<<(j+1);
				ssn1<<"X"<<i<<"."<<j;
				ssn2<<"X"<<i<<"."<<(j+1);
				ssf<<"F_X"<<i<<"."<<j<<"_"<<"X"<<i<<"."<<(j+1);
				string fid = ssf.str();
				string nid = ssn1.str();
				string pid = ssp.str();
				Belief* bptr = new Belief();
				LD sumP = 0.0;
				foreach(valIter, varVals) {
					foreach(valIter2, varVals) {
						vector<INT> sv;
						sv.push_back(*valIter);
						sv.push_back(*valIter2);
						StatePtr spf = fg->getStatePtr(sv);
						//    cout<< ssp.str()<<" state: "<<*spf<<"\n"; getchar();
						LD pVal = model->getP(ssp.str(), spf);
						sumP += pVal;
					}
				}

				foreach(valIter, varVals) {
					foreach(valIter2, varVals) {
						vector<INT> sv;
						sv.push_back(*valIter);
						sv.push_back(*valIter2);
						StatePtr spf = fg->getStatePtr(sv);
						//   cout<< ssp.str()<<" state: "<<*spf<<"\n"; getchar();
						LD pVal = model->getP(ssp.str(), spf);
						fg->getNode(ssf.str())->getBeliefPtr()->setPr(spf, pVal/sumP);
						//                bptr->setPr(spf, pVal/sumP);
					}
				}

				//            NodePtr xn1 = fg->getNode(ssn1.str());
				//            NodePtr xn2 = fg->getNode(ssn2.str());
				//            DynNodePtr fn = new DynNode(fid, 2, bptr);
				//            fg->addNode(fn);
				//            fg->addEdge(fn, xn1);
				//            fg->addEdge(fn, xn2);
			}
		}
	}
	for (int i=0; i < model->N1; i++) {
		for (int j=0; j < model->N2; j++) {
			if (i < (model->N1 - 1)) {
				// set singleton factors
				stringstream ssp(""), ssf(""), ssn1(""), ssn2("");
				ssp<<i<<"."<<j<<"_"<<(i+1)<<"."<<j;
				ssn1<<"X"<<i<<"."<<j;
				ssn2<<"X"<<(i+1)<<"."<<j;
				ssf<<"F_X"<<i<<"."<<j<<"_"<<"X"<<(i+1)<<"."<<j;
				string fid = ssf.str();
				string nid = ssn1.str();
				string pid = ssp.str();
				Belief* bptr = new Belief();
				LD sumP = 0.0;
				foreach(valIter, varVals) {
					foreach(valIter2, varVals) {
						vector<INT> sv;
						sv.push_back(*valIter);
						sv.push_back(*valIter2);

						StatePtr spf = fg->getStatePtr(sv);
						//      cout<< ssp.str()<<" state: "<<*spf<<"\n"; getchar();
						LD pVal = model->getP(ssp.str(), spf);
						sumP += pVal;
					}
				}

				foreach(valIter, varVals) {
					foreach(valIter2, varVals) {
						vector<INT> sv;
						sv.push_back(*valIter);
						sv.push_back(*valIter2);
						StatePtr spf = fg->getStatePtr(sv);
						//   cout<< ssp.str()<<" state: "<<*spf<<"\n"; getchar();
						LD pVal = model->getP(ssp.str(), spf);
						fg->getNode(ssf.str())->getBeliefPtr()->setPr(spf, pVal/sumP);
						//                bptr->setPr(spf, pVal/sumP);
					}
				}

				//            NodePtr xn1 = fg->getNode(ssn1.str());
				//            NodePtr xn2 = fg->getNode(ssn2.str());
				//            DynNodePtr fn = new DynNode(fid, 2, bptr);
				//            fg->addNode(fn);
				//      fg->addEdge(fn, xn1);
				//      fg->addEdge(fn, xn2);
			}
		}
	}
}

LD PottsGraphPPM::getVarMargPr(VideoNode* vn, int x, int y, INT val) {
	LD res = 0.0;
	map<StatePtr, LD, StatePtrCmp> mp = vn->getBeliefPtr()->getMp();
	vector<INT> L = VideoNode::getIntFromId(vn->getNodeId());
	int idx = (x-L[0])*L[3] + (y-L[1]);
	foreach(mpIter, mp) {
		StatePtr sptr = mpIter->first;
		if(sptr->getData(idx) == val) {
			res += mpIter->second;
		}
	}
	return res;
}

bool PottsGraphPPM::stoppingCriterion(int iter, LD PD_thres, LD blf_thres,
		int iter_thres) {
	if ( (maxChngInBlfs < blf_thres)|| (maxDevFromPars < PD_thres ) || (iter
			> iter_thres))
		return true;
	else
		return false;
}

LD PottsGraphPPM::updateTime(LD currDt) {
	currTime += currDt;
	foreach(iter, adj_map) {
		VideoNode* vn = (VideoNode*) iter->first;
		vn->getBeliefPtr()->clear();
		map<StatePtr, LD, StatePtrCmp> mp = vn->bNew->getMp();
		foreach(mpIter, mp) {
			StatePtr sptr = mpIter->first;
			LD pval = mpIter->second;
			vn->getBeliefPtr()->setPr(sptr, pval);

		}
	}
	return currTime;
}

void PottsGraphPPM::singleIter(LD theta_ppm, LD dt_ppm) {
	Timer timer;
	FPRINTF(stderr, "PottsGraphPPM::singleIter()\n");
	//PottsCB::setTheta(theta_ppm);
	//PottsCB::setDt(dt_ppm);
	foreach(iter, adj_map) {
		VideoNode* vn = (VideoNode*) iter->first;
		initCondP(vn, theta_ppm, dt_ppm);
		timer.tick();
	}

	MSG_STEP: PPM_step3();
	timer.tick();
	LOOP_STEP:

	for (int i=0; i < 2; i++) {
		PPM_inner_loop();
		// updateNodesInfo();
		timer.tick();
	}
	updateTime(dt_ppm);
	FPRINTF(stderr, "done\n");
}

vector<vector<INT> > PottsGraphPPM::greedyDecision() {
	vector<vector<INT> > decision;
	int N1 = model->N1;
	int N2 = model->N2;
	INT image[N1][N2];
	fi(0, N1) fj(0, N2) image[i][j] = 0;
	foreach(iter, adj_map) {
		VideoNode* vn = (VideoNode*) iter->first;
		if(vn->getNodeNumVar() == 1) {
			LD pMax=0.0;
			vector<INT> L = VideoNode::getIntFromId(vn->getNodeId());
			int x= L[0], y = L[1];
			map<StatePtr, LD, StatePtrCmp> mp = vn->bNew->getMp();
			foreach(mpIter, mp) {
				INT currPxl = mpIter->first->getData(0);
				LD pCurr = mpIter->second;
				if(pCurr > pMax) {
					pMax = pCurr;
					image[x][y] = currPxl;
				}
			}
		}
	}
	fi(0, N1) {
		vector<INT> cv;
		fj(0, N2) {
			cv.push_back(image[i][j]);
		}
		decision.push_back(cv);
	}
	return decision;
}

void PottsGraphPPM::initAllNodesOrigP(bool rndm) {
	foreach(iter, adj_map) {
		VideoNode* vn = (VideoNode*) iter->first;
		initOrigP(vn, rndm);
	}
}

void PottsGraphPPM::PPM_step2(LD theta, LD delta_t) {
	FPRINTF(stderr, "PPM_step2() ");
	Timer timer;
	foreach(iter, adj_map) {
		VideoNode* vn = (VideoNode*) iter->first;
		initCondP(vn, theta, delta_t);
		timer.tick();
	}
	FPRINTF(stderr, "done\n");
}

LD PottsGraphPPM::getGraphStBlf(int i) {
	int Q = model->Q;
	int N1 = model->N1;
	int N2 = model->N2;
	vector<INT> sCurr = model->getIthStVec(i, N1*N2);
	//	cout<<"getGraphStBlf("<<i<<") "<<getStr<INT>(sCurr)<<endl; 
	LD result = 1.0;
	foreach(iter, adj_map) {
		VideoNode* vn = (VideoNode*) iter->first;
		BeliefPtr bptr = vn->getBeliefPtr();
		vector<int> L = VideoNode::getIntFromId(vn->getNodeId());
		vector<INT> childStVec;
		for(int i=L[0]; i < L[0]+L[2]; i++) {
			for(int j=L[1]; j < L[1]+ L[3]; j++) {
				int idx = i*N2 + j;
				INT val = sCurr[idx];
				childStVec.push_back(val);
			}
		}
		ASSERT(childStVec.size() == vn->getNodeNumVar());
		State* sChild = new State(childStVec);
		LD pCurr = bptr->getPr(sChild);
		//		cout<<"vn : "<< vn->getNodeId()<<" state: "<<*sChild<<" p: "<<pCurr<<endl; 
		//		getchar(); 
		ASSERT(pCurr != BELIEF_NOT_FOUND);
		ASSERT(pCurr > 0.0);
		result = result * ((LD) pow((double) pCurr, (double) vn->getCR()));
		delete(sChild);
	}
	return result;
}

vector<LD> PottsGraphPPM::getNormGraphBlfs() {
	int Q = model->Q;
	int N1 = model->N1;
	int N2 = model->N2;
	//foreach(iter, adj_map) {
	//	VideoNode* vn = (VideoNode*) iter->first;
	//	cout<<vn->getNodeId()<<" => "<<vn->getCR()<<endl;
	//}
	
	int numStates = (int) MyDef::pow(Q, N1*N2);
	vector<LD> result;
	LD sumP = 0.0;
	for (int i=0; i < numStates; i++) {
		LD p = getGraphStBlf(i);
		sumP += p;
		result.push_back(p);
		// cout<<" st: "<<i<<" p: "<<p<<endl ; 
	}
	ASSERT(sumP > 0.0);
	for (int i=0; i < numStates; i++) {
		result[i] = result[i]/sumP;
	}
	return result;
}
/**
 * \arg modelFileName the AMPL model file generated
 * \arg dataFileName the AMPL data file generated
 * \arg commandFileName the file containing AMPL execution code
 *  \note The AMPL code has the following parameters:  	
 * \n	params: q, c_alpha(rgn), numVar(rgn), b0_rgn(1..q^n) , pc_rgn(1..q^n)
 * \n	vars: b2_rgn(1..q^n, 1..q^n)
 *  
 */ 
void PottsGraphPPM::getAMPL(string modelFileName, string dataFileName, string commandFileName) {
	FILE* modelFile = fopen(modelFileName.c_str(), "w"); 
	FILE* dataFile = fopen(modelFileName.c_str(), "w"); 
	FILE* amplFile = fopen(modelFileName.c_str(), "w"); 
	writeModelFile(modelFile); 
	//writeDataFile(dataFile); 
	//writeCommandFile(commandFile); 
	fclose(modelFile); 
	fclose(dataFile); 
	fclose(amplFile); 
}


/**
 * \arg modelFile
 */
 void PottsGraphPPM::writeModelFile(FILE* modelFile) {
 	fprintf(modelFile, "####################################################################\n");  
 	fprintf(modelFile, "# Q = The number of values variable can take (implemented directly)\n"); 
 	fprintf(modelFile, "# R = The set of regions\n"); 
 	fprintf(modelFile, "# cr= The counting number for each region\n"); 	
  	fprintf(modelFile, "####################################################################\n");  
 	// fprintf(modelFile, "param Q;\n");
 	// fprintf(modelFile, "set R;\n");
 	// fprintf(modelFile, "param cr {r in R};\n\n");
 	fprintf(modelFile, "####################################################################\n");  
 	fprintf(modelFile, "# The singleton original beliefs for each region\n");
 	fprintf(modelFile, "#\tb1_rgn = original state probability for region\n"); 
 	fprintf(modelFile, "#\tpc2_rgn = conditional p(t+delta t | t) probability for region\n"); 
 	fprintf(modelFile, "#\tb2_rgn = dual probability for region (variable)\n");
 	fprintf(modelFile, "####################################################################\n");
 	// total number of nodes  
 	int totNumOfNodes = 0; 
 	// the number of values each node can take 
 	int Q = this->model->Q;
 	foreach(iter, adj_map) { // iterate over each node 
 		NodePtr node = iter->first; 
 		string nodeName = node->getNodeId();
 		int nodeVar = node->getNodeNumVar();
 		int numStates = MyDef::pow(Q, nodeVar); 
 		fprintf(modelFile, "param cr_%s;\n", nodeName.c_str());
 		fprintf(modelFile, "param b1_%s {i in 1..%d};\n", nodeName.c_str(), numStates);
 		fprintf(modelFile, "param pc2_%s {i in 1..%d, j in 1..%d};\n", nodeName.c_str(), numStates, numStates);
 		fprintf(modelFile, "var b2_%s {i in 1..%d, j in 1..%d};\n", nodeName.c_str(), numStates, numStates);
 		fprintf(modelFile, "\n");
 		totNumOfNodes++;  
 	}
 	fprintf(modelFile, "\n"); 
 	fprintf(modelFile, "####################################################################\n");
 	fprintf(modelFile, "# Writing the cost\n"); 
 	fprintf(modelFile, "####################################################################\n");
 	fprintf(modelFile, "minimize VarFreeEnergy:\n"); 
 	foreach(iter, adj_map) {
 		RegionPPM* node = (RegionPPM*) iter->first;
 		string nodeName= node->getNodeId();
 		int nodeVar = node->getNodeNumVar();
 		int numStates = MyDef::pow(Q, nodeVar); 
 		fprintf(modelFile, "+sum{i in 1..%d, j in 1..%d} (",  numStates, numStates); 
 		fprintf(modelFile, " cr_%s * b2_%s[i, j] * ( log(b2_%s[i, j]) - log(b1_%s[i]) - log(pc2_%s[i,j]) ) ", 
 			nodeName.c_str(), nodeName.c_str(), nodeName.c_str(), nodeName.c_str(), nodeName.c_str()); 
 		fprintf(modelFile, ")\n"); 
 	}
 	fprintf(modelFile, ";\n\n"); 
 	
 
 	writeNormConstraints(modelFile); 
 }
 
 /**
  * \arg modelFile the AMPL model file which is being updated 
  * \note this function is to be modified in the GBP implementation
  */
 void PottsGraphPPM::writeNormConstraints(FILE* modelFile) {
 	int Q = this->model->Q;
 	fprintf(modelFile, "####################################################################\n");
 	fprintf(modelFile, "# Writing the normalization constraint \n"); 
 	fprintf(modelFile, "####################################################################\n");
 	foreach(iter, adj_map) {
 		RegionPPM* node = (RegionPPM*) iter->first;
 		string nodeName= node->getNodeId();
 		int nodeVar = node->getNodeNumVar();
 		int numStates = MyDef::pow(Q, nodeVar);
 		fprintf(modelFile, "s.t. norm_%s { i in 1..%d} : sum{ j in 1..%d } b2_%s[i,j] = b1_%s[i];\n", 
 			nodeName.c_str(), numStates, numStates, nodeName.c_str(), nodeName.c_str()
 			); 
 	}
 }
 


/*
void GraphMVI::initCondP(VideoNode* vn, LD theta, LD dt) {
	Timer timer;
	int N = vn->getNodeNumVar();
	vn->pHat->clear();
	int L = 1;
	int Q = model->Q;
	fi(0, N) L*=Q;
	FPRINTF(stderr, "initCondP(%s, theta=%g, dt=%g) num_states: %d", vn->getNodeId().c_str(), theta, dt, (L*L));
	LD sumP = 0.0;
	fi(0, L) {
		fj(0, L) {
			timer.tick();
			vector<INT> vi = getIthStVec(i, N);
			vector<INT> vo = getIthStVec(j, N);
			StatePtr si = getStatePtr(vi);
			StatePtr so = getStatePtr(vo);
			LD pCurr = get_p_hat(vn, si, so, theta, dt);
			//       cout<<" si: "<<*si<<" so: "<<*so<<" p: "<<pCurr<<" sumP: "<<sumP<<endl;   
			vn->pHat->setPr(si, so, pCurr);
			sumP += pCurr;
		}
	}
}
*/ 

//void GraphMVI::initOrigP(VideoNode* vn, bool rndm=true) {}
