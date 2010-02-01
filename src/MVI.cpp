#include "MVI.h" 
using namespace std; 

void GraphMVI::initAllNodesOrigP(bool rndm) {
	cout<<"GraphMVI::initAllNodesOrigP()\n"; 
	LD epsilon_zero = 0.001;
	foreach(iter, adj_map) {
		VideoNode* vn = (VideoNode*) iter->first; 
		int n = vn->getNodeNumVar(); 
		int ns = (int) MyDef::pow(model->Q, n); 
		LD p0 = 1.0 - (ns-1)*epsilon_zero; 
		fi(0, ns) {
			vector<INT> cval = getIthStVec(i, n);
			StatePtr sptr = getStatePtr(cval);
			if(i==0) {
				vn->getBeliefPtr()->setPr(sptr, p0); 
			} else {
				vn->getBeliefPtr()->setPr(sptr, epsilon_zero); 
			}
		} 
	}
	//getchar(); 
}

void GraphMVI::initCondP(VideoNode* vn, LD theta, LD dt) {
	//cout<<"GraphMVI::initCondP("<<vn->getNodeId()<<")\n"; 
	Timer timer;
	int N = vn->getNodeNumVar();
	vn->pHat->clear();
	int L = 1;
	int Q = model->Q;
	fi(0, N) L*=Q;
	//FPRINTF(stderr, "initCondP(%s, theta=%g, dt=%g) num_states: %d", vn->getNodeId().c_str(), theta, dt, (L*L));
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

	//FPRINTF(stderr, " done\n");
}


LD GraphMVI::get_p_t(int i, int j, int xi, int xo) {
	INT e = evidence[i][j];
	if(e == 0) { // difference wrt previous frame = 0
	if( xo == (INT) max((int) 0, (int) (xi-1)) ) {
	//	if(xo == xi) {
			return theta_t;
		} else {
			return epsilon_t; 
		}
	} else {
		if(xo == (C-1) ) {
			return (1.0-(C-1)*epsilon_zero);
		} else {
			return epsilon_zero; 
			}
		}	
}

LD GraphMVI::get_p_hat(VideoNode* vn, StatePtr si, StatePtr so, LD theta, LD dt) {
	//	printf("MVI::get_p_hat(%s)\n", vn->getNodeId().c_str()); 
		LD result; 
		if(vn->getNodeNumVar() == 1) { // model the temporal aspect
			vector<INT> idv = VideoNode::getIntFromId(vn->getNodeId() ); 
			// int idx = idv[0]*N2 + idv[1]; // corresponding to the evidence
			int i = idv[0], j = idv[1];
			INT xo = so->getData(0); 
			INT xi = si->getData(0);   
			return get_p_t(i, j, xi, xo); 				
		} else { // model the spatial characteristic
			INT xl = so->getData(0); 
			INT xr = so->getData(1);
			LD px = 1.0;  
			if(xl == xr) {
				 px = theta_r;
			} else {
				px = epsilon_r;
			}
			vector<INT> L = VideoNode::getIntFromId(vn->getNodeId());
			int i1 = L[0], j1 = L[1];
			int i2 = L[0], j2 = L[1]; 
			if(L[0] == 2) { i2++; }
			else { j2++; }
			px *= get_p_t(i1, j1, si->getData(0), so->getData(0));
			px *= get_p_t(i2, j2, si->getData(1), so->getData(1));
			return px; 
		}	
	}
	

vector< vector<INT> > GraphMVI::mmse_estimate() {
		cout<<"MVI::map_estimate()"<<endl; 
		vector<INT> row(model->N2, 0); 
		int N = model->N1*model->N2; 
		vector< vector<INT> > result; fi(0, model->N1) result.push_back(row);  
		foreach(iter, adj_map) {
			VideoNode* vn = (VideoNode*) iter->first;
			
			if(vn->getNodeNumVar() == 1) {
				vector<int> L = VideoNode::getIntFromId(vn->getNodeId());
				int idx = L[0]*N2 + L[1]; 
				LD xMMSE = 0.0;	
				LD pMap = 0.0;  
				for(int x=0; x < C; x++) { 
					StatePtr s = new State(1, &x); 
					LD pCurr = vn->getBeliefPtr()->getPr(s);
					if(pCurr > pMap) {
						xMMSE = ((LD) x)*pCurr;
						pMap = pCurr; 
					}
				}
				int i= L[0], j= L[1]; 
				if(xMMSE > ((LD) C)/2.0) {
					result[i][j] = C-1;
				} else {
					result[i][j] = 0; 
				}
			}
		}
		
		//vector<vector<INT> > decision; fi(0, model->N1) decision.push_back(row);  
		//fi(0, model->N1) fj(0, model->N2) { if(result[i][j] == (C-1) ) decision[i][j] = C-1; }
		return result; 
	}
	
	
vector< vector<INT> > GraphMVI::map_estimate() {
		cout<<"MVI::map_estimate()"<<endl; 
		vector<INT> row(model->N2, 0); 
		int N = model->N1*model->N2; 
		vector< vector<INT> > result; fi(0, model->N1) result.push_back(row);  
		foreach(iter, adj_map) {
			VideoNode* vn = (VideoNode*) iter->first;
			
			if(vn->getNodeNumVar() == 1) {
				vector<int> L = VideoNode::getIntFromId(vn->getNodeId());
				int idx = L[0]*N2 + L[1]; 
				INT xMap = 0;	
				LD pMap = 0.0;  
				for(int x=0; x < C; x++) { 
					StatePtr s = new State(1, &x); 
					LD pCurr = vn->getBeliefPtr()->getPr(s);
					if(pCurr > pMap) {
						xMap = x;
						pMap = pCurr; 
					}
				}
				int i= L[0], j= L[1]; 
				result[i][j] = xMap;
			}
		}
		
		vector<vector<INT> > decision; fi(0, model->N1) decision.push_back(row);  
		fi(0, model->N1) fj(0, model->N2) { if(result[i][j] == (C-1) ) decision[i][j] = (C-1); }
		return decision; 
	}
	
	
LD GraphMVI::updateTime(LD dt) {
		cout<<"MV::update time("<<dt<<")\n"; 
		LD psum = 0.0;
		currTime += dt;  
		foreach(iter, adj_map) {
			VideoNode* vn = (VideoNode*) iter->first;
		//	vn->getBeliefPtr()->clear();
			map<StatePtr, LD, StatePtrCmp> mp = vn->bNew->getMp();
			foreach(mpIter, mp) {
				StatePtr sptr = mpIter->first;
				LD pval = mpIter->second;
				psum += pval; 
		//		vn->getBeliefPtr()->setPr(sptr, pval);
			}
			if(psum > 0.0) {	
				foreach(mpIter, mp) {
				StatePtr sptr = mpIter->first;
				LD pval = mpIter->second;
				vn->getBeliefPtr()->setPr(sptr, pval/psum);
				}
			} else {
				cout<<"error in "<<*vn<<endl;	
				if(vn->getNodeNumVar() == 1) { // if
					vector<INT> L = VideoNode::getIntFromId(vn->getNodeId());
					if(evidence[(L[0]) ][(L[1])] == 1) { // moving object
						cout<<" is changing\n"; 
					}
				} else { // region approaches
				
				
				}
				getchar(); 
			}
		}
		return currTime; 
	}
void MVI::init() {
		vr = new VR(256, false); 
		vr->open_video(fname); 
		N1 = vr->getN1();
		N2 = vr->getN2(); 
		model = new ModelMVI(N1, N2, Q_potts);
		graph = new GraphMVI(TR, TT, model);
		vr->increment_frame(); 
	}
		
	
int MVI::run(bool ppm, int numLow) {
		vr->increment_frame();
	if(vr->frame_num > numLow) {
		vector<vector<INT> > fd = vr->frame_diff(vr->current_frame, vr->prev_frame, vr->diff_frame); 
		//foreach(iter, fd) { cout<<getStr<INT>(*iter)<<"\n"; }
		//BwImage img(vr->diff_frame);
		 
		//for(int i=0; i < N1; i++) {
		//	for(int j=0; j < N2; j++) {
		//		cout<<((INT)img[i][j])<<" "; 
		//	}
		//	cout<<"\n"; 
		//}
		// getchar(); 
			vector<vector<INT> > res; 
			if(ppm) {
				graph->update_evidence(fd); // update the evidence
				// graph->PPM_step2(0.1, 0.1); // update all cond Ps
				cout<<"PPM_step3: \n"; graph->PPM_step3(); 
				for(int i=0; i < 2; i++) { // 3 iterations of PPM_innerloop(); 
					cout<<" LOOP_"<<i<<endl; 
					graph->PPM_inner_loop();
			//	graph->updateNodesInfo(); 
			//	string stats = graph->get_status();
      		//	LD HR = graph->getHR();
		    //	LD lnZR = HR-graph->SR; // SR is a class variable
      		//	fprintf(stderr, "ITER  -lnZ_approx HR [SR  minM maxM <ParDiff> PD_max <PD_max> <Cblf> blfC_max <blfC_max>]\n");
      		//	fprintf(stderr, "%g: Z: %g H: %g [%s]\n", ctime, lnZR, HR, stats.c_str());
				}
				graph->updateTime(1.0);
				res =  graph->map_estimate();
				vr->get_frame(res, vr->res_frame, Q_potts);
				
			} else {
				res = fd; 
				for(int i=0; i < N1; i++) {
					for(int j=0; j < N2; j++) {
						if(fd[i][j] == (INT) 0) {
							res[i][j] = 0;
						} else {
							res[i][j] = 1; 
						}
					}
				}
				vr->get_frame(res, vr->res_frame, 2);
			}
		//	vr->get_frame(fd, vr->diff_frame, 2);
			
			int x = 0;
			//if( (vr->frame_num == 22) || (vr->frame_num == 100) || (vr->frame_num == 182) ) {
			vr->write_op();
			return vr->display(); 
			//} 
			// return x; 
		} else {
			return 0;
		}
}
