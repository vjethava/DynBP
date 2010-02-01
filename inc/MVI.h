#ifndef MVI_H_
#define MVI_H_
#include "PottsGraph.h"
#include "VR.h"
#include "Image.h" 
#include <string>
#include <iostream>
using namespace std; 

class ModelMVI: public PottsModel {
public: 
	ModelMVI(int N1, int N2, int Q): PottsModel(N1, N2, Q) { }
	LD getP(string fid, StatePtr sC) {
		vector<int> idx = VideoNode::getIntFromId(fid); 
		int numvar = idx[2]*idx[3]; 
		int numS = (int) MyDef::pow(Q, numvar); 
		LD L = 1.0/((LD) numS); 
		return L; 
	}
}; 

class GraphMVI : public PottsGraphPPM {
	LD theta_r;
	LD theta_t;
	LD epsilon_r; 
	LD epsilon_t; 
	// the constant to avoid algo trouble 
	LD epsilon_zero;
	vector<vector<INT> > evidence; 
	int C; 
	int N1, N2; 
	StatePtr ST[2]; 
	StatePtr SC; 
public:
	LD get_p_t(int i, int j, INT xi, INT xo);
	 
	LD getStatePr(VideoNode* vn, StatePtr sp) { 
		int numS = (int) MyDef::pow(model->Q, vn->getNodeNumVar()); 
		LD L = 1.0/((LD) numS); 
		return L; 
	}
	
	GraphMVI(LD S, LD T, ModelMVI* _M, int _rx=2, int _ry=2, int _cx=1, int _cy=1) :
		PottsGraphPPM(_M, _rx, _ry, _cx, _cy)
	{
		epsilon_zero = 0.0001; 
		N1 = model->N1; 
		N2 = model->N2; 
		C = model->Q; 
		theta_r = S; 
		theta_t = T; 
		epsilon_r = (1.0-theta_r)/((LD) (C-1));
		epsilon_t = (1.0-theta_t)/((LD) (C-1));
		genBETHE();
		cout<<" initializing origP for all nodes\n";
		initAllNodesOrigP(false); cout<<"done\n";
		int v0 = 0, v1 = 1, vc = (C-1); 
		ST[0] = new State(1, &v0); 
		ST[1] = new State(1, &v1); 
		SC = new State(1, &vc); 
		 
	}
	
	inline void update_evidence(vector<vector<INT> > _evidence) {
		Timer timer; 
		cout<<"update_evidence() & initAllNodesCondP()\n"; 
		evidence = _evidence;
		foreach(iter, adj_map) {
			VideoNode* vn = (VideoNode*) (iter->first); 
			initCondP(vn, 0.1, 0.1);
			timer.tick(); 
		}
		cout<<"done\n"; 
		
	}
	vector<vector<INT> > map_estimate(); 
	vector<vector<INT> > mmse_estimate();
	LD get_p_hat(VideoNode* vn, StatePtr si, StatePtr so, LD theta, LD dt); 
	//LD get_p_t(VideoNode* vn, INT si, INT so) {}
	void initCondP(VideoNode* vn, LD theta, LD dt);  
	// void initOrigP(VideoNode* vn, bool rndm=true); 
	LD updateTime(LD dt); 
	void initAllNodesOrigP(bool rndm=false);
};

class MVI {
	int N1, N2; 
	int Q_potts;
	LD TR, TT;  
	string fname;
	VR* vr;
	GraphMVI* graph;
	ModelMVI* model;
public:  
	MVI(string _fname,  int _Q_potts = 4, LD _TR=0.7, LD _TT=0.7) { 
		fname = _fname; 
		Q_potts = _Q_potts; 
		TR = _TR; 
		TT = _TT; 
	}  
	void init(); 
	int run(bool ppm=true, int numLow=0); 
	
}; 
#endif /*MVI_H_*/
