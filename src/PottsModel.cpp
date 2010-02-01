#include "PottsGraph.h"
#include "PottsModel.h"
#include "VideoNode.h"
#include <gsl/gsl_sf_exp.h>
using namespace std;

string PottsModel::getStats() {
	stringstream ss("");
	ss<<"{ size:"<<N1<<"x"<<N2<<" j: "<<jSigma<<" h: "<<hSigma<<" Q: "<<Q
			<<" seed: "<<seed<<"} ";
	return ss.str();
}

void PottsModel::genNodes() {
	for (int i=0; i < N1; i++) {
		for (int j=0; j < N2; j++) {
			if (isValid(i+1, j)) {
				stringstream ss("");
				ss<<i<<"."<<j<<"_"<<(i+1)<<"."<<j;
				string nid = ss.str();
				//          cout<<" node: "<<nid<<"\n";
				BeliefPtr bNew = new Belief();
				fieldMp.insert(make_pair<string, BeliefPtr>(nid, bNew));
			}

			if (isValid(i, j+1)) {
				stringstream ss("");
				ss<<i<<"."<<j<<"_"<<(i)<<"."<<(j+1);
				string nid = ss.str();
				//   cout<<" node: "<<nid<<"\n";
				BeliefPtr bNew = new Belief();
				fieldMp.insert(make_pair<string, BeliefPtr>(nid, bNew));
			}

			stringstream ss3("");
			ss3<<i<<"."<<j;
			string nid = ss3.str();
			BeliefPtr bNew = new Belief();
			fieldMp.insert(make_pair<string, BeliefPtr>(nid, bNew));
		}
	}
}

void PottsModel::updateBlfPtr(BeliefPtr bNew, int varNum) {
	if (varNum == 1) {
		Belief bTmp;
		LD sumP =0.0;
		for (int k=0; k < Q; k++) {
			LD hCurr = (LD) gsl_sf_exp( -gsl_ran_gaussian(myRNG, (double) hSigma)*GV(k));
			INT data[1];
			data[0] = ((INT)k);
			StatePtr sC = new State(1, data);
			bTmp.setPr(sC, hCurr);
			sumP += hCurr;
		}
		map<StatePtr, LD, StatePtrCmp> mp = bTmp.getMp();
		foreach(iter, mp) {
			bNew->setPr(iter->first, iter->second/sumP);
		}
	} else {
		Belief bTmp;
		LD sumP=0.0;
		for (int k=0; k < Q; k++) {
			for (int l=0; l < Q; l++) {
				LD jCurr = (LD)gsl_sf_exp(-gsl_ran_gaussian(myRNG,
						(double)jSigma)*GV(k)*GV(l));
				INT data[2];
				data[0] = ((INT)k);
				data[1] = ((INT)l);
				StatePtr sC = new State(2, data);
				bTmp.setPr(sC, jCurr);
				sumP += jCurr;

			}
		}
		map<StatePtr, LD, StatePtrCmp> mp = bTmp.getMp();
		foreach(iter, mp) {
			bNew->setPr(iter->first, iter->second/sumP);
		}

	}
}
void PottsModel::genFactors() {
	genNodes();
	foreach(iter, fieldMp) {
		string node_name = iter->first;
		BeliefPtr bNew = iter->second;
		vector<string> vs = split(node_name, '_');
		int varNum = vs.size();
		//  cout<<"node: "<<node_name<<" varNum: "<<varNum<<endl; getchar();
		updateBlfPtr(bNew, varNum);
	}
}

/**
 * Returns the proportionality constant - has to be divided by Z(t) to get 
 * true probability 
 *  
 */
LD PottsModel::getGlobalP(int sVal) {
	LD p = 1.0;
	int numStates = (int) MyDef::pow(Q, N1*N2);
	vector<INT> cSt = getIthStVec(sVal, N1*N2);
	// i - N1, j - N2
	foreach(iter, fieldMp) {
		string nId = iter->first;
		BeliefPtr nBlf = iter->second;
		vector<string> nodes = split(nId, '_');
		if(nodes.size() == 2) {
			vector<INT> v;
			foreach(niter, nodes) {
				int i, j;
				sscanf(niter->c_str(), "%d.%d", &i, &j);
				int idx = i*N2 + j;
				int val = cSt[idx];
				v.push_back(val);
			}
			StatePtr sptr = new State(v);
			LD pCurr = nBlf->getPr(sptr);
			ASSERT(pCurr != BELIEF_NOT_FOUND);
			p *= pCurr;
//			printf("getGlobalP(%d) nId: %s pCurr: %g p: %g\n", sVal, nId.c_str(), pCurr, p);
//			getchar(); 
			delete(sptr);
		} else {
			int i, j;
			sscanf(nId.c_str(), "%d.%d", &i, &j);
			int idx = i*N2 + j;
			int val = cSt[idx];
			vector<INT> v;
			v.push_back(val);
			StatePtr sptr = new State(v);
			LD pCurr = nBlf->getPr(sptr);
			ASSERT(pCurr != BELIEF_NOT_FOUND);
			p *= pCurr;
		//	printf("getGlobalP(%d) nId: %s pCurr: %g, p: %g\n", sVal, nId.c_str(), pCurr, p);
		//	getchar();
			delete(sptr);
		}
	}
	// printf("getGlobalP(%d) p: %g\n", sVal, p); 
	//getchar(); 
	return p;
}
/**
 * Returns the P_cond using Kikuchi model 
 */
LD PottsModel::getCondP(int vInit, int vNext, LD TD) {
	LD p = 1.0;
	vector<INT> sInit = getIthStVec(vInit, N1*N2);
	vector<INT>	sNext = getIthStVec(vNext, N1*N2);
	int diff = 0; 
	for(int i=0; i < N1*N2; i++) if(sInit[i] != sNext[i]) diff++; 
	LD pInit = getGlobalP(vInit);
	LD pNext = getGlobalP(vNext);
	LD Idiff = pow(TD, diff)*pow((1.0-TD), (N1*N2-diff));
	p = pNext/pInit*Idiff; 
/*	cout<<" diff: "<<diff<<" getCondP("<<vInit<<", "<<vNext<<") "; 
	cout<<" pOrig: "<<pInit<<" pNext: "<<pNext<<" Idiff: "<<Idiff<<" p: "<<p<<"\n"; 
	getchar(); */  
	return p;
}

/**
 * Returns the ith state given grid size of N
 */
vector<INT> PottsModel::getIthStVec(int idx, int N) {
	vector<INT> cval;
	for (int j=0; j < N; j++) {
		int id = idx%Q;
		int vl = id;
		cval.push_back(vl);
		idx = idx/Q;
	}
	return cval;
}
/*
void PottsModelBP::init_bp() {
	foreach(iter, adj_map) {
		VideoNode* vn = new VideoNode((string) iter->first); 
	}
	foreach(iter, adj_map) {
		
	}
	
}
*/ 
