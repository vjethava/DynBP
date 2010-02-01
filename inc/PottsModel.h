
#ifndef POTTSMODEL_H_
#define POTTSMODEL_H_
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_exp.h>
#include "belief.h"
struct PottsModel {
    /// the number of possible values each variable can take 1...Q
    int Q;
    /// variance for same node factors
    LD hSigma;
    /// variance for nearest node interactions
    LD jSigma;
    /// size of the grid
    int N1, N2;
    /// the see used for factor generation
    int seed;

    /// map containing singular interactions
    map<string, BeliefPtr> fieldMp;
    /// the random number generator
    gsl_rng* myRNG;

    PottsModel(int _N1, int _N2, int _Q=2, LD _jSigma=0.1, LD _hSigma=0.1, int _seed = 0) {
        N1 = _N1;
        N2 = _N2;
        jSigma = _jSigma;
        hSigma = _hSigma;
        seed = _seed;
        Q = _Q;
        const gsl_rng_type* T;
        T = gsl_rng_default;
        myRNG = gsl_rng_alloc(T);
        gsl_rng_set(myRNG, seed);
    }
    void genNodes();
    void updateBlfPtr(BeliefPtr bNew, int varNum);
    int GV(int k) {
        if(Q%2 == 0) {
            return 2*k-(Q-1);
        } else
            return (k-Q/2);
    }
    bool isValid(int r, int c) {
        return ( (r >=0) && (r < N1) && (c>= 0) && (c < N2) );
    }
    string getStats();
    void genFactors();

    LD getP(string fid, StatePtr sC) {
        ASSERT(fieldMp.find(fid) != fieldMp.end());
        ASSERT(fieldMp.find(fid)->second->getPr(sC) != BELIEF_NOT_FOUND);
        return fieldMp.find(fid)->second->getPr(sC);
    }
    LD getGlobalP(int sVal);
    LD getCondP(int sInit, int sNext, LD TD);
    vector<INT> getIthStVec(int i, int N);
};

/**
 * This class additionally provides the model init b(x_0) 
 * which is assumed as factored into singletons - (note)
 */
struct PottsModelInitP: public PottsModel {
    int init_p_seed;
    gsl_rng* pRNG;
    PottsModelInitP(int _N1, int _N2, int _Q=2, LD _jSigma=0.1, LD _hSigma=0.1, int _seed = 0):
    PottsModel(_N1, _N2, _Q, _jSigma, _hSigma, _seed) {
        init_p_seed = 0;
        const gsl_rng_type* T;
        T = gsl_rng_default;
        myRNG = gsl_rng_alloc(T);
        gsl_rng_set(myRNG, seed);

    }

    void setInitPseed(int _s) {
        init_p_seed = _s;
    }
	
	
};
/*
struct PottsModelBP: public PottsModel {
	
    vector<EdgePtr> M_bp; 
    vector<EdgePtr> N_bp;
    map<NodePtr, BeliefPtr, NodeCmp()> b_bp; 
    map<NodePtr, LD> Z;
    void init_bp();    
};
*/
#endif /*POTTSMODEL_H_*/
