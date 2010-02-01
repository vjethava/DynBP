#ifndef POTTSGRAPH_H_
#define POTTSGRAPH_H_

#include "main.h" 
#include <cassert>
#include <gsl/gsl_sf_exp.h>
#include "GraphPPM.h" 
using namespace std; 
#include "VideoNode.h"
#include "PottsModel.h" 

  
class PottsGraphPPM : public GraphPPM {
	
public: 
    PottsModel* model;
    int rx, ry, cx, cy; 
    FactorGraph* fg; 
    LD currTime; 
    LD currDt;
// public:     
    PottsGraphPPM(PottsModel* _M, int _rx, int _ry, int _cx, int _cy) {
        model = _M; rx = _rx; ry = _ry; cx = _cx; cy = _cy;
        fg = new FactorGraph();  
    }
    void genSubRegions(); 
    VideoNode* addNode(VideoNode*); 
    vector<INT> getIthStVec(int idx, int N); 
    void initOrigP(VideoNode* vn, bool rndm=true); 
    //void update_pHat(VideoNode* vn, LD theta=1.0, LD dt=0.5);
    /// Initializes the single belief and conditional belief values; 
    void PPM_step2(LD theta, LD delta_t); 
    LD getStatePr(VideoNode* vn, StatePtr sp);
    LD get_p_hat(VideoNode* vn, StatePtr si, StatePtr so, LD theta, LD dt);
    void initCondP(VideoNode* vn, LD theta, LD dt); 
    StatePtr getChildStatePtr(NodePtr vp, NodePtr vc, StatePtr parState); 
    void PPM_algo(LD theta, LD dt);
    void PPM_algo2(LD theta, LD dt);
    LD getHR(); 
    void initFG();
    LD getVarMargPr(VideoNode*, int x, int y, int val);
    LD updateTime(LD dt); 
    bool stoppingCriterion(int iter, LD PD_thres=0.05, LD blf_thres=0.05, int iter_thres=4); 
    vector<vector<INT> > greedyDecision(); 
    void singleIter(LD theta_ppm, LD dt_ppm); 
    void initAllNodesOrigP(bool rndm);
   bool isValid(int x, int y); 
   void genBETHE(); 
   /** 
    * Returns the normalized beliefs as a product of region beliefs with appropriate normalization
    */
    vector<LD> getNormGraphBlfs(); 
    LD getGraphStBlf(int i);
    /// generates the ampl model file
    void writeModelFile(FILE* modelFile);
    /// generates the normalization constraints     
    void writeNormConstraints(FILE* modelFile);
    /// generates the ampl files for simulating a single ppm iteration 
    void getAMPL(string modelFileName, string dataFileName, string amplFileName); 
}; 

#endif /*POTTSGRAPH_H_*/

/* REVISION LOG
 * ============
 * July 16 2008. Support addition for AMPL modeling language
 * 
 * 
 */
 