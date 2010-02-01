#ifndef VIDEOGRAPH_H_
#define VIDEOGRAPH_H_
#include <stdio.h>
#include <opencv/cv.h>
#include <math.h>
#include <opencv/highgui.h> 
#include <string>
#include <iostream>
#include "Image.h"
#include "state.h"
#include "belief.h" 
#include "mygraph.h" 
#include "main.h" 
#include "VideoNode.h"
#include "GraphPPM.h" 
using namespace std; 
class VideoReader; 
typedef pair<StatePtr, LD> SPLD; 

class VideoGraph : public GraphPPM
{
    /// the belief maps for number of variables under question
    map<int, ConditionalBelief* > blfMp;
    /// the belief maps for bOrig
    //map<pair<int, int>, BeliefPtr, PRCMP<int, int> > pOrigMp;  
    map<int, BeliefPtr> pOrigMp;  
    /// the levels
    int numLvls;
    /// whether is gray or not  
    bool isGray; 
    /// the levels that will be present in the video (0..255) 
    vector<int> levels; 
    /// the reader 
    VideoReader* reader;
    /// biggest region width
    int rx;
    /// biggest region size
    int ry;
    /// how much x-overlap between current and next region
    int cx;
    /// how much y-overlap between current and next region
    int cy;
    /// how many rows
    int N1;
    /// how many columns
    int N2;
    /// the number of cond beliefs that have been duplicated
    int hits; 
    /// the number of frames that have been  missed
    int misses;
     
    /// the sigma  value for gaussian around the mean value 
    LD sigma, nearestT; 
    /// the threshold for max-probability 
    LD selectThreshold; 
    int parAugmentNotFound, parAugmentFound;
    pair<StatePtr, LD> getBestState0(VideoNode* vn);  
public:
    /// the inpainting frame
    IplImage* inpainted_frame; 
    IplImage* local_mask; 
    /// pointer to reader output frame
    IplImage* output_frame; 
    /// pointer to the reader true frame
    IplImage* true_frame;
    BeliefPtr getOrigBlfFromMp(int key) {
        if(pOrigMp.find(key) != pOrigMp.end()) return pOrigMp.find(key)->second;
        BeliefPtr bptr = new Belief(true); 
        pOrigMp.insert(make_pair<int, Belief* >(key, bptr)); 
        return bptr;  
    }
    VideoNode* addNode(VideoNode* np) {
        VideoNode* nn = (VideoNode*) ((MyGraph*) this)->addNode(np);
        setNodeCondBlf(nn); 
        ASSERT(nn->getCPTR() != 0); 
        FILE* nodeFile=fopen("node_list", "a"); 
        fprintf(nodeFile, "%s\n", np->getNodeId().c_str());
        fclose(nodeFile);  
        return nn; 
    }
    
   
    map<VI, LD, VCMP<INT> > getGaussianNeighbours(VI v, LD NH, LD T); 
	VideoGraph(VideoReader* _reader, bool _isGray=true, int numL=4, int _rx=4, int _ry=4, int _cx=2, int _cy=2);
	virtual ~VideoGraph();
    void updateCondBlfMp(IplImage* last_frame, IplImage* curr_frame, IplImage* mask_frame);
    int addCondBlf(vector<INT> pastVals, vector<INT> nextVals);
    void setNodeCondBlf(VideoNode* vn);  
    void genSubRegions(); // generate the smaller regions while still possible
    vector<VideoNode*> getMaximalUnfilledRegions(); 
    StatePtr readNodeFromFrame(VideoNode* vn, IplImage* frame);
    //void getNextStateCands(VideoNode* vn, IplImage* frame, LD T=0.1, vector<StatePtr>* candidates, vector<LD>* probs); 
    vector<pair<StatePtr, LD> > getNextStateCands(VideoNode* vn, IplImage* frame, LD T=0.1); 
    int updateNodeBlf(VideoNode* vn, IplImage* frame); 
    int updateAllNodes(IplImage* frame, IplImage* mask_frame); 
    void getNextFramePPM(IplImage* frame, IplImage* holes); 
    LD getNodeMaskPct(VideoNode* vn, IplImage* mask_frame); 
    void initMessages(VideoNode* vnp, VideoNode* vnc, map<VI, LD, VCMP<INT> > choices);
    void updateMessages(VideoNode* vnp, VideoNode* vnc, map<VI, LD, VCMP<INT> > choices);
    void PPM_step2(IplImage* frame, LD candT); 
    bool isValid(int x, int y, int xi=0, int yi=0, int xl=-1, int yl=-1);
    bool getAbsolutePixel(string nid, int i, int j, StatePtr sptr, vector<INT>* result);
    StatePtr getChildStatePtr(string parId, string childId, StatePtr parState);
    void updateChildAugments(VideoNode* vn, StatePtr sOld);
    void updateParentAugments(VideoNode* vn, StatePtr oldState); 
    IplImage* inPaint(IplImage* input_frame, IplImage* mask_frame);
    vector<pair<int, int> > getNodeUnfilledPixels(VideoNode* vn, IplImage* mask_frame); 
    bool isMaskedPixel(int i, int j, IplImage* mask_frame); 
    vector<pair<StatePtr, LD> > getMatchingStates(VideoNode* vn, IplImage* frame, IplImage* mask_frame, bool approx=false);
    void updateMaskStatus(IplImage* mask_frame);  
    void fillFrame0();
    void fillFrame1(bool lookMax=false); 
    bool isPartialMatch(StatePtr s1, StatePtr s2, StatePtr mask);
    BeliefPtr getApproxBlf(VideoNode* vn, StatePtr sOld);
    vector<pair<VideoNode*, LD> > getIncompleteList(IplImage* mask_frame); 
    bool isClose(StatePtr s1, StatePtr s2, int NL, LD distS);
    vector<pair<StatePtr, LD> > getMaximalPixelPr(VideoNode* vn, int i, int j); 
    void clearAllNodes();
    void clearCondBlfMp();
    void initLocals(); 
    void getCandidates(IplImage* last_frame,IplImage* input_frame,  IplImage* mask_frame);
    StatePtr getChildStatePtr(NodePtr vp, NodePtr vc, StatePtr parState) {
        return getChildStatePtr(vp->getNodeId(), vc->getNodeId(), parState); 
    } 
};


#endif /*VIDEOGRAPH_H_*/
