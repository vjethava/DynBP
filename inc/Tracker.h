#ifndef TRACKER_H_
#define TRACKER_H_
#include "PottsGraph.h"
#include "VideoReader.h" 
using namespace std; 

class Tracker {
    int N2, N1; 
    int counter; 
    PottsModel* model; 
    PottsGraphPPM* graph;
    VideoReader* vreader;
    int C; // number of quantizations 
    LD epsilon; 
    LD theta; // the weight to similar states 
    IplImage* last_frame, *curr_frame; 
    IplImage* track_frame; 
public:
    Tracker(); 
    void initMain(string videoName, int NL, bool isGray, int C, LD theta, int rx, int ry, int cx, int cy); 
    void initGraph(int rx, int ry, int cx, int cy);
    void initVideo(string videoName, int NL, bool isGray);  
    void run(int priorNum , int maxNum , int frameJump);
    IplImage* getFrameDiff(IplImage* frame1, IplImage* frame2);
    void updateModel(IplImage* prev_frame, IplImage* next_frame);
    bool updateFrame(); 
    IplImage* greedyFillFrame();
    void run_step(LD theta, LD dt); 
};
#endif /*TRACKER_H_*/
