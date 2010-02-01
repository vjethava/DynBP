#include "Tracker.h"

Tracker::Tracker() {
    vreader = NULL;
    graph = NULL;
    model = NULL; 
    last_frame = NULL; 
    track_frame = NULL; 
    curr_frame = NULL; 
    C = 0; 
    epsilon = 0; 
    counter = 0; 
 
}

void Tracker::initVideo(string videoName, int NL, bool isGray) {
    FPRINTF(stderr, "Tracker::initVideo(NL: %d gray: %d)\n", NL, isGray); 
    vreader = new VideoReader(isGray, NL); 
    vreader->openVideo(videoName); 
    // get the last frame
    allocateOnDemand(&track_frame, vreader->frame_size, IPL_DEPTH_8U, 1); 
    allocateOnDemand(&(vreader->input_frame), vreader->frame_size, IPL_DEPTH_8U, 1); 
    allocateOnDemand(&(vreader->last_frame), vreader->frame_size, IPL_DEPTH_8U, 1); 
    last_frame = vreader->last_frame;
    curr_frame = vreader->input_frame; 
}

void Tracker::initGraph(int rx, int ry, int cx, int cy) {
    Timer timer; 
    FPRINTF(stderr, "Tracker::initGraph(rx: %d ry: %d cx: %d cy: %d)\n", rx, ry, cx, ry); 
    ASSERT(vreader != 0);
    LD pEq = 1.0/((LD) C);  
    N1 = vreader->getN1();
    N2 = vreader->getN2();  
    model = new PottsModel(N1, N2, C); 
    model->genFactors();
    // equalize the map here 
    // cout<< " mode->N1: "<<model->N1<<"  model->N2: "<<model->N2<<" model->Q: "<<model->Q<<endl; getchar(); 
    foreach(factorIter, model->fieldMp) {
        string node_name = factorIter->first;
        vector<string> name_parts = split(node_name, '_'); 
        bool isSingleton = true; 
        if( name_parts.size() == 2) { //"0.0_0.1" 
            isSingleton = false;
        }
      //  cout<<" node: "<<node_name<<" isSingleton: "<<isSingleton<<"\n" ; getchar(); 
        BeliefPtr bptr = factorIter->second;
        map<StatePtr, LD, StatePtrCmp> mp = bptr->getMp();
        bptr->clear();  
        foreach(stIter, mp) {
            StatePtr cSt = stIter->first;
            timer.tick(); 
            
            if(isSingleton) {
                bptr->setPr(cSt, pEq);
           
            } else {
                ASSERT(cSt->num == 2); 
                INT v1 = cSt->getData(0); 
                INT v2 = cSt->getData(1);
                timer.tick(); 
                if(v1 == v2) {
                    bptr->setPr(cSt, theta);
                     
                } else {
                    bptr->setPr(cSt, epsilon);
                }
            }
      //      LD pSt = model->getP(node_name, cSt); 
      //      cout<<" state: "<<*cSt<<" pSt: "<< pSt<<"\n";    
        }
    }
    graph = new PottsGraphPPM(model, rx, ry, cx, cy);
    graph->genBETHE(); 
    graph->initAllNodesOrigP(false); 
    FPRINTF(stderr, " done\n"); 
    
}

void Tracker::updateModel(IplImage* prev_frame, IplImage* next_frame) {
    Timer timer; 
    FPRINTF(stderr, "Tracker::updateModel() "); 
    BwImage currImg(prev_frame);
    BwImage nextImg(next_frame); 
    
    for(int i=0; i < N1; i++) {
        for(int j=0; j < N2; j++) {
            LD pEq = 1.0/((LD) C);  
            stringstream ssn(""); ssn<<i<<"."<<j;
            string nid = ssn.str(); 
            BeliefPtr bptr = model->fieldMp.find(ssn.str())->second;
            map<StatePtr, LD, StatePtrCmp> mp = bptr->getMp();
            bptr->clear(); 
            timer.tick(); 
            if(currImg[i][j] == nextImg[i][j]) {
                foreach(mpIter, mp) {
                    bptr->setPr(mpIter->first, pEq);                       
                }
            } else {
                foreach(mpIter, mp) {
                    if(mpIter->first->getData(0) == (C-1)) 
                        bptr->setPr(mpIter->first, 1.0);
                    else 
                        bptr->setPr(mpIter->first, 0.0);  
                }
            }
        }
    }
    FPRINTF(stderr, " done\n"); 
}

void Tracker::initMain(string videoName, int NL, bool isGray, int _C, LD _theta, int rx, int ry, int cx, int cy) {
    FPRINTF(stderr, "Tracker::initMain()\n"); 
    initVideo(videoName, NL, isGray); 
    theta = _theta;
    C = _C;
    epsilon =   ((LD)(1.0 - theta))/((LD) (C-1));
    initGraph(rx, ry, cx, cy); 
    
}

bool Tracker::updateFrame() { // goto the next frame
     FPRINTF(stderr, "Tracker::updateFrame()\n"); 
    counter++;
    bool res = vreader->incrementFrame(last_frame, curr_frame);
    cvNamedWindow("Tracker.updateFrame().last_frame", CV_WINDOW_AUTOSIZE);
    cvNamedWindow("Tracker.updateFrame().curr_frame", CV_WINDOW_AUTOSIZE);
    cvShowImage("Tracker.updateFrame().last_frame", last_frame);
    cvShowImage("Tracker.updateFrame().curr_frame", curr_frame);
    IplImage* diff_frame = NULL;
    int NP=0; 
    double diff = vreader->getFrameDifference(curr_frame, last_frame, diff_frame, &NP);  
    ASSERT(curr_frame != NULL);
    ASSERT(last_frame != NULL);
    ASSERT(curr_frame != last_frame); 
    cvNamedWindow("Tracker.updateFrame().diff_frame", CV_WINDOW_AUTOSIZE);
    cvShowImage("Tracker.updateFrame().diff_frame", diff_frame);
    cout<<" Tracker.updateFrame() diff = "<<diff<<endl; //getchar(); 
    return res;  
}

void Tracker::run_step(LD theta_ppm, LD dt_ppm) {
    Timer timer; 
     FPRINTF(stderr, "Tracker::run_step()\n"); 
      int local_count = 0;    
INIT_STEP:
    local_count++;  
    bool res = updateFrame();
    if(res == false) { cout<<" reached end at frame "<<counter<<"!\n"; exit(-1); } 
    if(local_count > 1) {
        
        ASSERT((last_frame != NULL) && (curr_frame != NULL) && (last_frame != curr_frame));
        updateModel(last_frame, curr_frame); 
        graph->singleIter(theta_ppm, dt_ppm); 
        vector<vector<INT> > decision = graph->greedyDecision();
        BwImage image(track_frame); 
        int N1 = vreader->getN1();
        int N2 = vreader->getN2() ; 
        fi(0, N1) fj(0, N2) image[i][j] = decision[i][j];
        stringstream ss(""); ss<<"track_"<<local_count<<".png";   
        cvSaveImage(ss.str().c_str(), track_frame);
        cvNamedWindow("TRACK", CV_WINDOW_AUTOSIZE); 
        cvShowImage("TRACK", track_frame); 
     //   cvWaitKey(1000);
        getchar(); 
    }
    goto INIT_STEP; 
   
}

/*
IplImage* Tracker::greedyFillFrame() {
    BwImage image(track_frame);
    int N1 = vreader->getN1();
    int N2 = vreader->getN1();
    
    bool unfilled[N1][N2][C];
    int stcount[N1][N2][C];
    LD stsum[N1][N2][C]; 
    fi(0, N1) { 
        fj(0, N2) { 
            image[i][j] = 0; 
            unfilled[i][j] = true;
            fk(0, C) {
                stcount[i][j][k] = 0;
                stsum[i][j][k] = 0.0;
            }      
        }
    }
    
   
    
//        if(vn->getNodeNumVar() == ((graph->rx)*(graph*ry))) {
//            vector<int> L = VideoNode::getIntFromId(vn->getNodeId()); 
//            for(int i=0; i < L[2]; i++) {
//                for(int j=0; j < L[3]; j++) {
//                    int x = L[0]+i;
//                    int y = L[1]+j;
//                    int idx = i*L[3] + j; 
//                    for(int k=0; k < C; k++) {
//                        
//                    
//                    }
//                }
//            }
//        }
    }
    return track_frame; 
}
*/ 

