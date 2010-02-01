#include "VideoReader.h"
#include "VideoGraph.h"
#include <vector>

using namespace std;
void allocateOnDemand( IplImage **img, CvSize size, int depth, int channels ) {
    if ( *img != NULL )
        return;

    *img = cvCreateImage( size, depth, channels );
    if ( *img == NULL ) {
        fprintf(stderr, "Error: Couldn't allocate image.  Out of memory?\n");
        exit(-1);
    }
}

VideoReader::VideoReader(bool _gray, int _numLvls) {
    frame_count=0;
    count11=0;
    input_video = NULL;
    input_frame = NULL;
    numLvls=2;
    input_frame= NULL;
    output_frame = NULL;
    mask_frame= NULL;
    input_window_handle = NULL;
    output_window_handle = NULL;
    mask_window_handle = NULL;
    frame = NULL;
    frame1 = NULL;
    frame1_1C = NULL;
    frame2_GS = NULL;
    grayscale = _gray;
    numLvls = _numLvls;
    last_frame = NULL;
    graph = NULL;
    diff_frame= NULL;
    output_video = NULL;
    difference_video = NULL;

    //output_video = cvCreateVideoWriter("output.avi", CV_FOURCC('D','I','B',' '), fps, frame_size, is_color);



}

VideoReader::~VideoReader() {
    if(input_video != NULL) {
        cvReleaseCapture(&input_video);
    }

    cvReleaseVideoWriter(&output_video);
    cvReleaseVideoWriter(&difference_video);


}

void VideoReader::openVideo(string videoFileName) {
#ifdef LAPTOP
    number_of_frames = 240;
    frames_per_second = 25;
    frame_size.height = 100;
    frame_size.width = 300;
#else

    input_video = cvCaptureFromFile( videoFileName.c_str() );
    if (input_video == NULL) {
        fprintf(stderr, "Error: Can't open video : %s\n", videoFileName.c_str());
        exit(-1);
    }
    IplImage* frame11 = cvQueryFrame( input_video );
    frame_size.height =
        (int) cvGetCaptureProperty( input_video, CV_CAP_PROP_FRAME_HEIGHT );
    frame_size.width =
        (int) cvGetCaptureProperty( input_video, CV_CAP_PROP_FRAME_WIDTH );

    cvSetCaptureProperty( input_video, CV_CAP_PROP_POS_AVI_RATIO, 1. );
    number_of_frames = (long) cvGetCaptureProperty( input_video, CV_CAP_PROP_POS_FRAMES );
    //cvSetCaptureProperty( input_video, CV_CAP_PROP_POS_FRAMES, 0. );

    cvSetCaptureProperty( input_video, CV_CAP_PROP_POS_FRAMES, 0. );

    frames_per_second = (int) cvGetCaptureProperty( input_video, CV_CAP_PROP_FPS);
    int depth = (int) frame11->depth;
    int nChannels = (int) frame11->nChannels; 
    cout<<"NUM: "<<number_of_frames<<"\n";
    printf(" video %s num_frames: %d  fps: %d fheight: %d fwidth: %d depth: %d channels: %d\n", videoFileName.c_str(), number_of_frames, frames_per_second, frame_size.height, frame_size.width, depth, nChannels );
    //fwidth = (int) cvGetCaptureProperty( input_video, CV_CAP_PROP_FRAME_HEIGHT );
    //fheight = (int) cvGetCaptureProperty( input_video, CV_CAP_PROP_FRAME_WIDTH );
#endif

    allocateOnDemand( &mask_frame, frame_size, IPL_DEPTH_8U, 1 );

    //  cvNamedWindow("input_video",  CV_WINDOW_AUTOSIZE);
    //  cvNamedWindow("last_frame", CV_WINDOW_AUTOSIZE);
    //  cvNamedWindow("ppm_pred", CV_WINDOW_AUTOSIZE);
    // cvNamedWindow("mask_frame", CV_WINDOW_AUTOSIZE);
    //   last_frame_window_handle = cvGetWindowHandle("last_frame");
    //   output_window_handle = cvGetWindowHandle("ppm_pred");
    //  input_window_handle = cvGetWindowHandle("input_video");
    // mask_window_handle = cvGetWindowHandle("mask_frame");
    BwImage img(mask_frame);
    /* for(int i=0; i < frame_size.height; i++) {
         for(int j=0; j < frame_size.width; j++) {
             img[i][j] = 255; 
         }
     }*/
    // cout<<" bw image pixel val: "<<((int) img[0][0])<<"\n"; getchar();
}


void VideoReader::display() {
    // display the input frame
    //const char* wname = cvGetWindowName(input_window_handle);
    const char* wname = "input_video"; 
    cvNamedWindow(wname, CV_WINDOW_AUTOSIZE); 
    cvShowImage(wname, frame);
    
    const char* wname2= "last_frame"; 
    cvNamedWindow(wname2, CV_WINDOW_AUTOSIZE); 
    //wname2 = cvGetWindowName(last_frame_window_handle);
    cvShowImage(wname2, last_frame);

    //  wname = cvGetWindowName(mask_window_handle);
    // cvShowImage(wname, mask_frame);
    //wname3 = cvGetWindowName(output_window_handle);
    //const char* wname3 = "output";
    //cvNamedWindow(wname3, CV_WINDOW_AUTOSIZE); 
    //cvShowImage(wname, output_frame);
 //   cvWaitKey(10);
}


bool VideoReader::convertInputFrame(bool gray, int numLvls) {

    if (frame == NULL) {
        fprintf(stderr, "Error: Hmm. The end came sooner than we thought.\n");
        return false;
    }

    allocateOnDemand( &frame1, frame_size, IPL_DEPTH_8U, 3 );
    cvConvertImage(frame, frame1, 0);
    if(gray) {
        allocateOnDemand( &frame1_1C, frame_size, IPL_DEPTH_8U, 1 );
        cvConvertImage(frame1, frame1_1C, 0);
        BwImage bwImg(frame1_1C);
        if(numLvls < 256) {
            int div = 256/numLvls;
            for(int i=0; i < frame_size.height; i++) {
                for(int j=0; j < frame_size.width; j++) {
                    unsigned char nval = (unsigned char) (((int) bwImg[i][j]/div)*div);
                    //    cout<<"orig: "<<((int)bwImg[i][j])<<" NL: "<<numLvls<<" val: "<<((int)nval)<<"\n";
                    bwImg[i][j] = nval;
                }
            }
        }
        input_frame = frame1_1C;
    } else {
        RgbImage rgbImg(frame1);
        if(numLvls < 256) {
            int div = 256/numLvls;
            for(int i=0; i < frame_size.height; i++) {
                for(int j=0; j < frame_size.width; j++) {
                    rgbImg[i][j].b = (unsigned char) (((int) rgbImg[i][j].b/div)*div);
                    rgbImg[i][j].g = (unsigned char) (((int) rgbImg[i][j].g/div)*div);
                    rgbImg[i][j].r = (unsigned char) (((int) rgbImg[i][j].r/div)*div);
                }
            }
        }
        input_frame = frame1;
    }
    if(graph != NULL) {
        graph->true_frame = input_frame;
    }
    return true;

}

vector<INT> VideoReader::getDataFromFrame(IplImage* frame, int x, int y, int r, int c) {
    vector<INT> vals;
    int nc = frame->nChannels;
    bool grayscale = false;
    if(nc == 1)
        grayscale = true;
    for(int i=x; i < (x+r); i++) {
        for(int j=y; j < (y+c); j++) {
            if(grayscale) {
                BwImage img(frame);
                if(isValidPixel(frame, i, j) ) {
                    vals.push_back((INT) img[i][j]);
                }
            } else {
                RgbImage img(frame);
                if(isValidPixel(frame, i, j)) {
                    vals.push_back((INT) img[i][j].r);
                    vals.push_back((INT) img[i][j].g);
                    vals.push_back((INT) img[i][j].b);
                }
            }
        }
    }
    //  cout<<"VideoReader::getDataFromFrame("<<x<<", "<<y<<", "<<r<<", "<<c<<" gray="<<grayscale<<") vals="<<vals.size()<<"\n";
    return vals;
}

inline bool VideoReader::isValidPixel(IplImage* frame, int x, int y) {
    int R = frame->height;
    int C = frame->width;
    return ( (x < R) && (y < C) && (x >=0) && (y >= 0) );
}

void VideoReader::updateMaskFrame(int x, int y, int l, int w) {
    BwImage image(mask_frame);
    for(int i=x; i < (x+l); i++) {
        for(int j=y; j < (y+w); j++) {
            image[i][j] = 255;
        }
    }
    // cvNamedWindow("mask_window", CV_WINDOW_AUTOSIZE);
    // cvShowImage("mask_window", mask_frame);
    // cvWaitKey(0);
}

void VideoReader::updateLastFrame(IplImage* frame) {
    IplImage* input_frame;
    if(frame == NULL) {
        input_frame = this->input_frame;
    } else {
        input_frame = frame;
    }
    ++frame_count;
    // cvReleaseImage(&last_frame);
    if(input_frame != NULL) {
        allocateOnDemand( &last_frame, frame_size, input_frame->depth, input_frame->nChannels);
        cvCopy(input_frame, last_frame, 0);
        allocateOnDemand( &diff_frame, frame_size, input_frame->depth, input_frame->nChannels);
        allocateOnDemand( &output_frame, frame_size, last_frame->depth, last_frame->nChannels);
        graph->output_frame = output_frame;
    }

    cout<<"updateLastFrame("<<frame_count<<")\n"<<flush;
}

double VideoReader::getFrameDifference(IplImage* frame1, IplImage* frame2, IplImage* resFrame, int *NP) {
    double res=0;
    int div = 256/numLvls;

    if(NP != NULL) {
        *NP = 0;
    }
    int N = 1;
    if(grayscale) {
        if(resFrame == NULL) {
            resFrame = cvCreateImage(cvSize(frame1->width, frame1->height), IPL_DEPTH_8U,1);
        }
        N = frame1->height*frame1->width;
        BwImage img1(frame1);
        BwImage img2(frame2);
        BwImage resImg(resFrame);
        fi(0, frame1->height) {
            fj(0, frame1->width) {
                int diff = (int) abs( ((int) (img1[i][j]))-((int) (img2[i][j])) );
                resImg[i][j] = (INT) diff;
                double p = ((double) diff)/((double) div);
                if(p > 0.0) {
                    *NP += 1;
                }
                res += p*p;
            }
        }
        res = sqrt(res/N);
    } else {
        if(resFrame == NULL) {
            resFrame = cvCreateImage(cvSize(frame1->width, frame1->height), frame1->depth, frame1->nChannels);
        }
        N = frame1->height*frame1->width;
        RgbImage img1(frame1);
        RgbImage img2(frame2);
        RgbImage resImg(resFrame);
        fi(0, frame1->height) {
            fj(0, frame1->width) {
                int diff[3];
                diff[0] = (int) abs( ((int) (img1[i][j].r))-((int) (img2[i][j].r)) );
                diff[1] = (int) abs( ((int) (img1[i][j].g))-((int) (img2[i][j].g)) );
                diff[2] = (int) abs( ((int) (img1[i][j].b))-((int) (img2[i][j].b)) );
                resImg[i][j].r = (INT) diff[0];
                resImg[i][j].g = (INT) diff[1];
                resImg[i][j].b = (INT) diff[2];
                double p = 0;
                fk(0, 3) p += ((LD) diff[i])/((LD) div);
                p = p/3.0;
                if(p > 0.0) {
                    *NP += 1;
                }
                res += p*p;
            }
        }
        res = sqrt(res/N);

    }
    cvNamedWindow("difference");
    cvShowImage("difference", resFrame);
    return res;
}



void VideoReader::runMain(bool gray, int numLvls, int maxNum, int priorNum, int frameJump) {
    int is_color = ((gray)?0:1);
    int fps = frames_per_second;
    output_video = cvCreateVideoWriter("output.avi", CV_FOURCC('I','4','2','0'), fps, frame_size, is_color);
    difference_video = cvCreateVideoWriter("difference.avi", CV_FOURCC('I','4','2','0'), fps, frame_size, is_color);
    FILE* rFile = fopen("run_main.log", "w");
    fprintf(rFile, "#############################################################\n");
    fprintf(rFile, "# gray: %d NL: %d priorNum: %d maxNum: %d\n", gray, numLvls, priorNum, maxNum);
    fprintf(rFile, "#############################################################\n");
    fclose(rFile);
    graph = new VideoGraph(this, gray, numLvls);
    bool nf=true;
    int key_pressed;
    int count=0;
    IplImage* local_frame=NULL;
    bool useOpFrame=false;
    double diff= 0.0;

    graph = new VideoGraph(this, gray, numLvls);
    IplImage* result_image = NULL;

    while(nf && (count < maxNum)) {
        frame = cvQueryFrame(input_video);
        nf = convertInputFrame(gray, numLvls);

        if(count == 0) {
            initPrevFrames(frameJump);
        } else {
            updatePrevFrames(local_frame, frameJump);
        }
        if(count < frameJump) {
            graph->updateAllNodes(input_frame, mask_frame);
            local_frame = input_frame;
        } else if(count < priorNum) {
            graph->updateAllNodes(input_frame, mask_frame);
            graph->updateCondBlfMp(input_frame, last_frame, mask_frame);
            local_frame = input_frame;
            display();
        } else {
            graph->updateAllNodes(local_frame, mask_frame);
            graph->updateCondBlfMp(local_frame, last_frame, mask_frame);
            graph->getNextFramePPM(last_frame, mask_frame);
            writeOp(count);
            local_frame = output_frame;
            display();
        }

        count++;
    }
}

void VideoReader::runFC(bool gray, int numLvls, int maxNum, int priorNum, int frameJump) {
    int is_color = ((gray)?0:1);
    int fps = frames_per_second;
    //gray = false;
    //numLvls = 32;
    //maxNum = 200;

    output_video = cvCreateVideoWriter("output.avi", CV_FOURCC('I','4','2','0'), fps, frame_size, is_color);
    difference_video = cvCreateVideoWriter("difference.avi", CV_FOURCC('I','4','2','0'), fps, frame_size, is_color);
    FILE* rFile = fopen("run_main.log", "w");
    fprintf(rFile, "#############################################################\n");
    fprintf(rFile, "# gray: %d NL: %d maxNum: %d\n", gray, numLvls, maxNum);
    fprintf(rFile, "#############################################################\n");
    fclose(rFile);
  
    bool nf=true;
    int key_pressed;
    int count=0;
    IplImage* local_frame=NULL;
    IplImage* inpainted_frame=NULL;
    bool useOpFrame=false;
    double diff= 0.0;
    updateMaskFrame(25, 120, 70, 60);
       
    graph = new VideoGraph(this, gray, numLvls);
    IplImage* result_image = NULL;
    IplImage *empty_mask = NULL;
    allocateOnDemand(&empty_mask, frame_size, IPL_DEPTH_8U, 1); 
    
    while(nf && (count < maxNum)) {
        frame = cvQueryFrame(input_video);
        nf = convertInputFrame(gray, numLvls);

        if(count == 0) {
            graph->initLocals(); 
         //   graph->updateAllNodes(input_frame, mask_frame); 
            initPrevFrames(frameJump);
          //  local_frame = graph->inPaint(input_frame, mask_frame);
            local_frame = cvLoadImage("frame32.png");  
        } else if( (count >= priorNum - frameJump-1) && (count > 0) ){
            updatePrevFrames(local_frame, frameJump);
        }
   
        if(count < priorNum) {
            overlapFrame(local_frame, input_frame, local_frame); 
        } else {
      
            updateBasedOnPrevFrames(empty_mask);
            updateInputFrame(input_frame, mask_frame);  
            graph->getCandidates(local_frame, input_frame, mask_frame); 
          
            graph->PPM_step3(); 
            for(int i=0; i < 3; i++) graph->PPM_inner_loop();
            graph->fillFrame0(); 
            local_frame = output_frame; 
            writeOp(count);
            display();  
        }

        count++;
    }
    

}

void VideoReader::overlapFrame(IplImage* inpainted, IplImage* input, IplImage* res_frame) {
    BwImage mask(mask_frame);
    if(grayscale) {
        BwImage pImg(inpainted);
        BwImage iImg(input);
        BwImage rImg(res_frame);
        for(int i=0; i < getN1(); i++) {
            for(int j=0; j < getN2(); j++) {
                if(mask[i][j] == 0) {
                    rImg[i][j] = iImg[i][j];
                } else {
                    rImg[i][j] = pImg[i][j];
                }
            }
        }
    } else {
        RgbImage pImg(inpainted);
        RgbImage iImg(input);
        RgbImage rImg(res_frame);
        for(int i=0; i < getN1(); i++) {
            for(int j=0; j < getN2(); j++) {
                if(mask[i][j] == 0) {
                    rImg[i][j].r = iImg[i][j].r;
                    rImg[i][j].g = iImg[i][j].g;
                    rImg[i][j].b = iImg[i][j].b;
                } else {
                    rImg[i][j].r = pImg[i][j].r;
                    rImg[i][j].g = pImg[i][j].g;
                    rImg[i][j].b = pImg[i][j].b;
                }
            }
        }
    }

}

void VideoReader::initPrevFrames(int N) {
    prev_frames.clear();
    allocateOnDemand( &diff_frame, frame_size, input_frame->depth, input_frame->nChannels);
    allocateOnDemand( &output_frame, frame_size, input_frame->depth, input_frame->nChannels);
    graph->output_frame = output_frame;
    fi(0, N) {
        IplImage* curr_f = NULL;
        allocateOnDemand( &curr_f, frame_size, input_frame->depth, input_frame->nChannels);
        prev_frames.push_back(curr_f);
    }

}

void VideoReader::runInPaint() {
    int count=1;
    bool gray = false;
    int NL = 32;
    string fname = "/home/vjethava/workspace/eclipse/GBP/data/jumping_girl_input.avi";
    openVideo(fname);
    updateMaskFrame(25, 120, 70, 60);
    //  cvNamedWindow("mask_frame", CV_WINDOW_AUTOSIZE);
    //   cvShowImage("mask_frame", mask_frame);
    //  cvWaitKey(0);
    graph = new VideoGraph(this, gray, NL);
    IplImage* result_image = NULL;
    graph->updateMaskStatus(mask_frame);
    for(int i=0; i < count; i++) {
        frame = cvQueryFrame(input_video);
        convertInputFrame(gray, NL);
        graph->updateAllNodes(input_frame, mask_frame);


        // getchar();
        // cvShowImage("ppm_pred", result_image);
        // cvWaitKey();
    }
    result_image = graph->inPaint(input_frame, mask_frame);

    cvSaveImage("frame.png", result_image);
}


void VideoReader::updateBasedOnPrevFrames(IplImage* mask_frame) {
  //  graph->clearAllNodes();
    graph->clearCondBlfMp(); 
    int N = prev_frames.size(); 
    graph->updateMaskStatus(mask_frame);
    for(int i=0; i < (N-1); i++) {
        cout<<"updateBasedOnPrevFrames() graph->updateCondBlfMp("<<i<<")\n"; 
        graph->updateCondBlfMp(prev_frames[i], prev_frames[i+1], mask_frame);  
    } 
    // graph->updateCondBlfMp(prev_frames[0], prev_frames[1], mask_frame); 
    foreach(iter, graph->adj_map) {
        VideoNode* vn = (VideoNode*) iter->first; 
        vn->localBlfPtr->clear();
        vn->localCondBlfPtr->clear(); 
        for(int i=0; i < (N); i++) {
            StatePtr currSt = graph->readNodeFromFrame(vn, prev_frames[i]);
            LD pastVal = vn->localBlfPtr->getPr(currSt); 
            vn->localBlfPtr->setPr(currSt, max(1.0, pastVal+1.0)); 
            if(i < (N-1)) {
                StatePtr prevSt = graph->readNodeFromFrame(vn, prev_frames[i+1]);
                LD pastVal = vn->localCondBlfPtr->getPr(prevSt, currSt);   
                vn->localCondBlfPtr->setPr(prevSt, currSt, max(1.0, pastVal+1.0) );  
            }
        }
    }
}

void VideoReader::getMissingMask(IplImage* mask_frame, LD probMissing) {
    BwImage mask(mask_frame); 
    srand(time(NULL)); 
    for(int i=0; i < getN1(); i++) {
        for(int j=0; j < getN2(); j++) {
            LD val = ((LD) (rand()%1000))/1000.0;
            if(val < probMissing) {
                mask[i][j] = 255; 
            } else {
                mask[i][j] = 0;
            }
        }
    }
}



void VideoReader::runNoisy(bool gray, int numLvls, int maxNum, int priorNum, int frameJump) {
    int is_color = ((gray)?0:1);
    int fps = frames_per_second;
    output_video = cvCreateVideoWriter("output.avi", CV_FOURCC('I','4','2','0'), fps, frame_size, is_color);
    difference_video = cvCreateVideoWriter("difference.avi", CV_FOURCC('I','4','2','0'), fps, frame_size, is_color);
    FILE* rFile = fopen("run_main.log", "w");
    fprintf(rFile, "#############################################################\n");
    fprintf(rFile, "# gray: %d NL: %d priorNum: %d maxNum: %d\n", gray, numLvls, priorNum, maxNum);
    fprintf(rFile, "#############################################################\n");
    fclose(rFile); 
    graph = new VideoGraph(this, gray, numLvls);
    bool nf=true;
    int key_pressed;
    int count=0;
    IplImage* local_frame=NULL;
    bool useOpFrame=false;
    double diff= 0.0;

   
    graph = new VideoGraph(this, gray, numLvls);
    IplImage* result_image = NULL;
    IplImage* curr_mask = NULL, *empty_mask = NULL;
    allocateOnDemand(&empty_mask, frame_size, IPL_DEPTH_8U, 1); 
    allocateOnDemand(&curr_mask, frame_size, IPL_DEPTH_8U, 1); 
    while(nf && (count < maxNum)) {
        frame = cvQueryFrame(input_video);
        nf = convertInputFrame(gray, numLvls);

        if(count == 0) {
            graph->initLocals(); 
            initPrevFrames(frameJump);
        } else {
            updatePrevFrames(local_frame, frameJump);
        }
   
        if(count < priorNum) {
           local_frame = input_frame; 
        } else {
            getMissingMask(curr_mask, 0.5); 
            updateBasedOnPrevFrames(empty_mask);
            updateInputFrame(input_frame, curr_mask);  
            graph->getCandidates(local_frame, input_frame, curr_mask); 
            graph->PPM_step3(); 
            graph->PPM_inner_loop();
            graph->fillFrame0(); 
            local_frame = output_frame; 
            writeOp(count);
            display();  
        }

        count++;
    }
    
}

void VideoReader::updateInputFrame(IplImage* input_frame, IplImage* curr_mask) {
    static int lcount=0; 
    lcount++; 
    BwImage mask(curr_mask); 
    if(grayscale) {
        BwImage ip(input_frame);
        for(int i=0; i < getN1(); i++) {
            for(int j=0; j < getN2(); j++) {
                if(mask[i][j] != 0) {
                    ip[i][j] = 0; 
                }
            }
        }
    } else {
        IplImage* my_frame= NULL; 
        allocateOnDemand( &my_frame, frame_size, input_frame->depth, input_frame->nChannels);
        
        cvCopy(input_frame, my_frame, 0);
        RgbImage limg(my_frame);  
        RgbImage ip(input_frame); 
        for(int i=0; i < getN1(); i++) {
            for(int j=0; j < getN2(); j++) {
                if(mask[i][j] != 0) {
                   // LD d = ((LD)(rand()%100))/100.0;
                   // if(d < 0.33) 
                    ip[i][j].r = 0;
                    //else if(d < 0.66) 
                    ip[i][j].g = 0; 
                    //else 
                    ip[i][j].b = 0; 
                }
                LD rv = rand()%2, gv = rand()%2, bv=rand()%2;
                if(rv)  limg[i][j].r = 0;
                if(gv)  limg[i][j].g = 0;
                if(bv)  limg[i][j].b = 0;
              
            }
        }
          stringstream ss(""); 
                ss<<"my_frame_"<<lcount<<".png";
                cvSaveImage(ss.str().c_str(), my_frame); 
                cvReleaseImage(&my_frame); 
    }
}  
