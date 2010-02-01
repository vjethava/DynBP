#ifndef VIDEOREADER_H_
#define VIDEOREADER_H_
#include <stdio.h>
#include <opencv/cv.h>
#include <math.h>
#include <opencv/highgui.h> 
#include <string>
#include <iostream>
#include "Image.h"
#include "state.h"
using namespace std; 
void allocateOnDemand( IplImage **img, CvSize size, int depth, int channels);
class VideoGraph; 
 
/**
 * This class handles all interactions with the video data 
 */
class VideoReader {

    /// frame count 
    int frame_count;
   
    /// whether to do the inference on bw or color
    bool grayscale;   
    /// the maximum value each pixel can take 0...maxgrayscale (default 255)
    int numLvls;  
    /// biggest region width 
    int rwidth;
    /// biggest region height
    int rheight;
    /// horizontal overlap between adjoining biggest regions
    int xoverlap;
    /// vertical overlap between adjoining biggest regions
    int yoverlap; 
   
    
    /// pointer to the difference frame
    IplImage* diff_frame; 
    void* last_frame_window_handle;
    /// pointer to the default window
    void* input_window_handle; 
    /// pointer to the output window 
    void* output_window_handle; 
    /// pointer to the mask frame
    void* mask_window_handle; 
    /// number of frames in the video sequence
    long number_of_frames; 
    ///number of frames per second
    int frames_per_second; 
   
   
    /// the graph used for inference 
    VideoGraph* graph; 
    CvVideoWriter* output_video;
    CvVideoWriter* difference_video;
    vector<IplImage* > prev_frames;  
    int count11; 
public:
 	/// temporary frames used for conversion
    IplImage *frame, *frame1 , *frame1_1C , *frame2_GS ; 
    /// frame size of the video 
    CvSize frame_size;
     /// the last frame 
    IplImage* last_frame;
    /// the default pointer to the current frame
    IplImage* input_frame;
    /// pointer to the output frame
    IplImage* output_frame;
   
     /// pointer to the mask frame
    IplImage* mask_frame;
     /// the default pointer to the current video 
    CvCapture* input_video; 
    inline void setGraph(VideoGraph* _g) { graph = _g; } 
    inline int getN1() { return frame_size.height; } 
    inline int getN2() { return frame_size.width; }
    void inline setGrayScale(bool _isGray) { grayscale = _isGray; }     
    void inline setNumLvls(int s) { numLvls=s; }
    VideoReader(bool _gray=true, int _numLvls=8); 
    /// returns the current frame
    IplImage* getCurrentFrame(CvCapture* currVideo);
    /// opens a video 
    void openVideo(string videoName);
    /// the main engine that updates all the three engines
    void runMain(bool gray=true, int numLvls=4, int maxNum=-1, int priorNum=5, int frameJump=1); 
    /// this function just displays the next frame of the input video until no more to display    
    bool convertInputFrame(bool gray=true, int numLvls=4); 
    /// this function is recursively called until all the video holes have been filled 
    void updateOutputVideo();  
    virtual ~VideoReader(); 
    /// get the data for block of pixels 
    vector<INT> getDataFromFrame(IplImage* frame, int x, int y, int r, int c);
    /// checks for validity of pixel in frame size 
    bool isValidPixel(IplImage* frame, int x, int y);
    /// updates the last frame
    void updateLastFrame(IplImage* frame = NULL); 
    void display(); 
    double getFrameDifference(IplImage* frame1, IplImage* frame2, IplImage* resFrame=NULL, int* NP=NULL);
    void initPrevFrames(int N); 
    void updateMaskFrame(int x, int y, int l, int w);
    void updatePrevFrames(IplImage* frame, int N) {
        for(int i=0; i < (N-1); i++) {
            cvCopy(prev_frames[i], prev_frames[i+1], 0);
        }
        cvCopy(frame, prev_frames[0], 0); 
        last_frame = prev_frames[N-1]; 
    }
    void overlapFrame(IplImage* inpainted, IplImage* input, IplImage* res_frame);
    void runFC(bool gray, int numLvls, int maxNum, int priorNum, int frameJump);
    // void overlapFrame(IplImage* inpainted, IplImage* input);
    void runInPaint();
    void getMissingMask(IplImage* mask_frame, LD probMissing);
    void updateBasedOnPrevFrames(IplImage* mask_frame); 
    void runNoisy(bool gray, int numLvls, int maxNum, int priorNum, int frameJump);
    void updateInputFrame(IplImage* input_frame, IplImage* curr_mask);
    bool incrementFrame(IplImage* pastFrame, IplImage* nextFrame) {
        allocateOnDemand(&last_frame, frame_size, IPL_DEPTH_8U, 1); 
        bool res = true; 
        if(count11 > 0) {
            cvCopy(input_frame, last_frame, 0);
            frame = cvQueryFrame(input_video);
            res = convertInputFrame(grayscale, numLvls);
            pastFrame = last_frame;
            nextFrame = input_frame; 
          
        } else {
            frame = cvQueryFrame(input_video);
            res = convertInputFrame(grayscale, numLvls);
        } 
        count11++; 
        return res; 
    }
    void writeOp(int count) {
        int NP = 0; 
        LD diff = getFrameDifference(input_frame, output_frame, diff_frame, &NP); 
        cvWriteFrame( output_video, output_frame );
        cvWriteFrame(difference_video, diff_frame ); 
        stringstream dss(""), oss(""), mss(""); 
        mss<<"input_frame"<<"_"<<count<<".png";
        dss<<"diff_frame"<<"_"<<count<<".png" ;
        oss<<"op_frame"<<"_"<<count<<".png" ;
        cvSaveImage(mss.str().c_str(), input_frame );
        cvSaveImage(dss.str().c_str(), diff_frame );
        cvSaveImage(oss.str().c_str(), output_frame );
        FILE* rFile2 = fopen("run_main.log", "a+"); 
        fprintf(rFile2, "frame_num: %d diff: %f pixel_diff: %d\n", count, diff, NP);  
        fclose(rFile2); 
    }
};

#endif /*VIDEOREADER_H_*/
