#ifndef VR_H_
#define VR_H_

#include <stdio.h>
#include <opencv/cv.h>
#include <math.h>
#include <opencv/highgui.h> 
#include <string>
#include <iostream>
#include "Image.h"
#include "state.h"
/**
 * Class implementation that handles Video I/O
 */
struct VR {
	/// temporary frames used for conversion
    IplImage *frame, *frame1 , *frame1_1C , *frame2_GS ; 
	IplImage* current_frame, *prev_frame, *diff_frame, *res_frame; 
	int number_of_frames;
	int frames_per_second; 
	int frame_num; // current frame number; 
	int Q; // number of quantizations - default 255
	bool is_grayscale; // whether current frame is gray scale or not
	CvSize frame_size; // default frame size under consideration
	CvCapture* input_video;
	int depth; 
	int nChannels;
	inline int getN1() { return frame_size.height; } 
    inline int getN2() { return frame_size.width; }
    CvVideoWriter* output_video;
    CvVideoWriter* difference_video;
    void open_video(string);
    void allocate_frame(IplImage** frame, bool is_grayscale=true);  
	void convert_frame(IplImage* ip_frame, IplImage* op_frame, int _Q, bool _grayscale=true); 
    vector<INT> getDataFromFrame(IplImage* frame, int x, int y, int r, int c);
	void get_frame(vector<vector<INT> > v, IplImage* op_frame, int Q);
	vector<vector<INT> > frame_diff(IplImage* frame_A, IplImage* frame_B, IplImage* diff_frame);    
	void increment_frame(); 
	void write_op(); 
	VR(int _Q=256, bool grayscale=true) {
		Q = _Q;
		is_grayscale =grayscale;  
		frame = NULL; frame1 = NULL;  frame1_1C = NULL; frame2_GS = NULL;  
	 	current_frame=NULL; prev_frame=NULL; diff_frame= NULL, res_frame = NULL;
		output_video = NULL; difference_video = NULL; 
	}
	
	~VR() {
    	if(input_video != NULL) {
        	cvReleaseCapture(&input_video);
    	}
    	cvReleaseVideoWriter(&output_video);
    	cvReleaseVideoWriter(&difference_video);
	}
	inline bool isValidPixel(IplImage* frame, int x, int y) {
    int R = frame->height;
    int C = frame->width;
    return ( (x < R) && (y < C) && (x >=0) && (y >= 0) );
}
	int display(); 
}; 
#endif /*VR_H_*/
