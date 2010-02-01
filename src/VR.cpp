#include "VR.h"
using namespace std; 

void VR::write_op() {
		cout<<"VR::write_OP()\n"; 
		int count = frame_num; 
	    cvWriteFrame( output_video, res_frame );
        cvWriteFrame(difference_video, diff_frame ); 
        stringstream dss(""), oss(""), mss(""); 
        mss<<"input_frame"<<"_"<<count<<".png";
        dss<<"diff_frame"<<"_"<<count<<".png" ;
        oss<<"op_frame"<<"_"<<count<<".png" ;
        cvSaveImage(mss.str().c_str(), current_frame );
        cvSaveImage(dss.str().c_str(), diff_frame );
        cvSaveImage(oss.str().c_str(), res_frame );
}

void VR::increment_frame() { 
	cout<<"VR::increment_frame()\n"; 
	if(frame_num > 0) { 
		cvCopy(current_frame, prev_frame, 0); 
	} 
	cout<<"\tquery frame => ";
	frame = cvQueryFrame(input_video);
	cout<<" done\n\tcvCopy => \n"; 
	ASSERT(current_frame != NULL); 
	// cvCopy(frame, current_frame, 0);
	convert_frame(frame, current_frame, 8, is_grayscale);  
	cout<<" done\n"; 
	frame_num++; 
}

void VR::allocate_frame(IplImage** frame, bool is_grayscale) {
	cout<<"VR::allocate_frame() graycale: "<<is_grayscale<<"\n";
	if(frame != NULL) return ;
	int channels; 
	int depth = IPL_DEPTH_8U; 
	if(is_grayscale) channels = 1; else channels = 3;
	*frame = cvCreateImage(frame_size, depth, channels); 
	if ( *frame == NULL ) {
        fprintf(stderr, "Error: Couldn't allocate image.  Out of memory?\n");
        exit(-1);
    }
}  

void VR::open_video(string videoFileName) {
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

 //   cvSetCaptureProperty( input_video, CV_CAP_PROP_POS_AVI_RATIO, 1. );
  //  number_of_frames = (long) cvGetCaptureProperty( input_video, CV_CAP_PROP_POS_FRAMES );
    //cvSetCaptureProperty( input_video, CV_CAP_PROP_POS_FRAMES, 0. );

    cvSetCaptureProperty( input_video, CV_CAP_PROP_POS_FRAMES, 0. );

    frames_per_second = (int) cvGetCaptureProperty( input_video, CV_CAP_PROP_FPS);
    depth = (int) frame11->depth;
    nChannels = (int) frame11->nChannels; 
    cout<<"NUM: "<<number_of_frames<<"\n";
    printf(" video %s num_frames: %d  fps: %d fheight: %d fwidth: %d depth: %d channels: %d\n", videoFileName.c_str(), number_of_frames, frames_per_second, frame_size.height, frame_size.width, depth, nChannels );
    //fwidth = (int) cvGetCaptureProperty( input_video, CV_CAP_PROP_FRAME_HEIGHT );
    //fheight = (int) cvGetCaptureProperty( input_video, CV_CAP_PROP_FRAME_WIDTH );
    frame_num=0; 
    if(nChannels == 1) is_grayscale = true; else is_grayscale = false;
    
    frame = cvCreateImage(frame_size, depth, nChannels); 
    current_frame = cvCreateImage(frame_size, depth, nChannels); 
    prev_frame = cvCreateImage(frame_size, depth, nChannels); 
    diff_frame = cvCreateImage(frame_size, depth, nChannels); 
    res_frame = cvCreateImage(frame_size, depth, 1); 
    //do frame allocations
    //allocate_frame(&frame, false);
    //allocate_frame(&current_frame, false); //is_grayscale); 
    //ASSERT(current_frame != NULL); 
	//allocate_frame(&prev_frame, false); //is_grayscale);
	//allocate_frame(&diff_frame, false); // is_grayscale); 
	//allocate_frame(&res_frame, true); 
	// open the writers
	output_video = cvCreateVideoWriter("output.avi", CV_FOURCC('I','4','2','0'), 25, frame_size, false);
    difference_video = cvCreateVideoWriter("difference.avi", CV_FOURCC('I','4','2','0'), 25, frame_size, true);
    
}

void VR::convert_frame(IplImage* frame, IplImage* op_frame,  int numLvls, bool gray) {
	if(frame1 == NULL) frame1 = cvCreateImage(frame_size, depth, nChannels);
    cvConvertImage(frame, frame1, 0);
   // allocate_frame(op_frame, _grayscale);
    if(gray) {
        if(frame1_1C == NULL) frame1_1C = cvCreateImage(frame_size, IPL_DEPTH_8U, 1 );
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
        // frame1_1C;
        
        cvCopy(frame1_1C, op_frame, 0); 
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
       cvCopy(frame1, op_frame, 0);
    }
}


vector<vector<INT> > VR::frame_diff(IplImage* frame1, IplImage* frame2, IplImage* resFrame) {
	cout<<"VR::frame_diff()\n"; 
	int N = frame1->width*frame->height;
	int N1 = frame1->height; 
	int N2 = frame1->width;  
	vector<vector<INT> >  result;
	vector<INT> row(N2, 1); fi(0, N1) {  result.push_back(row); }
	int div = 1; 
    if(is_grayscale) {
        if(resFrame == NULL) {
            resFrame = cvCreateImage(cvSize(frame1->width, frame1->height), IPL_DEPTH_8U,1);
        }
        N = frame1->height*frame1->width;
        BwImage img1(frame1);
        BwImage img2(frame2);
        BwImage resImg(resFrame);
        fi(0, frame1->height) {
            fj(0, frame1->width) {
            	int idx = i*getN2() + j; 
                int diff = (int) abs( ((int) (img1[i][j]))-((int) (img2[i][j])) );
                resImg[i][j] = (INT) diff;
  cout<<"resImg["<<i<<"]["<<j<<"] = "<<resImg[i][j]<<diff<<" = "<<diff<<"\n" ; getchar(); 
                if(img1[i][j] == img2[i][j]) {
                	result[i][j] = (INT) 0; 
                	resImg[i][j] = (INT) 0; 
                    //*NP += 1;
                } else {
                	resImg[i][j] = (INT) 255; 
                	result[i][j] = (INT) 1; 
                }
                //res += p*p;
            }
        }
        //res = sqrt(res/N);
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
                int tot = 0; 
                
             //    cout<<"resImg["<<i<<"]["<<j<<"] = "<<resImg[i][j].r<<diff<<" = "<<diff<<"\n" ; getchar(); 
                //int idx = i*getN2() + j;
                diff[0] = (int) abs( ((int) (img1[i][j].r))-((int) (img2[i][j].r)) );
                diff[1] = (int) abs( ((int) (img1[i][j].g))-((int) (img2[i][j].g)) );
                diff[2] = (int) abs( ((int) (img1[i][j].b))-((int) (img2[i][j].b)) );
                resImg[i][j].r = (INT) diff[0];
                resImg[i][j].g = (INT) diff[1];
                resImg[i][j].b = (INT) diff[2];
                tot = diff[0] + diff[1] + diff[2]; 
                if( (img1[i][j].r == img2[i][j].r) && (img1[i][j].g == img2[i][j].g) && (img1[i][j].b == img2[i][j].b)) {
                	result[i][j] = (INT) 0; 
                	resImg[i][j].r = (INT) 0; 
                	resImg[i][j].g = (INT) 0;
                	resImg[i][j].b = (INT) 0;
                } else {
                	resImg[i][j].r = (INT) 255; 
                	resImg[i][j].g = (INT) 255;
                	resImg[i][j].b = (INT) 255;
                	result[i][j] = (INT) 1;
                }
          //      res += p*p;
            }
        }
    //    res = sqrt(res/N);

    }
   	return result; 
}

vector<INT> VR::getDataFromFrame(IplImage* frame, int x, int y, int r, int c) {
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
	
void VR::get_frame(vector< vector<INT> > v, IplImage* op_frame,int Q) {
	int N1v = v.size(); 
	int N2v = v[0].size(); 
	ASSERT( (getN1() == N1v) && (getN2() == N2v) ); 
	cout<<"VR::get_frame()\n"; 
	int mult = 256/Q; 
	//if(isgrayscale) {
		BwImage img(op_frame); 
		for(int i=0; i < getN1(); i++) {
			for(int j=0; j < getN2(); j++) {
				int idx = i*getN2() + j;
				img[i][j] = v[i][j]*256/Q; 
			}
		}
//	}
}

int VR::display() {
	cout<<"VR::display()\n"; 
	const char* wname="input_frame"; 
	cvNamedWindow(wname, CV_WINDOW_AUTOSIZE); 
    cvShowImage(wname, current_frame);
   
    const char* wname2= "diff_frame"; 
    cvNamedWindow(wname2, CV_WINDOW_AUTOSIZE); 
    cvShowImage(wname2, diff_frame);
    
    const char* wname3= "res_frame"; 
    cvNamedWindow(wname3, CV_WINDOW_AUTOSIZE); 
    cvShowImage(wname3, res_frame);

	int x = 0; // x = cvWaitKey(10);
	return x;  
}
