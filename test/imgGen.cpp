#include <iostream>
#include <cv.h>
#include <highgui.h>
#include <cmath>
#include <sstream>
#include <fstream>
#include <cstdio>   
#include "Image.h" 
using namespace std; 

int main() {
	srand(0); 
	int N1 = 50, N2= 50;
	int P1 = 5, P2 = 5;  
	int NF = 25;
	int bk[N1][N2]; 
	for(int i=0; i < N1; i++) for(int j=0; j < N2; j++) bk[i][j] = rand()%256; 
	int pch[P1][P2]; 
	for(int i=0; i < P1; i++) for(int j=0; j < P2; j++) pch[i][j] = rand()%256; 
	CvSize frame_size; 
	frame_size.height = N1; 
	frame_size.width = N2;
	IplImage* frame = NULL;
	frame = cvCreateImage(frame_size, IPL_DEPTH_8U, 1); // grayscale
	cout<<" h: "<<frame->height<<" w: "<<frame->width<<"\n"; getchar(); 
	BwImage img(frame);
	int px = 5, py = 5;  
	for(int i=0; i < NF; i++) {
		stringstream ssf(""); 
		ssf<<"frame_"<<i<<".png"; 
		for(int j=0; j < N1; j++) {
			for(int k=0; k < N2; k++) {
				cout<<"j: "<<j<<" k: "<<k<<endl;
				img[j][k] = bk[j][k]; 
			}
		}
		for(int j=0; j < P1; j++) {
			for(int k=0; k < P2; k++) {
				int x = j + px; 
				int y = k + py; 
				img[x][y] = pch[j][k]; 
			}
		}
		px++;
		py++;
		cvSaveImage(ssf.str().c_str(), frame); 
	}
	return 0; 
}
