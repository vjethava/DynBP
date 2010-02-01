#include "VideoReader.h"
#include "VideoNode.h" 
#include "VideoGraph.h" 
#include "main.h"
#include <iostream>
using namespace std; 


int main(int argc, char* argv[]) {
    ios::sync_with_stdio(true);
   VideoReader vReader;
   //VideoNode("99_0_1_5"); getchar(); 
    string fname; 
    if (argc != 2)
    {
        fname = "/home/vjethava/workspace/eclipse/GBP/data/jumping_girl_input.avi"; 
    } else 
        fname = string(argv[1]);  
    bool gray=false; 
    int nl=32; 
    int priorNum=60;
    int maxNum=120;
    int frameJump=3; 
    
    cout<<"grayscale: "; cin>>gray; 
    cout<<"numLvls: "; cin>>nl;
    cout<<"priorNum: "; cin>>priorNum;  
    cout<<"maxNum: "; cin>>maxNum;  
    cout<<"frame_jump: "; cin>>frameJump; 
    
     
    vReader.openVideo(fname);
    // vReader.runMain(gray, nl, maxNum, priorNum, frameJump);  
    vReader.runNoisy(gray, nl, maxNum, priorNum, frameJump);  
    //vReader.runFC(gray, nl, maxNum, priorNum, frameJump); 
 
    return 0; 
}
 
/*
int main(int argc, char* argv[]) {
    ios::sync_with_stdio(true); 
    VideoReader vReader; 
     string fname; 
    if (argc != 2)
    {
        fname = "/home/vjethava/workspace/eclipse/GBP/data/jumping_girl_input.avi"; 
    } else 
        fname = string(argv[1]);
   vReader.openVideo(fname);
  // vReader.runInPaint(); 
    
}

 
*/ 
