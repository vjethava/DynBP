#include "Tracker.h"
#include <iostream> 
#include "main.h"

using namespace std; 

int main(int argc, char* argv[]) {
    ios::sync_with_stdio(true);
    string fname; 
    if (argc != 2)
    {
        fname = "/home/vjethava/workspace/eclipse/GBP/data/jumping_girl_input.avi"; 
    } else 
        fname = string(argv[1]);  
    bool gray=true; 
    int NL=4; 
    int priorNum=60;
    int maxNum=120;
    int frameJump=3; 
    int C = 2;
    LD theta = 0.5;  
    Tracker tracker; // = new Tracker();
    tracker.initMain(fname, NL, gray, C, theta, 2, 2, 1, 1);  
    tracker.run_step(1.0, 0.5); 
    return 0; 
}

