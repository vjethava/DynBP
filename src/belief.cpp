#include "belief.h"
#include <cmath>
using namespace std;

	
void Belief::setPr(StatePtr _s, LD _v) {
    //a cout<<"BLF.setPr("<<_s<<", "<<_v<<")\n";
    map<StatePtr, LD>::iterator itCurr = bMp.find(_s);
    
    if(INT_MODE) { // save frequency
        if(itCurr == bMp.end()) {
            bMp.insert(make_pair<StatePtr, LD>(_s, 1.0));
        } else
            itCurr->second = itCurr->second + 1.0;
        incrTot(); // update total frequency
    } else {

        if(itCurr != bMp.end()) {
            //        cout<<" reassigning ";
            //        StatePtr s = (StatePtr) itCurr->first;
            //         cout<<s;
            //         cout<<"\n";
            allocPr -= itCurr->second;
            bMp.erase(itCurr);
        }
        bMp.insert(make_pair<StatePtr, LD>(_s, _v));
        allocPr += _v;
        //       assert(allocPr <= (1.0 + 1E-9));
    }
}

LD Belief::getPr(StatePtr s) {
    LD result=0.0;
    map<StatePtr, LD>::iterator itCurr = bMp.find(s);
    if(itCurr != bMp.end()) {
        result = itCurr->second;
        if(INT_MODE) { // get probability from frequency
            result = result/((LD) (totalFreq));
        } else {
            result = result;
        }
    } else
        result = BELIEF_NOT_FOUND;

    return result;
}

void ConditionalBelief::setPr(StatePtr sC, StatePtr sN, LD _cbp, LD _origp) {
    //  cout<<"CB::setPr("<<sC<<", "<<sN<<", "<<_cbp<<", "<<_origp<<")\n";
    if(_origp > 0.0) {
        Belief* pe= this;
        pe->setPr(sC, _origp);
    }
    if(cbMp.find(sC) == cbMp.end()) {
        //    cout<<"\tInserted new Belief";
        BeliefPtr bN = new Belief();
        cbMp.insert(make_pair<StatePtr, Belief>(sC, *bN));
    }
    cbMp.find(sC)->second.setPr(sN, _cbp);
}

LD ConditionalBelief::getPr(StatePtr sC, StatePtr sN) {
    LD result;
    map<StatePtr, Belief>::iterator itCurr = cbMp.find(sC);
    if(itCurr != cbMp.end()) {
        result = itCurr->second.getPr(sN);
    } else
        result = BELIEF_NOT_FOUND;
    return result;
}

LD ConditionalBelief::getAllocPr(StatePtr sC) {
    LD result = BELIEF_NOT_FOUND;
    map<StatePtr, Belief>::iterator itCurr = cbMp.find(sC);
    if(itCurr != cbMp.end()) {
        result = itCurr->second.getAllocPr();
    }
    return result;
}

LD ConditionalBelief::getTotalAlloc() {
    LD result = 0.0;
    foreach(it, cbMp) {
        result += it->second.getAllocPr();
    }
    return result;
}

/*
void Belief::clear() {
    totalFreq = 0;
    bMp.clear(); 
    allocPr = 0.0; 
}
*/
