#ifndef BELIEF_H_
#define BELIEF_H_

#include <iostream>
#include <map>
#include <string>
#include <iomanip>
#include <sstream>
#include <cassert>
#include "state.h"
#include "main.h"
//typedef std::string State;
using namespace std;
const LD BELIEF_NOT_FOUND=(LD)-1.0;

class Belief;
typedef Belief* BeliefPtr;

/// class that stores all possible the beliefs
class Belief {

    /// the total probability mass allocated in bMp
    LD allocPr;
    /// map having actual belief - probability pairs
    //std::map<std::string, LD> bMp;
    map<StatePtr, LD, StatePtrCmp> bMp;
    /// variable used in tracking true_pr = act_pr/totalFreq ( to make compatible with int mode )
    int totalFreq;
    inline void incrTot() {
        totalFreq = totalFreq + 1;
    }
public:
    // void clear();
    map<StatePtr, LD, StatePtrCmp> getMp() const {
        return bMp;
    }
    /*map<StatePtr, LD, StatePtrCmp>* getMpPtr() const {
        return (&(bMp));
    }*/
    /// variable that is true when the map stores the frequency
    bool INT_MODE;
    /// the default constructor - note that INT_MODE must be initialized before call to this
    Belief(bool intMode=false) {
        INT_MODE=intMode;

        allocPr = (LD)0.0;
        if(!INT_MODE) {
            totalFreq = 1;
        } else {
            totalFreq = 0;
        }
    }

	//Belief(int Q, int numVar, bool isRandom, int seed=0); 
    /// returns the total allocated probability ( relevant to LD_MODE )
    inline LD getAllocPr() const {
        return allocPr;
    }
    /// returns the totalFreq (INT_MODE)
    inline int getTotalFreq() const {
        return totalFreq;
    }
    /** sets the probability of given state ( if INT_MODE = false)
     * increments the frequency of given state and the totalFreqency ( if INT_MODE = true )
     */
    void setPr(StatePtr _s, LD _v=-1.0);
    /// returns the probability of given state(if found) or BELIEF_NOT_FOUND=-1.0
    LD getPr(StatePtr s);
    ~Belief() { }
    void clear() {
        allocPr = (LD) 0.0;
        totalFreq = 0;
        bMp.clear();
    }
    /// normalizes allocPr = 1.0
    void normalize() {
        ASSERT(INT_MODE == false);
        ASSERT(allocPr != 0.0);
        if(allocPr != 1.0) {
            foreach(iter, bMp) {
                iter->second = iter->second/allocPr;
            }
        }
        allocPr = (LD) 1.0;
    }
    friend ostream& operator<<(ostream& os, Belief& b) {
        std::map<StatePtr, LD>::iterator itCurr=b.bMp.begin();
        std::map<StatePtr, LD>::iterator itEnd=b.bMp.end();
        os<<"{ ";
        while(itCurr != itEnd) {
            State s = (State) *(itCurr->first);
            os<<s;
            os<<"="<<b.getPr(itCurr->first)<<" ";
            itCurr++;
        }
        if(b.INT_MODE) {
            os<<" tf: "<<b.totalFreq<<"} ";
        } else {
            os<<" ap: "<<b.allocPr<<"} ";
        }
        return os;
    }
};

/// class that stores the conditional beliefs
class ConditionalBelief: public Belief {

    std::map<StatePtr, Belief, StatePtrCmp> cbMp;
public:
    void clear() {
        ((BeliefPtr) this)->clear(); 
        foreach(iter, cbMp) {
            iter->second.clear(); 
        }
        cbMp.clear(); 
    }
    std::map<StatePtr, Belief, StatePtrCmp>* getCbMpPtr() { return (&cbMp); } 
    BeliefPtr getBeliefForState(StatePtr sptr) {
        std::map<StatePtr, Belief, StatePtrCmp>::iterator miter = cbMp.find(sptr);
        if(miter == cbMp.end())
            return 0;
        return &(miter->second);
    }
    /** sets the probability of finding P(next state | current state)
     *  P(sN | sC) = _cbp
     *  P(sC) = _origp
     */
    void setPr(StatePtr sC, StatePtr sN, LD _cbp, LD _origp=-1.0);
    /// returns the probability p(sN | sC) or BELIEF_NOT_FOUND
    LD getPr(StatePtr sC, StatePtr sN);
    /// get the allocation corresponding to belief starting at sC or BELIEF_NOT_FOUND
    LD getAllocPr(StatePtr sC);
    /// the total allocatio present
    LD getTotalAlloc();
    friend ostream& operator<<(ostream& os, ConditionalBelief& b) {
        std::map<StatePtr, Belief>::iterator itCurr=b.cbMp.begin();
        std::map<StatePtr, Belief>::iterator itEnd=b.cbMp.end();
        os<<"{ ";
        while(itCurr != itEnd) {
            State s = *(itCurr->first);
            os<<"\t";
            os<<s;
            os<<" => "<<itCurr->second<<"\n";
            itCurr++;
        }
        os<<"}\n";
        return os;
    }
};

typedef ConditionalBelief* ConditionalBeliefPtr;
#endif /*BELIEF_H_*/

