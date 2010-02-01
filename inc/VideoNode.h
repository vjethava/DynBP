#ifndef VIDEONODE_H_    
#define VIDEONODE_H_
#include "node.h" 
#include "GraphPPM.h"
#include <string>
#include <sstream>
using namespace std;
#define NEXT_CAND_MP vector<pair<StatePtr, LD> >
#define CAND_MP_PUSH push_back
class VideoNode : public RegionPPM 
{
    
    /// the global conditional belief pointer of appropriate size
    ConditionalBeliefPtr globalCondPtr; 
    /// the block counting number
    int cr; 
    /// dimensions of the block
    int x, y, l, w;
    /// the global belief pointer of appropriate size;  
    BeliefPtr globalBlfPtr; 
    /// the original state
    StatePtr sOld; 
public:
	//BeliefPtr M, N;
	//BeliefPtr b_bp;  
    BeliefPtr localBlfPtr;
    ConditionalBeliefPtr localCondBlfPtr; 
    /// is masked
    bool isMasked; 
    BeliefPtr getGlobalBlfPtr() { return globalBlfPtr; }
    void setGlobalBlfPtr(BeliefPtr bptr) { globalBlfPtr = bptr; }   
    void setOldState(StatePtr s) { sOld = s; }
    StatePtr getOldState() { return sOld; }
    map<StatePtr, LD, StatePtrCmp> nextStateMp;   
    map<StatePtr, LD, StatePtrCmp> augmentsMp;
    bool addAugmentCand(StatePtr sptr, LD l) {
        if(nextStateMp.find(sptr) == nextStateMp.end()) {
            augmentsMp.insert(make_pair<StatePtr, LD>(sptr, l));
            return true;
        } else
            return false; 
    }
         
    void finalizeCands(StatePtr sO); 
        
    void setNextStateCandMp(NEXT_CAND_MP mp) {
        nextStateMp.clear(); 
        foreach(iter, mp) {
            StatePtr sp  = iter->first; 
            LD val = iter->second; 
            nextStateMp.insert(make_pair<StatePtr, LD>(sp, val));
        } 
        /*foreach(iter, nextStateMp) { 
            cout<<"state: ";
            foreach(it2, iter->first) {
                cout<<" "<<*it2; 
            }
            cout<<" p: "<<iter->second<<"\n";
        }*/ 
    }
    NEXT_CAND_MP getNextStateCandMp() { 
        vector<pair<StatePtr, LD> > result; 
        foreach(iter, nextStateMp) {
            result.push_back(make_pair<StatePtr, LD>(iter->first, iter->second)); 
        }
        return result;     
    }
  
    inline void setCPTR(ConditionalBeliefPtr _cptr) { globalCondPtr = _cptr; }
    inline ConditionalBeliefPtr getCPTR() { return globalCondPtr; } 
	VideoNode(string id);
    virtual ~VideoNode();
    static vector<int> getIntFromId(string s); 
    static bool isInterior(int x, int y, string sBig); 
    static bool isSubset(string sSmall, string sBig);
    static string intersection(string s1, string s2);  
    static string getIdFromInt(int x, int y, int l, int w) {
        stringstream ss(""); 
        ss<<x<<"_"<<y<<"_"<<l<<"_"<<w; 
        return ss.str(); 
    }  
};

#endif /*VIDEONODE_H_*/
