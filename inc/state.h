#ifndef STATE_H_
#define STATE_H_
#include <iostream>
#include "main.h"
#include <vector>
using namespace std;
class State;
typedef State* StatePtr; 

class State {
public:
    State(int _n, INT* _d) { 
        num=_n; 
        data= new INT[num]; 
        for(int i=0; i < num; i++) {
            *(data+i) = *(_d+i);     
        } 
    }
    State(vector<INT> v) {
    	num = v.size(); 
    	data = new INT[num];
    	for(int i=0; i < num; i++) {
    		*(data+i) = v[i]; 
    	}
    }
    State(const State& st2) {
        num = st2.num;
        data = new INT[num]; 
        for(int i=0; i < num; i++) {
           // *(data+i) = *(st2.data+i);
           data[i] = st2.data[i];  
        }
    }
    static LD SSD(State& s1, State& s2);
    int num;
    INT *data;
    bool operator==(const State& s2) const {
        bool res = true;
        if(num != s2.num)
            res = false;
        else {
            int i=0;
            while(res && (i < num)) {
                res = (*(data+i) == *(s2.data+i));
                i++;
            }
        }
        return res;
    }
    inline INT getData(int idx) const { return *(data+idx); } 
       
    friend ostream& operator<<(ostream& os, State& s) {
        //INT s0 = 0; 
        //cout<<"s.num: "<<s.num<<"\n"; fflush(stdout); 
        os<<"[";
        int i=0; 
        while(i < s.num) {
           // os<<((int) (*(s.data+i));
            os<<((int) s.data[i]);
            if(i < (s.num-1)) os<<" ";
            i++; 
        }
        os<<"]";
        return os;
    }
    State operator+(const State& s2) const {
        int n = num + s2.num;
        INT* data = new INT[n];  
        for(int i=0; i < n; i++) {
            if(i < num) {
                *(data+i) = getData(i); 
            } else {
                *(data+i) = s2.getData(i-num); 
            }
        }
        State st = State(n, data);
        delete(data); 
        return st;  
    }
    ~State() {
        delete(data); 
    }
    
    VI getVec() { 
        VI result; 
        
        for(int i=0; i < num; i++) {
            result.push_back(data[i]);
        }
        return result; 
    }
};

struct StatePtrCmp {
    bool operator()(const StatePtr& sp1, const StatePtr& sp2) const {
        State s1 = *(sp1); 
        State s2 = *(sp2); 
        bool res=false; 
        if(s1.num < s2.num) res=true;
        else if(s1.num > s2.num) res = false;
        else {
            bool found = false;
            int i=0; 
            while((!found) && (i < s1.num)) {
                if(*(s1.data + i) < *(s2.data + i) ) {
                    res = true; found = true; 
                } else if(*(s1.data + i) > *(s2.data + i) ) {
                    res = false; found = true; 
                } 
                i++; 
            } 
        }
        return res; 
    }
};

struct SPLD_Cmp {
    bool choice; 
    SPLD_Cmp(bool ch=true) {
        choice = ch; 
    }
    bool operator()(const pair<StatePtr, LD>& p1,  const pair<StatePtr, LD>& p2) const {
        if(choice) return (p1.second > p2.second);
        else return (p1.second < p2.second); 
    }
};
#endif /*STATE_H_*/

