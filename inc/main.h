#ifndef MAIN_H_
#define MAIN_H_
#include <ctime>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <fstream>
#include <cassert>
#include <cstdio> 
#include <algorithm>
using namespace std; 
 
#define fi(a, b) for(int i=a; i < b; i++)
#define fj(a, b) for(int j=a; j < b; j++)
#define fk(a, b) for(int k=a; k < b; k++)
#define foreach(it,c) for(typeof((c).begin()) it=(c).begin();it!=(c).end();++it)
#define VS vector<string>
#define VPII vector<pair<int, int> >
#ifdef DEBUG
    #define PRINTF printf
    #define FPRINTF fprintf
    #define ASSERT assert
#else
    #define PRINTF //
    #define FPRINTF //
    #define ASSERT //
#endif    
typedef double LD; 
typedef int INT;
typedef vector<INT> VI;

namespace MyDef {
	long pow(int Q, int N);
}

template<typename T> string getStr(vector<T> v) {
    stringstream ss("");
    ss<<"{"; 
    foreach(iter, v) {
        ss<<" "<<*iter;  
    }
    ss<<"}"; 
    return ss.str(); 
}

struct Cluster {
    vector<string> first; // vars
    vector<string> second; // factors
    Cluster(vector<string> vars, vector<string> fctrs) {
        sort(first.begin(), first.end()); 
        sort(second.begin(), second.end()); 
        first = vars; 
        second = fctrs; 
    
    } 
    
    Cluster() { ;  }
    
    friend ostream& operator<<(ostream& os, const Cluster& c) {
        os<<"[var: {"; foreach(iter, c.first) { os<<" "<<(*iter); }
        os<<"} fct: {"; foreach(iter, c.second) { os<<" "<<(*iter); }
        os<<"}] "; 
        return os; 
    }
    /// returns c1 subset c2
    static bool subset(const Cluster& c1, const Cluster& c2);
    bool operator==(const Cluster& c) const {
        if(first.size() != c.first.size()) return false; 
        else if(second.size() != c.second.size()) return false; 
        for(int i=0; i < (int)first.size(); i++) {
            if(first[i] != c.first[i]) return false; 
        }
        for(int i=0; i < (int)second.size(); i++) {
            if(second[i] != c.second[i]) return false; 
        }
        return true;
    }
    string getRegionId() const {
        stringstream ss(""); 
        ss<<"R"; 
        foreach(iter, second) {
            ss<<" "<<*iter; 
        }
        return ss.str(); 
    }
    bool operator!=(const Cluster& c) const {
        if( (*this) == c ) return false; 
        else return true; 
    }
    static bool intersection(Cluster* cc, const Cluster& c1, const Cluster& c2); 
    
};
/// splits the string based on delimiters
vector<string> split(string ip, char delim);

struct ClusterCmp {
    bool operator()(const Cluster& c1, const Cluster& c2) const {
        /*cout<<"clusterCmp()\n\t"<<c1;
        cout<<"\n\t"<<c2;
        cout<<"\n"; */
        vector<string> f1 = c1.first;
        vector<string> f2 = c2.first;
        if(f1.size() == f2.size()) {
            for(int i=0; i < (int) f1.size(); i++) {
                if(f1[i] != f2[i]) return (f1[i] < f2[i]); 
            }
        } else 
            return (f1.size() < f2.size()); 
        return false; 
    }

    
};



template <typename T, typename C> struct PRCMP2 {
    bool choice; 
    PRCMP2(bool _choice = true) { choice = _choice; }  
    bool operator() (const pair<T, C>& p1, const pair<T, C>& p2) const {
        if(choice) return p1.second < p2.second;
        else return p1.second > p2.second;  
    }
};

template <typename T, typename C> struct PRCMP {
    bool operator() (const pair<T, C>& p1, const pair<T, C>& p2) const {
        if(p1.first == p2.first) return p1.second < p2.second;
        else return p1.first < p2.first; 
    }
};

template <typename T> struct VCMP {
    bool operator ()(const vector<T>& v1, const vector<T>& v2) const {
        if(v1.size() != v2.size())
            return (v1.size() < v2.size());
        else
            for(int i=0; i < v1.size(); i++) {
                if(v1[i] != v2[i])
                    return v1[i] < v2[i];
            }
        return false;
    }
};

class Timer {
    time_t tb;
    double last;
    int incr;
    FILE* fout   ;
public: 
    Timer(int _incr=1, FILE* _fout = stderr) {
        tb = time(NULL); 
        incr = _incr; 
        fout = _fout; 
    }
    void tick() {
        time_t tc = time(NULL); 
        double elapsed = difftime(tc, tb) - last; 
        if(elapsed - last > incr) {
            last += incr;
            FPRINTF(fout,"."); 
        }
    }
    double getElapsed() { return last; } 
}; 
     
#endif /*MAIN_H_*/
