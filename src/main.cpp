#include "main.h"
//namespace common {
using namespace std; 
/**
 * Return Q^N
 */
long MyDef::pow(int Q,int N) {
		long res = 1;
		for(int i=0; i < N; i++) { res *= (long) Q; }
		return res; 
}

vector<string> split(string ip, char delim) {
    vector<string> result;
    int pos = 0;
    while(ip.find(delim) != string::npos) {
        pos = ip.find(delim);
        string curr = ip.substr(0, pos);
        ip = ip.substr(pos+1);
        result.push_back(curr);
    }
    ASSERT(ip.size() > 0); // last char is not a delim
    result.push_back(ip); // last part
    return result;
}

bool Cluster::subset(const Cluster& c1, const Cluster& c2) {
    vector<string> f1 = c1.first;
    vector<string> f2 = c2.first;
    if(f1.size() > f2.size())
        return false;
    else {
        foreach(iter, f1) {
            int i=0;
            bool cfound = false;
            while( (i < (int) f2.size()) && (!cfound)) {
                if(*iter == f2[i])
                    cfound = true;
                i++;
            }
            if(cfound == false)
                return false;
        }
    }
    return true;
}

bool Cluster::intersection(Cluster* cc, const Cluster& c1, const Cluster& c2) {
    vector<string> cvars;
    vector<string> cfctrs;
    // cout<<c1;
    // cout<<c2;
    foreach(iter1, c1.first) { // look at c1.vars
        foreach(iter2, c2.first) {
            if((*iter1) == (*iter2)) {
                cvars.push_back((*iter1));
                //            PRINTF("common::intersection(varMatch) %s %s\n", (*iter1).c_str(), (*iter2).c_str());
            }
        }
    }

    if(cvars.size() == 0)
        return false;
    foreach(iter1, c1.second) { // look at c1.factors
        foreach(iter2, c2.second) {
            if((*iter1) == (*iter2)) {
                cfctrs.push_back((*iter1));
            }
        }
    }
    sort(cvars.begin(), cvars.end());
    sort(cfctrs.begin(), cfctrs.end());
    cc->first = cvars;
    cc->second = cfctrs;
    return true;
}

//};
