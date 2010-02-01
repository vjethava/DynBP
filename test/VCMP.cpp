#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
using namespace std;
#define foreach(it,c) for(typeof((c).begin()) it=(c).begin();it!=(c).end();++it)
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

int main() {
    srand(time(NULL)); 
    vector<int> v1, v2;
    for(int i=0; i < 10; i++) {
        v1.push_back(rand()%100);
        v2.push_back(rand()%100);
    }
    
    sort(v1.begin(), v1.end()); 
    foreach(iter, v1) cout<<*iter<<" "; 
    /*
    map<vector<int>, double, VCMP<int> > mp;
    mp.insert(make_pair<vector<int>, double>(v1, 0.9));
    mp.insert(make_pair<vector<int>, double>(v2, 0.9));
    foreach(miter, mp) {
        cout<<"v: {";
        foreach(viter, miter->first) cout<<" "<<*viter;
        cout<<"} p: "<<miter->second<<"\n"<<flush;
    }*/ 
    return 0;
}
