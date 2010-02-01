#include <cmath>
#include <iostream>
#include <cassert>
#include "state.h"
LD State::SSD(State& s1, State& s2) {
    assert(s1.num == s2.num);
    LD res = 0.0;
    for(int i=0; i < s1.num; i++) {
        int diff = ((int) s1.getData(i)) - ((int) s2.getData(i)); 
        res += diff*diff; 
    }
    res = sqrt(res/((LD) s1.num));
    return res;
}

