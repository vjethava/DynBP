#include "VideoNode.h"
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <cmath>
using namespace std;
VideoNode::VideoNode(string id) :
	RegionPPM(id, 0, 0) {
	BeliefPtr b = new Belief();
	setBelief(b);
	globalBlfPtr = 0;
	globalCondPtr = 0;
	vector<int> nums = VideoNode::getIntFromId(id);
	x = nums[0];
	y = nums[1];
	l = nums[2];
	w = nums[3];
	sOld = 0;
	isMasked = false;
	setNodeNumVar(l*w);
	int N = l*w;
	// ncBlfPtr = new BeliefPtr();
	// algoBlfPtr = new BeliefPtr();
	// cout<<"x: "<<x<<" Y: "<<y<<" l: "<<l<<" w: "<<w<<"N: "<<N<<"\n";
	// M = 0; N = 0; b_bp = 0; 
}

VideoNode::~VideoNode() {
}

vector<int> VideoNode::getIntFromId(string id) {
	vector<string> nums = split(id, '_');
	assert(nums.size() == 4);
	int x;
	vector<int> result;
	fi(0, 4) {
		stringstream ss(nums[i]);
		ss>>x;
		result.push_back(x);
	}
	return result;
}

bool VideoNode::isSubset(string sSmall, string sBig) {
	vector<int> bounds = getIntFromId(sSmall);
	bool res;
	if (VideoNode::isInterior(bounds[0], bounds[1], sBig) && isInterior(
			(bounds[0]+bounds[2]-1), (bounds[1]+bounds[3]-1), sBig) )
		res = true;
	else
		res = false;
	return res;
}

/**
 * returns the block that is common between s1 and s2
 */
string VideoNode::intersection(string s1, string s2) {
//	cout<<"intersection( s1: "<<s1<<", s2: "<<s2<<")\n";
	vector<int> a = getIntFromId(s1);
	vector<int> b = getIntFromId(s2);
	string result="";
	vector<int> c(4,0);
	if(isInterior(a[0], a[1], s2)) {
		c[0] = a[0]; 
		c[1] = a[1]; 
		int x = a[0], y = a[1]; 
		int l=0, w=0;  
		while(isInterior(x+l, y, s2) && isInterior(x+l, y, s1) ) { l++; }
		while(isInterior(x, y+w, s2) && isInterior(x, y+w, s1) ) { w++; }
		c[2] = l;
		c[3] = w; 
		result = getIdFromInt(c[0], c[1], l, w);
		
	//	cout<<"FW: s1: "<<s1<<" s2: "<<s2<<" intersect: "<<result<<" c: "<<getStr<int>(c)<<endl;
	//	getchar(); 
	} else if(isInterior(a[0]+a[2]-1, a[1]+a[3]-1, s2) ) {
		int l=0, w=0;
		int x=a[0]+a[2]-1;
		int y=a[1]+a[3]-1;
		while(isInterior(x-l, y, s2) && isInterior(x-l, y, s1) ) { l++; }
		while(isInterior(x, y-w, s2) && isInterior(x, y-w, s1) ) { w++; }
		c[0] = x-l+1;
		c[1] = y-w+1;
		c[2] = l;
		c[3] = w; 
		result= getIdFromInt(c[0], c[1], l, w); 
	//	cout<<"BW: s1: "<<s1<<" s2: "<<s2<<" intersect: "<<result<<" c: "<<endl;
	} else if(isInterior(b[0], b[1], s1) || isInterior(b[0]+b[2]-1, b[1]+b[3]-1, s1) ) {
		result = intersection(s2, s1); 
	} else {
		result = ""; 
	}
	
	return result; 
//	vector<int> bounds = getIntFromId(s1);
//	vector<int> bounds2 = getIntFromId(s2);
//	if (bounds2[2]*bounds[3] < bounds[2]*bounds[3]) {
//		return intersection(s2, s1);
//	}
//	vector<int> result(4, 0);
//	int xState=0, yState=0; // 0 not begun, 1 begun, 2 end
//	for (int i=bounds[0]; i < bounds[2]+bounds[0]; i++) {
//		for (int j=bounds[1]; j < bounds[3]+bounds[1]; j++) {
//			bool res = isInterior(i, j, s2);
//			if (xState == 2)
//				goto X_FOUND;
//			else if ( (xState == 0) && (res)) {
//				xState = 1;
//				result[0] = i;
//				result[2] = 1;
//			} else if (res) { // xState ==1 || !res || (xState==0 && !res)
//				result[2]++;
//			} else if (xState == 1) {
//				xState = 2;
//			}
//		}
//	}
//	X_FOUND: 
//	for (int j=bounds[1]; j < bounds[1]+bounds[3]; j++) {
//		for (int i=bounds[0]; i < bounds[2]+bounds[0]; i++) {
//			bool res = isInterior(i, j, s2);
//			if (yState == 2)
//				goto Y_FOUND;
//			else if ( (yState == 0) && (res)) {
//				yState = 1;
//				result[1] = j;
//				result[3] = 1;
//			} else if (res) { // xState ==1 || !res || (xState==0 && !res)
//				result[3]++;
//			} else if (yState == 1) {
//				yState = 2;
//			}
//		}
//	}
//Y_FOUND: 
//	string nid = ""; 
//	if( (result[2] > 0) && (result[3]> 0)) {
//		nid = getIdFromInt(result[0], result[1], result[2], result[3]); 
//	//	cout<<"intersects: "<<s1<<" with "<<s2<<" => "<<nid<<"\n" ; getchar(); 
//	} 
//	return nid; 
	
	//    if(isInterior(bounds2[0], bounds2[1], s1) || (bounds[1] < bounds2[1] ))
	//        return intersection(s2, s1);
	//    int l=-1, w=-1;
	//    int xn=0, yn=0;
	//    for(int i=bounds[0]; i <= min( (bounds[0]+bounds[2]), (bounds2[0]+bounds2[2]) ); i++) {
	//        if( VideoNode::isInterior(i, bounds[1], s2) )//  || isInterior(i, (bounds[1]+bounds[3]), s2) )
	//            ++l;
	//        if( l == 0) {
	//            xn = i;
	//            yn = bounds[1];
	//        }
	//        //      printf("\t l: %d xn: %d yn: %d i: %d\n", l, xn, yn, i);
	//    }
	//    if( l > 0) {
	//        for(int i=bounds[1]; i <= min( (bounds[1]+bounds[3]), (bounds2[1]+bounds2[3]) ); i++) {
	//            if( VideoNode::isInterior(xn, i, s2) )
	//                //                || VideoNode::isInterior((bounds[0]+bounds[2]), i, s2) )
	//                w++;
	//            //    if(w == 0) yn = i;
	//        }
	//    }
	//    string res = "";
	//    if( (l > 0) && (w > 0) ) {
	//        stringstream ss("");
	//        ss<<xn<<"_"<<yn<<"_"<<l<<"_"<<w;
	//        res = ss.str();
	//    }
	//    //cout<<"s1: "<<s1<<" s2: "<<s2<<" res: "<<res<<"\n";
	//    //  PRINTF("Vnode::intersections1() s1: %s s2: %s res: %s xn: %d yn: %d\n", s1.c_str(), s2.c_str(), res.c_str(), xn, yn);
	//	return res;
}

/// checks whether block is interior
bool VideoNode::isInterior(int x, int y, string sBig) {
	vector<int> nums = VideoNode::getIntFromId(sBig);
		bool ret;
	if ( (x >= nums[0]) && (y >= nums[1]) 
			&& (x < (nums[0] + nums[2])) 
			&& (y < (nums[1] + nums[3]))
			)
		ret = true;
	else
		ret = false;
	//    PRINTF("Vnode::isInterior(%d, %d, %s) = %d\n" , x, y, sBig.c_str(), ret);
	// cout<<"isInterior: ("<<x<<", "<<y<<") s: "<<sBig<<" num: "<<getStr<int>(nums)<<" = "<<ret<<endl; 
	return ret;
}

void VideoNode::finalizeCands(StatePtr _sOld) {
	sOld = _sOld;
	nextStateMp.insert(augmentsMp.begin(), augmentsMp.end());
	augmentsMp.clear();
	clear();
	getBeliefPtr()->setPr(sOld, 1.0);
	ASSERT(sOld!=0);
	lambdaR->setPr(sOld, 1.0);
	LD pAlloc = 0.0;
	int stCount=nextStateMp.size();
	LD eqProb = 1.0/((LD) stCount);
	foreach(iter, nextStateMp) {
		pAlloc += iter->second;

	}
	ASSERT(pAlloc > 0.0);
	ASSERT(sOld!=0);
	foreach(iter, nextStateMp) {
		StatePtr sNew = iter->first;
		pHat->setPr(sOld, sNew, iter->second/pAlloc);
		bJointRev->setPr(sNew, sOld, eqProb);
		lambdaPar->setPr(sOld, sNew, 1.0);
		lambdaChild->setPr(sOld, sNew, 1.0);
		bNew->setPr(sNew, eqProb);
	}

}

