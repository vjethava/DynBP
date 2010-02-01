/**
 * This file tracks the original probability distribution with time 
 */
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <pthread.h>
#include "main.h" 
#include "GraphPPM.h"
#include "PottsModel.h"  
#include "PottsGraph.h" 
using namespace std;
// DVP_steady.avg R00.avg - key files 
/// Readme file generator 
void readme(string str="") {
  FILE* file = fopen("README.Porig", "w");
  fprintf(file,
          "pX_N1xN2_sY.log has the true probability at time X for initial seed Y\n");
  fprintf(
          file,
          "bX_N1xN2_sY.log has the approximated probability at time X for initial seed Y\n");
  fprintf(file,
          "\nSR.log has the comparison between S({b_approx}) and SR({b_r})\n");
  fprintf(file,
          "\n P00.txt has the true probability for 0x0 node being 0\n");
  fprintf(
          file,
          " B00.txt has the approximate probability for 0x0 node being 0 by looking at b_approx and at region 0_0_rx_ry respectively as [time b_approx b_rgn] \n");
  fprintf(file, "\nDVP.txt has the KL-distance: [time D(b_approx||p)]\n\n");
  if (str != "")
    fprintf(file, "%s\n", str.c_str());
  fclose(file);
}
/**
 * File Name given the parameters
 */
string getFileName(int N1, int N2, int seed, LD time=0) {
  stringstream ss("");
  ss<<"p"<<time<<"_"<<N1<<"x"<<N2<<"_s"<<seed<<".log";
  //cout<<"ss: "<<ss.str()<<"\n"; 
  return ss.str();
}

vector<LD> getCondP(const int& N1, const int& N2, const int& Q, 
                    const LD& theta, const LD& delta_t,  const LD& J, const LD& H) {
  fprintf(stderr, "getCondP(N1 = %d, N2 = %d, Q = %d, theta=%g, delta_t=%g, J=%g, H=%g)\n", 
          N1, N2, Q, theta, delta_t, J, H); 
  Timer timer(1);
  int numStates= MyDef::pow(Q, N1*N2);
  PottsModel* model = new PottsModel(N1, N2, Q, J, H);
  model->genFactors();
  vector<LD> pCond(numStates*numStates, 1.0);
  for (int i=0; i < numStates; i++) {
    LD pSum = 0.0;
    cout<<"\tinit state: "<<i<<"\n";
    for (int j=0; j < numStates; j++) {
      timer.tick();
      vector<INT> sInit = model->getIthStVec(i, N1*N2);
      vector<INT> sNext = model->getIthStVec(j, N1*N2);
      // cout<<"sInit: "<<getStr<INT>(sInit)<<"\n";
      // cout<<"sNext: "<<getStr<INT>(sNext)<<"\n";
      LD pc = model->getCondP(i, j, theta*delta_t);
      pCond[i*numStates + j] = pc;
      pSum = pSum + pc;
    }
    ASSERT(pSum > 0.0);
    for (int j=0; j < numStates; j++) { // normalization
      pCond[i*numStates+j] = pCond[i*numStates+j]/pSum;
    }
  }
  timer.tick();
  cout<<"done\n";
  return pCond; 
}
int N1=3, N2=3, Q=2; LD maxTime=40.0; 
pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
vector<LD> pCond = getCondP(N1, N2, Q, 0.1, 1.0, 0.1, 0.1); 
LD delta_t = 1.0; 
/**
 * Computes the true P 
 */
vector<LD> computeTrueP(const int& N1, const int& N2, const int& Q, const int& seed, const LD& maxTime, const vector<LD>& pCond) {
 vector<LD> result; 
  Timer timer(2);
  LD ctime = 0.0;
  int numStates= MyDef::pow(Q, N1*N2);
  LD L = 1.0/((LD) numStates);
  vector<LD> pCurr(numStates, L), pNext(numStates, 0.0);
  ASSERT(pCond.size() == numStates*numStates); 
  stringstream sp00("");
  sp00<<N1<<"x"<<N2<<"_s"<<seed<<".p00";
  FILE* fileP00 = fopen(sp00.str().c_str(), "w"); //fopen(sp00.str().c_str(), "w");
  while (ctime < maxTime) {
    string fname = getFileName(N1, N2, seed, ctime);
    //  printf("time: %g => writing to %s\n", ctime, fname.c_str());
    FILE* file = fopen(fname.c_str(), "w");
    LD p00 = 0.0;
    for (int i=0; i < numStates; i++) {
      timer.tick();
      if (i%Q == 0) {
        p00 += pCurr[i];
      }
      fprintf(file, "%g\n", pCurr[i]);
    }
    fprintf(fileP00, "%g\n", p00);
    result.push_back(p00); 
    fflush(fileP00);
    fclose(file);
    // printf("computing next p ") ;
    for (int i=0; i < numStates; i++) { // orig state loop
      for (int j=0; j < numStates; j++) { // next state loop
        pNext[j] += pCurr[i]*pCond[i*numStates+j];
      }
      timer.tick();
    }
    // printf("done\n");
    LD newZ=0.0;
    for (int i=0; i < numStates; i++) {
      pCurr[i] = pNext[i];
      pNext[i] = 0;
      newZ += pCurr[i];
    }
    ctime += delta_t;
    // printf("ctime: %g newZ: %g\n", ctime, newZ);
  }
  fclose(fileP00);
  return result; 
}



void* simulate(void * i) { 
  stringstream args("");
  //int N1 = 3, N2 = 3, Q = 2;
  LD theta = 0.1, ctime=0.0;
  LD J=0.1, H=0.1;
  int *ms; ms = (int*) i;
  cout<<"ms = "<<*ms<<endl; 
  int seed= *ms;     
  // cout<<" N1: "; cin>>N1;
  // cout<<" N2: "; cin>>N2; 
  // cout<<" theta: "; cin>>theta;
  // cout<<" seed: ";  cin>>seed;
  // cout<<" H(=0.1): "; cin>>H;
  int numStates = (int) MyDef::pow(Q, N1*N2); 
  //args<<"var_params: N1: "<<N1<<" N2: "<<N2<<" Q: "<<Q<<" seed: "<<seed<<" theta: "<<theta;
  //   args<<" delta_t: "<<delta_t<<" maxTime: "<<maxTime<<"\n";
  //   args<<"fixed_params: "<<"J: "<<J<<" H: "<<H<<" fixed in computeTrueP()";
  //  readme(args.str());

  vector<LD> p00vec = computeTrueP(N1, N2, Q, seed, maxTime, pCond);
   
  PottsModel* model = new PottsModel(N1, N2, Q, J, H, seed);
  int rx=2, ry=2, cx=1, cy=1;
  model->genFactors();
  PottsGraphPPM* graph = new PottsGraphPPM(model, rx, ry, cx, cy);
  graph->genSubRegions();
  graph->initAllNodesOrigP(false);
  graph->PPM_step2(theta, delta_t);
     
  stringstream sr_ss(""); sr_ss<<"SR_"<<seed<<".txt";
  FILE* fileS = fopen(sr_ss.str().c_str(), "w");
  stringstream sb00("");
  sb00<<N1<<"x"<<N2<<"_s"<<seed<<".b00";
  vector<LD> b00vec; 
  FILE* fileB00 = fopen(sb00.str().c_str(), "w"); // fopen(sb00.str().c_str(), "w");
  while (ctime < maxTime) {
    graph->PPM_step3();
    for (int i=0; i < 3; i++) {
      graph->PPM_inner_loop(); // update bNew based on bOld
      graph->updateNodesInfo();
      string stats = graph->get_status();
      LD HR = graph->getHR();
      LD lnZR = HR-graph->SR; // SR is a class variable
      fprintf(stderr, "ITER  -lnZ_approx HR [SR  minM maxM <ParDiff> PD_max <PD_max> <Cblf> blfC_max <blfC_max>]\n");
      fprintf(stderr, "%g: Z: %g H: %g [%s]\n", ctime, lnZR, HR, stats.c_str());

      //getchar(); 
    }

    vector<LD> blfs = graph->getNormGraphBlfs();
    LD bSum = 0.0;
    LD sApprox = 0.0;

    stringstream ss("");
    ss<<"b"<<ctime<<"_"<<N1<<"x"<<N2<<"_s"<<seed<<".log";
    string fname2 = ss.str();
    FILE* file2 = fopen(fname2.c_str(), "w");
    LD b00_approx = 0.0;
    for (int i=0; i < blfs.size(); i++) {
      bSum += blfs[i];
      if (blfs[i] > 0.0)
        sApprox += -blfs[i] * log(blfs[i]);
      fprintf(file2, "%g\n", blfs[i]);
      // code for finding b00approx
      if (i%Q == 0) {
        b00_approx += blfs[i];
      }
    }
    cout<<" bSum : "<<bSum<<endl;
    fprintf(fileS, "%g %g\n", graph->SR, sApprox);
    fflush(fileS);
    // look at the first big region
    stringstream nss("");
    nss<<"0_0_"<<rx<<"_"<<ry;
    NodePtr np = graph->getNode(nss.str());
    ASSERT(np != 0);
    VideoNode* n00 = (VideoNode*) np;
    LD b00_rgn = graph->getVarMargPr(n00, 0, 0, 0);
    b00vec.push_back(b00_rgn); 
    fprintf(fileB00, "%g %g %g\n", ctime, b00_approx, b00_rgn);
    graph->updateTime(delta_t); // set bNew as actual blfs 
    ctime += delta_t;
    fclose(file2);
    // cout<<" ctime: "<<ctime<<endl; getchar(); 
  }
  fclose(fileS);
  fclose(fileB00);
  stringstream ss_dvp("");
  ss_dvp<<N1<<"x"<<N2<<"_s"<<seed<<".DVP";
  FILE* fileDvsP = fopen(ss_dvp.str().c_str(), "w"); //fopen(ss_dvp.str().c_str(), "w");
  // computing D(b||p) 
  // int numStates = (int) MyDef::pow(Q, N1*N2);
  for (LD ctime=0.0; ctime < maxTime; ctime += delta_t) {
    vector<LD> b, p;
    stringstream ss_p(""), ss_b("");
    ss_p<<"p"<<ctime<<"_"<<N1<<"x"<<N2<<"_s"<<seed<<".log";
    ss_b<<"b"<<ctime<<"_"<<N1<<"x"<<N2<<"_s"<<seed<<".log";

    ifstream pStream, bStream;
    pStream.open(ss_p.str().c_str());
    bStream.open(ss_b.str().c_str());
    //  FILE* pFile = fopen(ss_p.str().c_str(), "r");
    //  FILE* bFile = fopen(ss_b.str().c_str(), "r");
    LD dCurr = 0.0;
  
   
    for (int i=0; i < numStates; i++) {
      LD bCurr, pCurr;
      pStream>>pCurr;
      bStream>>bCurr;
      //        fscanf(pFile, "%g\n", &pCurr);
      //        fscanf(bFile, "%g\n", &bCurr); 
      //        ASSERT(pCurr > 0.0);
      if (bCurr > 0.0) {
        dCurr += (LD) ((double)bCurr)*log( ((double) bCurr)/((double) pCurr)) ;
      }
     
      pStream.close(); 
      bStream.close(); 
    }
    fprintf(fileDvsP, "%g %g\n", ctime, dCurr);
    //  getchar(); 
  }
  fclose(fileDvsP);
     
  

 
  // code for generating b/p results (assumes b & p logs are present in cwd)
  ofstream opStream;
  stringstream ss_bbyp(""); ss_bbyp<<"BbyP_s"<<seed<<".txt"; 
  opStream.open(ss_bbyp.str().c_str());
  vector<LD> vr_avg, vr_min, vr_max; 
  int NT=0; 
  for(LD ctime=0.0; ctime < maxTime; ctime+=delta_t) {
    stringstream ss_p(""), ss_b("");
    ss_p<<"p"<<ctime<<"_"<<N1<<"x"<<N2<<"_s"<<seed<<".log";
    ss_b<<"b"<<ctime<<"_"<<N1<<"x"<<N2<<"_s"<<seed<<".log";
    ifstream pStream, bStream;
    pStream.open(ss_p.str().c_str());
    bStream.open(ss_b.str().c_str());
    LD r_min=10000.0, r_max=0.0, r_avg=0.0; 
    LD diff_min=1.0, diff_max=0.0, diff_avg=0.0;
    for(int i=0; i < numStates; i++) {
      LD bCurr, pCurr;
      pStream>>pCurr;
      bStream>>bCurr;
      LD diff = abs(bCurr - pCurr); 
      diff_min= min(diff_min, diff);
      diff_max= max(diff_max, diff);
      diff_avg += diff;
      ASSERT(pCurr > 0.0); 
      LD ratio = diff/pCurr; 
      r_min = (LD) min((double) r_min, (double) ratio); 
      r_max = (LD) max((double) r_max, (double) ratio); 
      r_avg += ratio; 
  
    }
    diff_avg = diff_avg/((LD) numStates); 
    r_avg = r_avg/((LD) numStates); 
    vr_avg.push_back(r_avg); 
    vr_min.push_back(r_min); 
    vr_max.push_back(r_max); 
    opStream<<ctime<<" "<<r_min<<" "<<r_avg<<" "<<r_max<<" "<<diff_min<<" "<<diff_avg<<" "<<diff_max<<"\n"; 
    pStream.close();  
    bStream.close();
    NT++;  
  }
  
  // steady state ratios
  int NC = 10; 
  LD r_avg_ss=0.0, r_min_ss=0.0, r_max_ss=0.0;  
  for(int i=NT-1; i >= NT-NC; i--) {
    r_avg_ss += vr_avg[i];
    r_min_ss += vr_min[i];
    r_max_ss += vr_max[i];
  }
  r_avg_ss = r_avg_ss/((LD) NC); 
  r_min_ss = r_min_ss/((LD) NC); 
  r_max_ss = r_max_ss/((LD) NC); 
  opStream.close(); 

  stringstream ss_p00(""); ss_p00<<N1<<"x"<<N2<<"_s"<<seed<<".p00"; 
  stringstream ss_b00(""); ss_p00<<N1<<"x"<<N2<<"_s"<<seed<<".b00";
  ifstream pStrf, bStrf;
  pStrf.open(ss_p00.str().c_str()); 
  bStrf.open(ss_b00.str().c_str()); 
  //FILE* pFile = fopen(ss_p00.str().c_str(), "r");
  //FILE* bFile = fopen(ss_b00.str().c_str(), "r");
  printf(" b_vs_p.avg computation\n"); 
  LD r_avg=0.0, r_min=100.0, r_max=0.0;
  int ii; 
  ASSERT(b00vec.size() == p00vec.size()); 
  for(int i=0; i < b00vec.size(); i++) {
  	LD p, b_rgn, b_approx;
  	b_rgn = b00vec[i];  
  	p = p00vec[i]; 
 // 	fscanf(pFile, "%f\n", &p); 
 // 	fscanf(bFile, "%d %f %f\n", &ii, &b_approx, &b_rgn); 
  //	pStrf>>p;
 // 	bStrf>>ii>>b_approx>>b_rgn; 
  	printf("p: %g b: %g\n", p, b_rgn);
  	ASSERT(p > 0.0);  
  	LD r = abs(b_rgn-p)/p; 
  	r_avg += r;
  	r_max = max(r, r_max); 
  	r_min = min(r, r_min); 
  }
  r_avg = r_avg/b00vec.size(); 
  // fclose(pFile); fclose(bFile); 
  
  // synchronized code - append to the global files 
  pthread_mutex_lock( &mutex1 );
  FILE* dvpGlobalFile = fopen("DVP_steady.avg", "a+"); 
  fprintf(dvpGlobalFile, "%d %g %g %g\n", seed, r_avg_ss, r_min_ss, r_max_ss); 
  fclose(dvpGlobalFile);
  
  FILE* rGlobal = fopen("R00.avg", "a+"); 
  fprintf(rGlobal, "%d %g %g %g\n", seed, r_avg, r_min, r_max); 
  fclose(rGlobal);
  
  pthread_mutex_unlock( &mutex1 );

}

int main() {
  int num_threads = 1; 
  int total_threads = 1; 
  pthread_t thread[num_threads]; 
  int seed[num_threads];
  int iret[num_threads]; 
  FILE* dvpFile = fopen("DVP_steady.avg", "w"); fclose(dvpFile); 
  FILE* rGlobal = fopen("R00.avg", "w"); fclose(rGlobal); 
  fj(0, total_threads) { 
  	fi(0, num_threads) {
    	seed[i] = num_threads*j + i; 
    	iret[i] = pthread_create(&thread[i], NULL, simulate, (void*) (&seed[i])); 
  	}
  	fi(0, num_threads) {
    	pthread_join( thread[i], NULL);
  	}
  }
  return 0; 
}




