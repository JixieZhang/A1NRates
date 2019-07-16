#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iomanip>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "ACCTools.h"
#include "ACCInc.h"
//#define ACCEPTANCE_DEBUG 4

namespace ACCEPTANCE {
//put them as global variables to avoid frequently memory allocation

  //these files are too big to include in this way
  //I will read them from files later
  /*
  #include "ACCEPTANCE_SHMS_1.inc"
  #include "ACCEPTANCE_SHMS_2.inc"
  #include "ACCEPTANCE_HMS_1.inc"
  #include "ACCEPTANCE_HMS_2.inc"
  */
  double**** ACCEPTANCE_SHMS_1;
  double**** ACCEPTANCE_SHMS_2;
  double**** ACCEPTANCE_HMS_1;
  double**** ACCEPTANCE_HMS_2;
  bool IsAcceptanceLoaded = false;

  void InitAcceptance()
  {
    if (IsAcceptanceLoaded) return;
      
    int bin1=ACCEPTANCE::kYtarBinNum; 
    int bin2=ACCEPTANCE::kDeltaBinNum; 
    int bin3=ACCEPTANCE::kThetaBinNum; 
    int bin4=ACCEPTANCE::kPhiBinNum;
    double dDefVal=0;

    ACCEPTANCE_SHMS_1 =ACCEPTANCE::IniDynamicArray(bin1,bin2,bin3,bin4,dDefVal);
    ACCEPTANCE_SHMS_2 =ACCEPTANCE::IniDynamicArray(bin1,bin2,bin3,bin4,dDefVal);
    ACCEPTANCE_HMS_1 =ACCEPTANCE::IniDynamicArray(bin1,bin2,bin3,bin4,dDefVal);
    ACCEPTANCE_HMS_2 =ACCEPTANCE::IniDynamicArray(bin1,bin2,bin3,bin4,dDefVal);

    ReadAcc("ACCEPTANCE_SHMS_1.root",ACCEPTANCE_SHMS_1);
    ReadAcc("ACCEPTANCE_SHMS_2.root",ACCEPTANCE_SHMS_2);
    ReadAcc("ACCEPTANCE_HMS_1.root",ACCEPTANCE_HMS_1);
    ReadAcc("ACCEPTANCE_HMS_2.root",ACCEPTANCE_HMS_2);
    IsAcceptanceLoaded = true;
  }

  void FreeAcceptance()
  {
    if (!IsAcceptanceLoaded) return;
    
    int bin1=ACCEPTANCE::kYtarBinNum; 
    int bin2=ACCEPTANCE::kDeltaBinNum; 
    int bin3=ACCEPTANCE::kThetaBinNum; 

    ACCEPTANCE::FreeDynArray(ACCEPTANCE_SHMS_1,bin1,bin2,bin3);
    ACCEPTANCE::FreeDynArray(ACCEPTANCE_SHMS_2,bin1,bin2,bin3);
    ACCEPTANCE::FreeDynArray(ACCEPTANCE_HMS_1,bin1,bin2,bin3);
    ACCEPTANCE::FreeDynArray(ACCEPTANCE_HMS_2,bin1,bin2,bin3);
  }

  //return exclusive acceptance D(e,e'pi- p)p for the given bin, if no found, return 0;	
  //input type=3/4/5/6 are corresponsed to MM3 MM4 MM5 and all Exc
  double GetAcceptance(int run, double pYtar, double pDelta,double pTheta, double pPhi, int type,
      int &pYtarIdx, int &pDeltaIdx, int &pThetaIdx, int &pPhiIdx)
  {
    if(!IsAcceptanceLoaded) InitAcceptance();
    
    int ret=GetIndex(pYtar, pDelta, pTheta, pPhi, pYtarIdx, pDeltaIdx, pThetaIdx, pPhiIdx);
    if(ret<0) return 0.0;

    double acc=0.0;
    if(run==1) 
    {
      if(type==1)
        acc=ACCEPTANCE_HMS_1[pYtarIdx][pDeltaIdx][pThetaIdx][pPhiIdx];
      else
        acc=ACCEPTANCE_HMS_2[pYtarIdx][pDeltaIdx][pThetaIdx][pPhiIdx];            
    }
    else 
    {
      if(type==1)
        acc=ACCEPTANCE_SHMS_1[pYtarIdx][pDeltaIdx][pThetaIdx][pPhiIdx];
      else
        acc=ACCEPTANCE_SHMS_2[pYtarIdx][pDeltaIdx][pThetaIdx][pPhiIdx];
    }
    return acc;
  }

  double GetAcceptance(int run, float pYtar, float pDelta, float pTheta, float pPhi, int type,
      int &pYtarIdx, int &pDeltaIdx, int &pThetaIdx, int &pPhiIdx)
  {
    double dYtar=pYtar, dDelta=pDelta, dTheta=pTheta, dPhi=pPhi;
    return GetAcceptance(run,dYtar,dDelta,dTheta,dPhi,type,pYtarIdx,pDeltaIdx,pThetaIdx,pPhiIdx);

  }

  double GetAcceptance(int run, double pYtar, double pDelta, double pTheta, double pPhi, int type)
  {    
    int pYtarIdx,pDeltaIdx,pThetaIdx,pPhiIdx;
    return  GetAcceptance(run,pYtar,pDelta,pTheta,pPhi,type,pYtarIdx,pDeltaIdx,pThetaIdx,pPhiIdx);
  }

  double GetAcceptance(int run,float pYtar, float pDelta,float pTheta,  float pPhi, int type)
  {
    double dYtar=pYtar, dDelta=pDelta, dTheta=pTheta,  dPhi=pPhi;
    return GetAcceptance(run,dYtar,dDelta,dTheta,dPhi,type);
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  //reurn -9 if not found otherwise return 0
  int  GetYtarIndex(double pYtar, int &pYtarIdx)
  {
    pYtarIdx=abs(BinarySearch(kYtarBin,0, kYtarBinNum+1,pYtar));
    return(pYtarIdx>=kYtarBinNum)? -9:0;
  }
  //reurn -9 if not found otherwise return 0
  int  GetDeltaIndex(double pDelta, int &pDeltaIdx)
  {
    pDeltaIdx=abs(BinarySearch(kDeltaBin,0, kDeltaBinNum+1,pDelta));
    return(pDeltaIdx>=kDeltaBinNum)? -9:0;
  }

  //reurn -9 if not found otherwise return 0
  int  GetThetaStarIndex(double pTheta,  int &pThetaIdx)
  {
    pThetaIdx=abs(BinarySearch(kThetaBin,0, kThetaBinNum+1,cos(pTheta*3.14159/180.)));        
    return(pThetaIdx>=kThetaBinNum)? -9:0;
  }    

  //reurn -9 if not found otherwise return 0
  int  GetCosThetaStarIndex(double pCosThetaStar, int &pThetaIdx)
  {
    pThetaIdx=abs(BinarySearch(kThetaBin,0, kThetaBinNum+1,pCosThetaStar));        
    return(pThetaIdx>=kThetaBinNum)? -9:0;
  }

  int  GetPhiStarIndex(double pPhi, int &pPhiIdx)
  {
    if(pPhi< 0.) pPhi+=360.;
    if(pPhi>360.) pPhi-=360.;
    pPhiIdx=abs(BinarySearch(kPhiBin,0, kPhiBinNum+1,pPhi));
    return(pPhiIdx>=kPhiBinNum)? -9:0;
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  //for exclusive, use kYtarBin array 
  //reurn -9 if not found otherwise return 0
  int  GetIndex(double pYtar, double pDelta,double pTheta,  double pPhi,
      int &pYtarIdx, int &pDeltaIdx, int &pThetaIdx, int &pPhiIdx)
  {
#if defined(ACCEPTANCE_DEBUG) &&  (ACCEPTANCE_DEBUG >=4)
    std::cout<<"ACCEPTANCE::GetIndex():"<<"\t pYtar="<<pYtar<<"\t pDelta="<<pDelta<<"\t  Theta="
             <<pTheta<<"\t pPhi="<<pPhi <<std::endl;
#endif

    /*
    //this part is used to test the BinarySearch
    static const int kIntBin[]={ 0, 30, 60, 90, 120, 150, 180};
    pThetaIdx=BinarySearch(kIntBin,0,7,int(pTheta));
    pThetaIdx=BinarySearch(kIntBin,0,7,60);
    pThetaIdx=BinarySearch(kIntBin,0,7,180);
    pThetaIdx=BinarySearch(kIntBin,0,7,280);
    */
    pYtarIdx=abs(BinarySearch(kYtarBin,0, kYtarBinNum+1,pYtar));
    pDeltaIdx=abs(BinarySearch(kDeltaBin,0, kDeltaBinNum+1,pDelta));
    pThetaIdx=abs(BinarySearch(kThetaBin,0, kThetaBinNum+1,pTheta));
    pPhiIdx=abs(BinarySearch(kPhiBin,0, kPhiBinNum+1,pPhi));

    if(pYtarIdx>=kYtarBinNum || pDeltaIdx>=kDeltaBinNum || pThetaIdx>=kThetaBinNum || pPhiIdx>=kPhiBinNum)
    {
#if defined(ACCEPTANCE_DEBUG) &&  (ACCEPTANCE_DEBUG >=4)
      std::cout<<"ACCEPTANCE::GetIndex()**Bad Event:  Ytar_idx="<<pYtarIdx<<"  Delta_idx="<<pDeltaIdx
               <<"  Theta_idx="<<pThetaIdx<<"  Phi_idx="<<pPhiIdx <<std::endl;
#endif
      return -9;
    }
    return 0;  
  }

///////////////////////////////////////////////////////////////////////////////////////////////////////
/*
int BinarySearch(const int sortedArray[], int first, int last, int key) 
{
  // function:
  //   Searches sortedArray[first]..sortedArray[last-1] for key.  
  // returns: index of the matching element if it finds key, 
  //          otherwise  -index (after which it could be inserted).
  //          which satisfy sortedArray[index]< key < sortedArray[index+1].
  // parameters:
  //   sortedArray in  array of sorted (ascending) values.
  //   first, last in  lower and upper subscript bounds
  //   key         in  value to search for.
  // returns:
  //   index of key, or -insertion_position if key is not 
  //                 in the array.
  if(key<sortedArray[0] || key>sortedArray[last-1]) return -last+1;	   
  while (first <= last) {
    int mid = (first + last) / 2;  // compute mid point.
    if (key > sortedArray[mid]) 
      first = mid + 1;  // repeat search in top half.
    else if (key < sortedArray[mid]) 
      last = mid - 1; // repeat search in bottom half.
    else
      return mid;     // found it. return position /////
  }
  return -first+1;    // failed to find key
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  int BinarySearch(const double sortedArray[], int first, int last, double key) 
  {
  // function:
  //   Searches sortedArray[first]..sortedArray[last-1] for key.  
  // returns: index of the matching element if it finds key, 
  //          otherwise  -index (after which it could be inserted).
  //          which satisfy sortedArray[index]< key < sortedArray[index+1].
  // parameters:
  //   sortedArray in  array of sorted (ascending) values.
  //   first, last in  lower and upper subscript bounds
  //   key         in  value to search for.
  // returns:
  //   index of key, or -insertion_position if key is not 
  //                 in the array. 
  if(key<sortedArray[0] || key>sortedArray[last-1]) return -last+1;	 	   
  while (first <= last) 
  {
    int mid = (first + last) / 2;  // compute mid point.
    if (key > sortedArray[mid]) 
      first = mid + 1;  // repeat search in top half.
    else if (key < sortedArray[mid]) 
      last = mid - 1;  // repeat search in bottom half.
    else
      return mid;      // found it. return position /////
  }
  return -first+1;     // failed to find key
}
*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////

  void ReadAcc(const char *filename, double**** AccArray)
  {
    //std::cout<<" Loading acceptance file '"<<filename<<"' ... "<<std::endl;
    
    TFile *file = new TFile(filename);
    TTree *tree = (TTree*)gDirectory->Get("acc");

    //Declaration of leaves types
    Int_t           iyt;
    Int_t           ide;
    Int_t           ith;
    Int_t           iph;
    Float_t         Ytar;
    Float_t         Delta;
    Float_t         Theta;
    Float_t         Phi;
    Float_t         acc;
    Int_t           n;
    Int_t           n_true;

    // Set branch addresses.
    tree->SetBranchAddress("iyt",&iyt);
    tree->SetBranchAddress("ide",&ide);
    tree->SetBranchAddress("ith",&ith);
    tree->SetBranchAddress("iph",&iph);
    tree->SetBranchAddress("Ytar",&Ytar);
    tree->SetBranchAddress("Delta",&Delta);
    tree->SetBranchAddress("Theta",&Theta);
    tree->SetBranchAddress("Phi",&Phi);
    tree->SetBranchAddress("acc",&acc);
    tree->SetBranchAddress("n",&n);
    tree->SetBranchAddress("n_true",&n_true);

    Long64_t nentries = tree->GetEntries();

    Long64_t nbytes = 0;
    for (Long64_t i=0; i<nentries;i++) {
      nbytes += tree->GetEntry(i);
      AccArray[iyt][ide][ith][iph]=acc;
    }
    file->Close();
    file->Delete();
    //std::cout<<" Acceptance file '"<<filename<<"' loaded, "<<nentries<<" records stored"<<std::endl;
  }

}; //end of namespace ACCEPTANCE

