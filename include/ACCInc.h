
#ifndef _ACCEPTANCE_EXC_NAMESPACE_
#define _ACCEPTANCE_EXC_NAMESPACE_

#include "ACCIncBin.h"
namespace ACCEPTANCE {
    
  void InitAcceptance();
  void FreeAcceptance();

  /////////////////////////////////////////////////////////////////////////////
  //reurn -9 if not found otherwise return 0
  int  GetYtarIndex(double pYtar, int &pYtarIdx);
  int  GetDeltaIndex(double pDelta, int &pDeltaIdx);
  int  GetThetaIndex(double pTheta, int &pThetaIdx);
  int  GetPhiIndex(double pPhi, int &pPhiIdx);

  //reurn -9 if not found otherwise return 0
  int GetIndex(double pYtar, double pDelta, double pTheta, double pPhi,
      int &pYtarIdx, int &pDeltaIdx, int &pThetaIdx, int &pPhiIdx);

  //return the acceptance for the given bins, if not found, return 0;
  double GetAcceptance(int run,double pYtar, double pDelta, double pTheta, double pPhi, int type,
         int &pWIdx,int &pDeltaIdx,int &pThetaIdx,int &pPhiIdx);
  double GetAcceptance(int run, float pYtar, float pDelta, float pTheta, float pPhi, int type,
         int &pYtarIdx, int &pDeltaIdx, int &pThetaIdx, int &pPhiIdx);
  double GetAcceptance(int run,double pYtar,double pDelta,double pTheta, double pPhi,int type);
  double GetAcceptance(int run,float pYtar, float pDelta,float pTheta,  float pPhi, int type);
     
  /////////////////////////////////////////////////////////////////////////////
  void ReadAcc(const char *filename, double**** AccArray);
};

#endif //_ACCEPTANCE_EXC_NAMESPACE_

