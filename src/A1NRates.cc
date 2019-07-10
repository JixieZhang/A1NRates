#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <time.h>

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TVector3.h>
#include <TRandom.h>

#include "XSModel.hh"          //XSModel return DiffXS in ub/MeV-Sr for inelas or ub for elas.
#include "G2PRand.hh"
#include "HRSTransform_TCSNHCS.hh"
#include "ACCInc.h"

using namespace std;

//uncomment the next line if you want to compare these 2 function: the old one is binned in Theta, the new one is binnedin cosTheta 
//#define CMP_GETINTEXS

//DEBUG level will be used to control what message to print:
// >=1: will print pressure curve rates for H2, N2 and 3He target
// >=2: will print rates for each x bin for 3He target
// >=3: add VZ bin rates
//      if CMP_GETINTEXS is defined, will also print the comparison between old and new GetInteXS()
// >=4: print a lot of debug information in GetInteXS()     
// >=5: print more debug information in GetInteXS()      
// >=6: print even more debug information in GetInteXS() and GetXS()     
#define DEBUG 2

static const double deg = acos(0.0)/90.0;
double GetXS(float Ei, float Ef, float Theta, float Tb, float Ta, int ElasOnly=0);
double GetXS_GE180(float Ei, float Ef, float Theta, float Tb, float Ta, int ElasOnly=0);

////////////////////////////////////////////////////////////////////////////
//return in unit of 10^30, corresponsed to microbarn, number of nuclei per second per cm2
double GetLumi10pow30(double current_na, int nuclei_per_molecule, double mol_mass_gpermol, double thickXdens_gpercm2)
{
  double I=current_na * 1.0e-9;
  double T=thickXdens_gpercm2 ;

  const double N_A = 6.022e+23; //Avagadro's number
  const double q_e = 1.602e-19; //electron charge

  double Lumi=0;
  Lumi = nuclei_per_molecule * ( I/q_e *N_A* T / mol_mass_gpermol ) /1.0e30;

  //L=3.759*1.0e33 * current_na*thickXdens_gpercm2/mol_mass_gpermol;
  if(DEBUG>=6) cout<<"Luminosity="<<Lumi<<" X 1.0E+30 #/s/cm2"<<endl;
  return Lumi;
}

////////////////////////////////////////////////////////////////////////////
//return elastic kinimatics for given beam energy and theta
double GetElasKin(double beam_gev, double theta_e_rad, double M_gev, 
                  double &Ef, double &ptot_p, double & theta_p_rad)
{
  double Ei=beam_gev, t=theta_e_rad, M=M_gev;
  Ef = Ei / (1+2*Ei/M*pow(sin(t/2),2.0));

  double Pperp_p,Pz_p;
  Pperp_p = -Ef * sin(t);
  Pz_p = Ei - Ef * cos(t);
  ptot_p = sqrt(Pperp_p*Pperp_p+Pz_p*Pz_p);
  theta_p_rad = asin(Pperp_p/ptot_p);

  if(DEBUG>=6) {
    cout<< "Ei="<<Ei<< " Theta_e="<<t/deg<<" M="<<M<<"  ==>  Ef=" <<Ef
        <<"  P_p=" <<ptot_p<<"   Theta_p=" <<theta_p_rad/deg<<endl;
  }

  return Ef; 
}

double GetElasEprime(double beam_gev, double theta_e_rad, double M_gev)
{
  double Ei=beam_gev, t=theta_e_rad, M=M_gev;
  double Ef = Ei / (1+2*Ei/M*pow(sin(t/2),2.0));
  return Ef; 
}

//xs wrapper
double GetXS(int Z, int N, float Ei, float Ef, float Theta, float Tb, float Ta, int ElasOnly)
{
  if(Z<0) return GetXS_GE180(Ei,Ef,Theta,Tb,Ta,ElasOnly);
  
  double pXs=0.0;
  if(ElasOnly) {
    pXs = ElasModel::GetXS(Z, N, (double)Ei, (double)Theta);
    if(isnanf(pXs)) {
      pXs=0.0;
#ifdef DEBUG 
      if(DEBUG>=4) cout<<" isnan  from ElasModel::GetXS(Z="<<Z<<", N="<<N<<", Ei="<<Ei<<", Theta="<<Theta<<")\n";
#endif
    }
    pXs *= 1.0;   //already in ub/Sr
  } else {
    //note that PBosted::GetXS() will not work for Q2>11, I give up these points
    pXs = PBosted::GetXS(Z, N, (double)Ei, (double)Ef, (double)Theta, (double)Tb, (double)Ta);
    if(isnanf(pXs)) {
      pXs=0.0;
#ifdef DEBUG 
      if(DEBUG>=4) cout<<" isnan  from PBosted::GetXS(Z="<<Z<<", N="<<N<<", Ei="<<Ei<<", Ef="<<Ef<<", Theta="<<Theta<<", Tb="<<Tb<<", Ta="<<Ta<<")\n";
      //char cc[100];cout<<"\nPress any key to continue ...";cin>>cc;
#endif
    }
    pXs *= 1000.0;   //turn from ub/MeV/Sr into ub/GeV/Sr
  }
  return pXs;
}


//XS for GE180 compound
double GetXS_GE180(float Ei, float Ef, float Theta, float Tb, float Ta, int ElasOnly)
{
  //http://galileo.phys.virginia.edu/research/groups/spinphysics/glass_properties.html
  //GE180 glass composition:
  //SiO2  60.3%
  //BaO   18.2%
  //Al2O3 14.3%
  //CaO   6.5%
  //SrO   0.25%
  // Z/A = 0.4829
  // The total of the above is 99.55%, what is the rest? I treat it as H
  //H    0.45%
  
  double MolMass_SiO2=(28.0856*1+15.9994*2);  //in g/mol   
  double MolMass_BaO =(137.327*1+15.9994*1);  //in g/mol
  double MolMass_Al2O3=(26.982*2+15.9994*3);  //in g/mol
  double MolMass_CaO =(40.078*1+15.9994*1);   //in g/mol
  double MolMass_SrO =(87.62*1+15.9994*1);    //in g/mol
  double MolMass_pr=1.00794;                  //in g/mol
    
  //const char*  kGE180_MName[]={"SiO2","BaO","Al2O3","CaO","SrO","H"};
  const double kGE180_MFraction[]={0.603,0.182,0.143,0.065,0.0025,0.0045};
  const double kGE180_MMolMass[]={MolMass_SiO2,MolMass_BaO,MolMass_Al2O3,MolMass_CaO,MolMass_SrO,MolMass_pr};
        
  double pMolMass_aver_ge180_rev = 0;
  for(int i=0;i<6;i++) {
    pMolMass_aver_ge180_rev += kGE180_MFraction[i]/kGE180_MMolMass[i];
  }
  double pMolMass_aver_ge180 = 1.0 / pMolMass_aver_ge180_rev;
  int pZ_aver_ge180 = lround(pMolMass_aver_ge180 * 0.4829);
  int pN_aver_ge180 = lround(pMolMass_aver_ge180-pZ_aver_ge180);  //round up to nearest int
    

  //comput xs using average Z and N
  double pXs=0.0;
  pXs = GetXS(pZ_aver_ge180, pN_aver_ge180, Ei, Ef, Theta, Tb, Ta, ElasOnly);
    
  //comput xs for each element
  double pXs_O  = GetXS( 8,  8, Ei, Ef, Theta, Tb, Ta, ElasOnly);
  double pXs_Si = GetXS(14, 14, Ei, Ef, Theta, Tb, Ta, ElasOnly);
  double pXs_Ba = GetXS(56, 82, Ei, Ef, Theta, Tb, Ta, ElasOnly);
  double pXs_Al = GetXS(13, 14, Ei, Ef, Theta, Tb, Ta, ElasOnly);
  double pXs_Ca = GetXS(20, 20, Ei, Ef, Theta, Tb, Ta, ElasOnly);
  double pXs_Sr = GetXS(38, 50, Ei, Ef, Theta, Tb, Ta, ElasOnly);
  double pXs_H  = GetXS( 1,  0, Ei, Ef, Theta, Tb, Ta, ElasOnly);

  const double kGE180_MXs[] = {pXs_Si+2*pXs_O, pXs_Ba+pXs_O, 2*pXs_Al+3*pXs_O,pXs_Ca+pXs_O, pXs_Sr+pXs_O, pXs_H};
    
  double pXs_ge180=0.0;
  for(int i=0;i<6;i++) {
    pXs_ge180 += pMolMass_aver_ge180*kGE180_MFraction[i]/kGE180_MMolMass[i]*kGE180_MXs[i];
  }

#ifdef DEBUG 
  if(DEBUG>=5) {
    cout<<" using average Z and N, XS_GE180 = "<<pXs
        <<"  sum up all elements, XS_GE180 = "<<pXs_ge180<<endl;
  }
#endif
  return pXs_ge180;
}


////////////////////////////////////////////////////////////////////////////
//old version is binned in Theta-Phi-P, new version is binned in CosTheta-Phi-P
//get integrated xs for given element at given vertex z,
//applying x, q2 and w cuts if their upper limits are larger than zero
//ElasOnly=-1: pure inelastic for full acceptance
//ElasOnly=0:  inelastic + elastic for full acceptance
//ElasOnly=1:  pure elastic for full acceptance
//ElasOnly=-30: pure inelastic for 2-SC-Bar acceptance
//ElasOnly=30: inelastic + elastic for 2-SC-Bar acceptance
//ElasOnly=31: pure elastic for 2-SC-Bar acceptance
double GetInteXS(double pBeamE, double pAngle, double pMomentum, double Z, int N, double VZ, string Name, int ElasOnly=0,
                 double pXmin=-1.0, double pXmax=-1.0, double pQ2min=-1.0, double pQ2max=-1.0,double pWmin=-1.0, double pWmax=-1.0)
{
  int pre = cout.precision();   //get the original precision of cout
  cout.precision(4);
  //for HMS
  int run = 1;
  double pTheta_tg_max=0.058;   //to better match method2
  double pTheta_tg_min=-0.054;
  double pPhi_tg_max=0.031;    
  double pPhi_tg_min=-0.031; 
  double pEprime_max=pMomentum*1.13;
  double pEprime_min=pMomentum*0.89;
  if(Name=="SHMS") {
    run = 2;
    pTheta_tg_max=0.046;
    pTheta_tg_min=-0.046;    //to better match method2
    pPhi_tg_max=0.025;
    pPhi_tg_min=-0.028;
    pEprime_max=pMomentum*1.27;
    pEprime_min=pMomentum*0.82;
  }
  if(pEprime_max>pBeamE) pEprime_max=pBeamE;
    
  const double AMU = 0.9314941;
  double pMtg=(fabs(Z)+fabs(N))*AMU;
  double pP_elas=0.0;

  double pTheta_tg,pPhi_tg=0;
  double pTheta=-999.,pCosTh=-999.,pPhi=-999.,pEprime=-999.;
  double pXs_elas=0.0,pXs=0.0,pInteXs=0.0;

  double pTheta_min = pAngle - 20*deg;
  double pTheta_max = pAngle + 20*deg;
  if(pTheta_min < 5*deg) pTheta_min = 5*deg;
  double pCosTh_min = cos(pTheta_max);
  double pCosTh_max = cos(pTheta_min);  
  
  double pPhi_min = -80 *deg ;
  double pPhi_max =  80 *deg ;

  //double dCosTh = 0.0005;  //xs changes rapidly in theta, therefore put tiny bin width here
  //double dPhi = 0.002;
  //double dEprime = 0.002;
  //If Change dPhi from 0.002 to 0.01, geometry acceptance change is <0.3% for all targets
  //If change dCosTh from 0.0005 to 0.001, 3He geometry acceptance will change 4%
  //but geometry acceptance  for GE180(glass) at vz=-20cm will change ~7%
  //If change dEprime from 0.002 to 0.005, geometry acceptance change is <0.3% for all targets
  //except geometry acceptance for GE180(glass) at vz=-20cm will change ~7%
  double dCosTh = 0.0005;  //xs changes rapidly in theta, therefore put tiny bin width here
  double dPhi = 0.01;
  double dEprime = 0.002;
  
  double x_tg,y_tg,z_tg=0;

  pCosTh = pCosTh_min - 0.5*dCosTh;
  while (pCosTh < pCosTh_max) {
    
    pCosTh += dCosTh;
    pTheta = acos(pCosTh);
    //get elastic Eprime 
    pP_elas = GetElasEprime(pBeamE, pTheta, pMtg);
    double dOmega = dCosTh*dPhi;

    //reset pEprime_max for each theta since pP_elas is the maximum momentum of a physics event
    double thisEprime_max = pEprime_max;
    if(thisEprime_max>pP_elas) thisEprime_max=pP_elas;
    if(DEBUG>=4) {
      cout<<" pTheta="<<pTheta/deg<<"(deg),  pEprime_min="<<pEprime_min<<"  pEprime_max="<<pEprime_max<<"  thisEprime_max="<<thisEprime_max<<endl;
    }

    pPhi = pPhi_min - 0.5*dPhi;
    while (pPhi < pPhi_max) {
      pPhi += dPhi;

      //convert Lab frame to Tranportation frame
      //void  P_HCS2TCS(double Theta_hall, double Phi_hall, double EndPlaneTheta_hall, double &Theta_tr, double &Phi_tr)
      Transform::P_HCS2TCS(pTheta, pPhi, pAngle, pTheta_tg, pPhi_tg);
      //Jixie: do not need this cut any more, but I found that method 1 is 10% more than method 2
      //I am using these 2 cuts to make them match ( I could be wrong, since method 2 might underestimate rates by itself ...)
      //conclusion:
      //1) Theta_tg cut will reduce up to 15% geometry acceptance for HMS, up to 7% geometry acceptance for SHMS
      //2) Phi_tg cut will reduce ~3.5% geometry acceptance.
      //3) I would like to cut only onto Theta_tg, since Phi_tg has strong VZ dependence.  If I cut Phi_tg, GE180 calculation
      //   will not be reliable.
      if(pTheta_tg>pTheta_tg_max || pTheta_tg<pTheta_tg_min)  continue;
      //if(pPhi_tg>pPhi_tg_max || pPhi_tg<pPhi_tg_min)  continue;  
      
      //project to target plane     
      //void Project(double x,double y,double z,double z_drift,double theta,double phi,double &x_out, double &y_out, double &z_out)
      x_tg=y_tg=z_tg=0;
      
      //void  X_HCS2TCS(double x, double y, double z,double EndPlaneTheta_hall, double &x_tr, double &y_tr, double &z_tr)
      double x_tr=0, y_tr=0, z_tr=0;
      Transform::X_HCS2TCS(0,0,VZ,pAngle,x_tr,y_tr,z_tr);
      Transform::Project(x_tr, y_tr, z_tr, -z_tr,pTheta_tg,pPhi_tg,x_tg,y_tg,z_tg);
      if(DEBUG>=4) {
        cout<<" x_tr="<<x_tr<<"  y_tr="<<y_tr<<"  z_tr="<<z_tr<<"  theta_tr="<<pTheta_tg<<"  phi_tr="<<pPhi_tg
            <<" ==> x_tg="<<x_tg<<"  y_tg="<<y_tg<<"  z_tg="<<z_tg<<endl;
      }

      pEprime = pEprime_min - 0.5*dEprime;
      double pEprime_up = pEprime + 0.5*dEprime;
      while (pEprime < thisEprime_max) {
        //check if this is the last bin
        double deltaEprime = 0.0;
        double pEprimeLeft = pEprime_max - pEprime_up;   //how much Eprime has not been integrated yet
        if(pEprimeLeft <= 0.0) break;
        if(pEprimeLeft > dEprime) {
          pEprime += dEprime;
          pEprime_up += dEprime;
          deltaEprime = dEprime;
        } else {
          pEprime = pEprime_up + pEprimeLeft/2.0;
          pEprime_up += pEprimeLeft;
          deltaEprime = pEprimeLeft;
        }
        
        ////////////////////////////////////////////////////////////////////////
        //now apply acc cut
        //double ACCEPTANCE::GetAcceptance(int run,double pYtar,double pDelta,double pTheta, double pPhi,int type);
        int pType = (ElasOnly==30 || ElasOnly==31)?2:1;   //if (ElasOnly==3), only turn on 2 SC bars
        double pDelta = (pEprime - pMomentum) / pMomentum * 100.;
        double pAcc = ACCEPTANCE::GetAcceptance(run,y_tg,pDelta,pTheta_tg,pPhi_tg,pType);
        if(DEBUG>=4) cout<<" y_tg="<<y_tg<<" theta_tg="<<pTheta_tg<<"  phi_tg="<<pPhi_tg<<"  delta="<<pDelta<<"  ==> pAcc="<<pAcc<<endl;
        if(pAcc<0.1) continue;
        
        ////////////////////////////////////////////////////////////////////////
        double pQ2=2.0*pBeamE*pEprime*(1.0-cos(pTheta));
        double pXbj=pQ2/(2.0*0.9383*(pBeamE-pEprime));
        double M_p=0.9383;
        double nu=pBeamE-pEprime;
        double pW=sqrt(pow(M_p,2.0)+2*M_p*nu-pQ2);
        //applying X cuts
        if(pXmax>0.0 && pXmax>pXmin) {
          if(pXbj<pXmin || pXbj>pXmax) continue;
        }
        //applying Q2 cuts
        if(pQ2max>0.0 && pQ2max>pQ2min) {
          if(pQ2<pQ2min || pQ2>pQ2max) continue;
        }
        //applying W cuts: for Delta transverse Asymmetry, need to apply W cut for 3He Delta resonance peak        
        if(pWmax>0.0 && pWmax>pWmin) {
          if(pW<pWmin || pW>pWmax) continue;//continue for jump the rest of the commands
        }
        
        if(DEBUG>=4) cout<<" pTheta="<<pTheta/deg<<"  pPhi="<<pPhi/deg<<"  pEprime="<<pEprime<<"  deltaEprime="<<deltaEprime<<"  pP_elas="<<pP_elas<<endl;
        
        if(ElasOnly!=1 && ElasOnly!=31) {
          pXs = 0.0;
          if(pQ2 < 11.0) pXs = GetXS(Z, N, pBeamE, pEprime, pTheta, 0.000, 0.000, 0);
          if(DEBUG>=4) cout<<" Xbj = "<< pXbj<<"  Q2 = "<<pQ2<<"  PBosted::GetXS() = "<<pXs*1.0E3<<" (nb/GeV/Sr)"<<endl;
          pInteXs += pXs*deltaEprime*dOmega*pAcc;
        }
        
        //check if elastic events are accepted or not, add only once per dOmega
        //please note that elas xs is for per nucleus
        if(fabs(pP_elas-pEprime)<0.51*deltaEprime && ElasOnly>=0) {
          pXs_elas = GetXS(Z, N, pBeamE, pEprime, pTheta, 0.000, 0.000, 1);
          if(DEBUG>=4) cout<<" pP_elas="<<pP_elas<<" deltaEprime="<<deltaEprime<<",  ElasModel::GetXS() = "<<pXs_elas*1.0E3<<" (nb/Sr)"<<endl;
          pInteXs += pXs_elas*dOmega*pAcc;
          if(DEBUG>=7) {char cc[100];cout<<"\nPress any key to continue ...";cin>>cc;}
        } 
        
      }
    }    
  }
  
  cout.precision(pre);    //recover cout default precision
  return pInteXs;
}

////////////////////////////////////////////////////////////////////////////
//old version is binned in Theta-Phi-P, new version is binned in CosTheta-Phi-P
//get integrated xs for given element at given vertex z,
//applying x and q2 cuts if their upper limits are larger than zero
//ElasOnly=-1: pure inelastic for full acceptance
//ElasOnly=0:  inelastic + elastic for full acceptance
//ElasOnly=1:  pure elastic for full acceptance
//ElasOnly=-30: pure inelastic for 2-SC-Bar acceptance
//ElasOnly=30: inelastic + elastic for 2-SC-Bar acceptance
//ElasOnly=31: pure elastic for 2-SC-Bar acceptance
double GetInteXS_old(double pBeamE, double pAngle, double pMomentum, double Z, int N, double VZ, string Name, int ElasOnly=0,
                     double pXmin=-1.0, double pXmax=-1.0, double pQ2min=-1.0, double pQ2max=-1.0,double pWmin=-1.0, double pWmax=-1.0)
{
  int pre = cout.precision();   //get the original precision of cout
  cout.precision(4);
  //for HMS
  int run = 1;
  double pTheta_tg_max=0.058;    //to better match method2
  double pTheta_tg_min=-0.054;
  double pPhi_tg_max=0.031;
  double pPhi_tg_min=-0.031; 
  double pEprime_max=pMomentum*1.13;
  double pEprime_min=pMomentum*0.89;
  if(Name=="SHMS") {
    run = 2;
    pTheta_tg_max=0.046;
    pTheta_tg_min=-0.046;
    pPhi_tg_max=0.025;
    pPhi_tg_min=-0.028;
    pEprime_max=pMomentum*1.27;
    pEprime_min=pMomentum*0.82;
  }
  if(pEprime_max>pBeamE) pEprime_max=pBeamE;
    
  const double AMU = 0.9314941;
  double pMtg=(fabs(Z)+fabs(N))*AMU;
  double pP_elas=0.0;

  double pTheta_tg,pPhi_tg=0;
  double pTheta=-999.,pPhi=-999.,pEprime=-999.;
  double pXs_elas=0.0,pXs=0.0,pInteXs=0.0;

  double pTheta_min = pAngle - 20*deg;
  double pTheta_max = pAngle + 20*deg;
  if(pTheta_min < 5*deg) pTheta_min = 5*deg;
  
  double pPhi_min = -80 *deg ;
  double pPhi_max =  80 *deg ;

  double dTheta = 0.2*deg;  //xs changes rapidly in theta, therefore put tiny bin width here
  double dPhi = 0.5*deg;
  double dEprime = 0.005;
  //double dTheta = 0.5*deg;
  //double dPhi = 1.0*deg;
  //double dEprime = 0.01;
  
  double x_tg,y_tg,z_tg=0;

  pTheta = pTheta_min - 0.5*dTheta;
  while (pTheta < pTheta_max) {
    
    pTheta += dTheta;
    //get elastic Eprime 
    pP_elas = GetElasEprime(pBeamE, pTheta, pMtg);
    double dOmega = (cos(pTheta-0.5*dTheta)-cos(pTheta+0.5*dTheta))*dPhi;

    //reset pEprime_max for each theta since pP_elas is the maximum momentum of a physics event
    double thisEprime_max = pEprime_max;
    if(thisEprime_max>pP_elas) thisEprime_max=pP_elas;
    if(DEBUG>=4) {
      cout<<" pTheta="<<pTheta/deg<<"(deg),  pEprime_min="<<pEprime_min<<"  pEprime_max="<<pEprime_max<<"  thisEprime_max="<<thisEprime_max<<endl;
    }

    pPhi = pPhi_min - 0.5*dPhi;
    while (pPhi < pPhi_max) {
      pPhi += dPhi;

      //convert Lab frame to Tranportation frame
      //void  P_HCS2TCS(double Theta_hall, double Phi_hall, double EndPlaneTheta_hall, double &Theta_tr, double &Phi_tr)
      Transform::P_HCS2TCS(pTheta, pPhi, pAngle, pTheta_tg, pPhi_tg);
      //Jixie: do not need this cut any more
      if(pTheta_tg>pTheta_tg_max || pTheta_tg<pTheta_tg_min)  continue;
      //if(pPhi_tg>pPhi_tg_max || pPhi_tg<pPhi_tg_min)  continue;
      
      //project to target plane     
      //void Project(double x,double y,double z,double z_drift,double theta,double phi,double &x_out, double &y_out, double &z_out)
      x_tg=y_tg=z_tg=0;
      
      //void  X_HCS2TCS(double x, double y, double z,double EndPlaneTheta_hall, double &x_tr, double &y_tr, double &z_tr)
      double x_tr=0, y_tr=0, z_tr=0;
      Transform::X_HCS2TCS(0,0,VZ,pAngle,x_tr,y_tr,z_tr);
      Transform::Project(x_tr, y_tr, z_tr, -z_tr,pTheta_tg,pPhi_tg,x_tg,y_tg,z_tg);
      if(DEBUG>=4) {
        cout<<" x_tr="<<x_tr<<"  y_tr="<<y_tr<<"  z_tr="<<z_tr<<"  theta_tr="<<pTheta_tg<<"  phi_tr="<<pPhi_tg
            <<" ==> x_tg="<<x_tg<<"  y_tg="<<y_tg<<"  z_tg="<<z_tg<<endl;
      }

      pEprime = pEprime_min - 0.5*dEprime;
      double pEprime_up = pEprime + 0.5*dEprime;
      while (pEprime < thisEprime_max) {
        //check if this is the last bin
        double deltaEprime = 0.0;
        double pEprimeLeft = pEprime_max - pEprime_up;   //how much Eprime has not been integrated yet
        if(pEprimeLeft <= 0.0) break;
        if(pEprimeLeft > dEprime) {
          pEprime += dEprime;
          pEprime_up += dEprime;
          deltaEprime = dEprime;
        } else {
          pEprime = pEprime_up + pEprimeLeft/2.0;
          pEprime_up += pEprimeLeft;
          deltaEprime = pEprimeLeft;
        }
        
        ////////////////////////////////////////////////////////////////////////
        //now apply acc cut
        //double ACCEPTANCE::GetAcceptance(int run,double pYtar,double pDelta,double pTheta, double pPhi,int type);
        int pType = (ElasOnly==30 || ElasOnly==31)?2:1;   //if (ElasOnly==3), only turn on 2 SC bars
        double pDelta = (pEprime - pMomentum) / pMomentum * 100.;
        double pAcc = ACCEPTANCE::GetAcceptance(run,y_tg,pDelta,pTheta_tg,pPhi_tg,pType);
        if(DEBUG>=4) cout<<" y_tg="<<y_tg<<" theta_tg="<<pTheta_tg<<"  phi_tg="<<pPhi_tg<<"  delta="<<pDelta<<"  ==> pAcc="<<pAcc<<endl;
        if(pAcc<0.1) continue;
        
        ////////////////////////////////////////////////////////////////////////
        double pQ2=2.0*pBeamE*pEprime*(1.0-cos(pTheta));
        double pXbj=pQ2/(2.0*0.9383*(pBeamE-pEprime));
        double M_p=0.9383;
        double nu=pBeamE-pEprime;
        double pW=sqrt(pow(M_p,2.0)+2*M_p*nu-pQ2);
        //applying X cuts
        if(pXmax>0.0 && pXmax>pXmin) {
          if(pXbj<pXmin || pXbj>pXmax) continue;
        }
        //applying Q2 cuts
        if(pQ2max>0.0 && pQ2max>pQ2min) {
          if(pQ2<pQ2min || pQ2>pQ2max) continue;
        }
        //applying W cuts: for Delta transverse Asymmetry, need to apply W cut for 3He Delta resonance peak        
        if(pWmax>0.0 && pWmax>pWmin) {
          if(pW<pWmin || pW>pWmax) continue;//continue for jump the rest of the commands
        }
        
        if(DEBUG>=4) cout<<" pTheta="<<pTheta/deg<<"  pPhi="<<pPhi/deg<<"  pEprime="<<pEprime<<"  deltaEprime="<<deltaEprime<<"  pP_elas="<<pP_elas<<endl;
        
        if(ElasOnly!=1 && ElasOnly!=31) {
          pXs = 0.0;
          if(pQ2 < 11.0) pXs = GetXS(Z, N, pBeamE, pEprime, pTheta, 0.000, 0.000, 0);
          if(DEBUG>=4) cout<<" Xbj = "<< pXbj<<"  Q2 = "<<pQ2<<"  PBosted::GetXS() = "<<pXs*1.0E3<<" (nb/GeV/Sr)"<<endl;
          pInteXs += pXs*deltaEprime*dOmega*pAcc;
        }
        
        //check if elastic events are accepted or not, add only once per dOmega
        //please note that elas xs is for per nucleus
        if(fabs(pP_elas-pEprime)<0.51*deltaEprime && ElasOnly>=0) {
          pXs_elas = GetXS(Z, N, pBeamE, pEprime, pTheta, 0.000, 0.000, 1);
          if(DEBUG>=4) cout<<" pP_elas="<<pP_elas<<" deltaEprime="<<deltaEprime<<",  ElasModel::GetXS() = "<<pXs_elas*1.0E3<<" (nb/Sr)"<<endl;
          pInteXs += pXs_elas*dOmega*pAcc;
          if(DEBUG>=7) {char cc[100];cout<<"\nPress any key to continue ...";cin>>cc;}
        } 
        
      }
    }    
  }
  
  cout.precision(pre);    //recover cout default precision
  return pInteXs;
}

////////////////////////////////////////////////////////////////////////////
//Beam Current in uA, All energies are in GeV unit. All angles are in radian unit.
double GetRate(double pBeamCurrent, double pBeamE, double pDetectorAngle, double pDetectorMomentum, string pDetectorName, int pElasOnly=0)
{
  if(pElasOnly==-1)  cout<<"\n===================Full acceptance:  Pure Inelastic =======================\n";
  if(pElasOnly==0)   cout<<"\n===================Full acceptance:  Inelastic + Elastic ==================\n";
  if(pElasOnly==1)   cout<<"\n===================Full acceptance:  Pure elastic =========================\n";
  if(pElasOnly==2)   cout<<"\n===================Full acceptance:  With Delta (1.1<W<1.35) Cut ==========\n";
  if(pElasOnly==4)   cout<<"\n===================Full acceptance:  With DIS (2.0<W<100.0) Cut ===========\n";
  if(pElasOnly==-30) cout<<"\n===================Only 2 SC bars:  Pure Inelastic ========================\n";
  if(pElasOnly==30)  cout<<"\n===================Only 2 SC bars:  Inelastic + Elastic ===================\n";
  if(pElasOnly==31)  cout<<"\n===================Only 2 SC bars:  Pure elastic ==========================\n";
  int pre = cout.precision();   //get the original precision of cout
  //define material 

  //const double amu = 0.93149403;//molar mass constant in GeV/c2

  double MolMass_pr=1.00794;    //in g/mol
  double MolMass_he3=3.01603;   //in g/mol
  double MolMass_c12=12.01078;  //in g/mol
  double MolMass_n14=14.0067;   //in g/mol
  double MolMass_ge180 = 54.7251;
  
  double Z_ge180 = MolMass_ge180 * 0.4829;
  int N_ge180 = lround(MolMass_ge180-Z_ge180);  //round up to nearest int
  Z_ge180 = MolMass_ge180 - N_ge180;
  
  if(DEBUG>=3) cout<<"  Z_ge180="<<Z_ge180<<"  N_ge180="<<N_ge180<<"    MolMass_ge180="<<MolMass_ge180<<endl;
  
  //double DH2=8.349E-4; // g/cm3; H2 gas density at 21.1°C (70°F) 10atm; this is correct
  //double DN2=0.0116;   // g/cm3; N2 gas density at 21.1°C (70°F) 10atm;
  double Dc12=2.267;     // g/cm3;
  double Dhe3=1.6147E-3; // g/cm3;  //Density at 12 amg, 
  double Dge180 = 2.77;  // g/cm3;

  //comput density using ideal gas law
  double R_const=82.057; //calculate gas density from ideal gas law, R=82.057 cm^3*atm*k^(-1)*mol^(-1)
  double Temp=21.1+273.15;//Density at 21.1°C   
  
  double DH2=10*MolMass_pr*2.0/(Temp*R_const);       // g/cm3; H2 gas density at 21.1°C (70°F) 10atm;
  double DN2=10*MolMass_n14*2.0/(Temp*R_const);      // g/cm3; N2 gas density at 21.1°C (70°F) 10atm; 
  double D3He=10*MolMass_he3/(Temp*R_const);         // g/cm3; 3He gas density at 21.1°C (70°F) 10atm;
  
  if(DEBUG>=3) cout<<"  density_he3_@_12amg="<<Dhe3<<" g/cm^3"<<", (according to ideal gas law, density_@_12_atm_294.25k="<<D3He*1.2<<")"<<endl;
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
  double pI_na=pBeamCurrent*1000.0;
  double pThickXDens;// in g/cm2
  double pLumi=0, pInteXs, pInteRate;

  cout <<"\n Detector= "<<pDetectorName<<",  Beam_current= "<<pBeamCurrent<<" uA,  Angle= "<<pDetectorAngle/deg<<" deg"<<endl;

  //Loop1, 
  const char *L1Name[]={"C12","He3","GE-180","GE-180","H2","N2","He3"};
  double L1Z[]={6., 2., double(-Z_ge180), double(-Z_ge180), 1.0, 7.0, 2.0};  //number of protons in one atom, for GE180, use negative z
  double L1N[]={6., 1., double( N_ge180), double( N_ge180), 0.0, 7.0, 1.0};  //number of neutrons in one atom
  int L1Nuclei_per_molecule[]={1, 1, 1, 1, 2, 2, 1};  //number of atoms per molecule, use 1 for all compounds, n for Hn, Dn
  double L1Dens[]={Dc12, Dhe3, Dge180, Dge180, DH2, DN2, D3He};  //g/cm3
  double L1MolMass[]={MolMass_c12, MolMass_he3, MolMass_ge180, MolMass_ge180, MolMass_pr*2, MolMass_n14*2, MolMass_he3};  //g
  double L1Thick[]={0.254, 40.0, 0.014, 0.014, 40.0, 40.0, 40.0}; //cm
  double L1VZ[]={0.0, 0.0, -20.0, 20.0, 0.0, 0.0, 0.0};  // vertex location in cm

  //x an Q2 binning
  //v/cr vxcen(15) R 0.250 0.300 0.350 0.400 0.450 0.500 0.550 0.600 0.650 0.705 0.765 0.825 0.885 0.945 1.005
  //v/cr vxbin(15) R 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.025 0.030 0.030 0.030 0.030 0.030 0.030
  double L1Xbj[]={0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.735, 0.795, 0.855, 0.915, 0.975};  // x boundary
  double L1Q2[]={2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0};  //Q2 boundary


  cout.setf(ios_base::fixed); 
  cout.precision(4);
  cout<<setw(4)<<"Det"<<setw(8)<<"Target"
      <<setw(10)<<"BeamE"<<setw(10)<<"Current"<<setw(10)<<"DetP0"
      <<setw(10)<<"DetAngle"<<setw(10)<<"VZ(cm)"<<setw(10)<<" Thick(cm)"
      <<setw(12)<<"Lumi(10^33)"<<setw(12)<<"InteXS(pb)"<<setw(12)<<"Rate(Hz)"
      <<endl;
  //just for compare to method2, in this case: pElasOnly<0, inelastic only, only calculate rates for i>=4
  int istart = (pElasOnly<0)?4:0;
  for(int i=istart;i<7;i++) {
    cout.precision(3);
    pThickXDens=L1Dens[i]*L1Thick[i]; // (g/cm2);
    pLumi=GetLumi10pow30(pI_na,L1Nuclei_per_molecule[i],L1MolMass[i],pThickXDens);
    
    pInteXs = 0.0;
    
    //do C12 or pressure curve only for 1-pass elastic kinematic points
    if((i==0 || i>=4)  && !(fabs(pBeamE-2.1)<0.1 && pDetectorMomentum>2.0)) continue;
    
    //I will also do x binning for He3 DIS run
    if(i==1 && pBeamE>8.0) {
      //He3, calculate xs for each x-z bin
      for(int ix=0;ix<14;ix++) {
      
        int nTry = 20;
        double pXS_he3 = 0, pXS_tot_he3 = 0;
        for(int iz=0;iz<nTry;iz++) {
          double pVZ = -0.5*L1Thick[i] + (iz+0.5)* (L1Thick[i]/nTry);
          
          if(pElasOnly==2)
            pXS_he3=GetInteXS(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], pVZ, pDetectorName, pElasOnly,L1Xbj[ix],L1Xbj[ix+1],-1.,-1.,1.10,1.35);
          else if(pElasOnly==4)
            pXS_he3=GetInteXS(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], pVZ, pDetectorName, pElasOnly,L1Xbj[ix],L1Xbj[ix+1],-1.,-1.,2.00,100.0);
          else
            pXS_he3=GetInteXS(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], pVZ, pDetectorName, pElasOnly, L1Xbj[ix],L1Xbj[ix+1]);

          if(DEBUG>=3) cout<<"  x = "<<(L1Xbj[ix]+L1Xbj[ix+1])/2.0<<" +/- "<<(L1Xbj[ix+1]-L1Xbj[ix])/2.0<<"  VZ = "<<setw(6)<<pVZ<<"  InteXS_he3 = "<<setw(8)<<pXS_he3*1.0E6<<" (pb)"<<endl;
          pXS_tot_he3 += pXS_he3;
        }
        if(DEBUG>=3) cout<<"  Averaged  InteXS_he3 = "<<pXS_tot_he3/nTry*1.0E6<<" (pb)"<<endl;
      
        double pInteXs_xq=pXS_tot_he3/nTry;  //Need to get average, or change the luminosity during integration
        double pRate_xq=pLumi*pInteXs_xq;
        if(DEBUG>=2) cout<<"  x="<<(L1Xbj[ix]+L1Xbj[ix+1])/2.0<<" +/- "<<(L1Xbj[ix+1]-L1Xbj[ix])/2.0<<"  pInteXs_xq = "<<setw(8)<<pInteXs_xq*1.0E6<<" (pb)  Rate="<<pRate_xq<<" (Hz)"<<endl;

        pInteXs += pInteXs_xq;  
      }
    }
    else if (i==1 || i>=4) {
      //long target, calculate xs for each z bin     
      int nTry = 20;
      double pXS_long = 0, pXS_tot_long = 0;
#ifdef CMP_GETINTEXS        
      double pXS_long_old, pXS_tot_long_old=0;
#endif      
      //for He3 target, calculate rate for each x bin
      for(int iz=0;iz<nTry;iz++) {
        double pVZ = -0.5*L1Thick[i] + (iz+0.5)* (L1Thick[i]/nTry);
        if(pElasOnly==2)
          pXS_long=GetInteXS(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], pVZ, pDetectorName, pElasOnly,-1.,-1.,-1.,-1.,1.10,1.35);
        else if(pElasOnly==4)
          pXS_long=GetInteXS(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], pVZ, pDetectorName, pElasOnly,-1.,-1.,-1.,-1.,2.00,100.0);
        else
          pXS_long=GetInteXS(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], pVZ, pDetectorName, pElasOnly);
        if(DEBUG>=3) cout<<"  VZ = "<<setw(6)<<pVZ<<"  InteXS_"<<L1Name[i]<<" = "<<setw(8)<<pXS_long*1.0E6<<" (pb)"<<endl;
        pXS_tot_long += pXS_long;

        //,,,,,,,,,,,,,,,,compare old and new GetInteXS() start ,,,,,,,,,,,,,,,,,
#ifdef CMP_GETINTEXS        
        if(DEBUG>=3) {
          if(pElasOnly==2)
            pXS_long_old=GetInteXS_old(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], pVZ, pDetectorName, pElasOnly,-1.,-1.,-1.,-1.,1.10,1.35);
          else if(pElasOnly==4)
            pXS_long_old=GetInteXS_old(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], pVZ, pDetectorName, pElasOnly,-1.,-1.,-1.,-1.,2.00,100.0);
          else
            pXS_long_old=GetInteXS_old(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], pVZ, pDetectorName, pElasOnly);
          pXS_tot_long_old += pXS_long_old;
          cout<<" New - OLD:  VZ = "<<setw(6)<<pVZ<<"  InteXS_"<<L1Name[i]<<" = "<<setw(8)<<(pXS_long-pXS_long_old)*1.0E6<<" (pb)"<<endl;
        }
#endif        
        //,,,,,,,,,,,,,,,,,,compare old and new GetInteXS() end,,,,,,,,,,,,,,,,,,
      }
      if(DEBUG>=3) cout<<"  Averaged  InteXS_"<<L1Name[i]<<" = "<<pXS_tot_long/nTry*1.0E6<<" (pb)"<<endl;
#ifdef CMP_GETINTEXS        
      if(DEBUG>=3) cout<<"  OLD  Averaged  InteXS_"<<L1Name[i]<<" = "<<pXS_tot_long_old/nTry*1.0E6<<" (pb)"<<endl;
#endif      
    
      double pInteXs_long=pXS_tot_long/nTry;  //Need to get average, or change the luminosity during integration
      double pRate_long=pLumi*pInteXs_long;
      if(DEBUG>=3) cout<<"  pInteXs_"<<L1Name[i]<<" = "<<setw(8)<<pInteXs_long*1.0E6<<" (pb)  Rate="<<pRate_long<<" (Hz)"<<endl;

      pInteXs += pInteXs_long;
    }
    else {
      if(pElasOnly==2)
        pInteXs=GetInteXS(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], L1VZ[i], pDetectorName, pElasOnly,-1.,-1.,-1.,-1.,1.10,1.35);
      else if(pElasOnly==4)
        pInteXs=GetInteXS(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], L1VZ[i], pDetectorName, pElasOnly,-1.,-1.,-1.,-1.,2.00,100.0);
      else
        pInteXs=GetInteXS(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], L1VZ[i], pDetectorName, pElasOnly);
    }
    pInteRate=pLumi*pInteXs;
  
    cout<<setw(4)<<pDetectorName<<setw(8)<<L1Name[i]
        <<setprecision(4)<<setw(10)<<pBeamE
        <<setw(10)<<setprecision(2)<<pBeamCurrent<<setw(10)<<setprecision(3)<<pDetectorMomentum
        <<setprecision(2)<<setw(10)<<pDetectorAngle/deg<<setw(10)<<L1VZ[i]
        <<setw(10)<<setprecision(3)<<L1Thick[i]
        <<setprecision(2)<<setw(12)<<pLumi/1000.<<" "<<setw(11)<<pInteXs*1.0E6<<" "<<setw(11)<<pInteRate
        <<endl;

    //print out the pressure curve rates, just scale them
    if(DEBUG>=1) {
      if(i>=4 && pDetectorMomentum>2.0 && pDetectorMomentum<2.2) {
        // Elastic (P_bP_t) pressure curve
        for(int iden=1;iden<=5;iden++) {
          cout.precision(3);
          cout<<"  ["<<L1Name[i]<<"] Pressure ="<<setw(6)<<2.0*iden<<" (atm)"<<" Lumi="<<setw(8)<<pLumi/1000.*0.2*iden<<" (10^33)"
              <<" pInteXs="<<setw(8)<<pInteXs*1.0E6<<" (pb)"<<" Rate="<<setw(8)<<pInteXs*pLumi*0.2*iden<<" (Hz)"<<endl;
        }
      }
    } 
  }   
  cout<<endl;
  cout.precision(pre);    //recover cout default precision
  cout.unsetf(std::ios::floatfield);  //recover cout floatfield properties
  
  return 0;
}


////////////////////////////////////////////////////////////////////////////
void A1NRates()
{
  const double kBeamI[]={1.0,1.0,30.0,30.0,30.0,30.0};  //in uA
  const double kBeamE[]={2.1,2.1,10.5,10.5,10.5,10.5};
  const double kHMSAngle[]={11.7,11.7,12.5,12.5,30.0,30.0};
  const double kHMSP0[]={2.068,1.682,5.7,6.8,2.9,3.5};
  const double kSHMSAngle[]={8.5,8.5,12.5,12.5,30.0,30.0};
  const double kSHMSP0[]={2.083,1.718,5.8,7.5,2.4,3.4};
  
  double pRateHMS=0.0,pRateSHMS=0.0;
  
  for(int j=0;j<6;j++) {
    //double GetRate(double pBeamE, double pDetectorAngle, double pDetectorMomentum, string pDetectorName)
    pRateHMS = GetRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"HMS",0);  
    if(j!=1) pRateSHMS = GetRate(kBeamI[j],kBeamE[j],kSHMSAngle[j]*deg,kSHMSP0[j],"SHMS",0);
    if(j==0) {
      //pure elas
      pRateSHMS = GetRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"HMS",1);
      pRateSHMS = GetRate(kBeamI[j],kBeamE[j],kSHMSAngle[j]*deg,kSHMSP0[j],"SHMS",1);
      //pure inelastic
      //pRateHMS = GetRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"HMS",-1);  
      //pRateSHMS = GetRate(kBeamI[j],kBeamE[j],kSHMSAngle[j]*deg,kSHMSP0[j],"SHMS",-1);    
      //2-sc-bar only
      pRateSHMS = GetRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"HMS",30);
      pRateSHMS = GetRate(kBeamI[j],kBeamE[j],kSHMSAngle[j]*deg,kSHMSP0[j],"SHMS",30);
      //pure elas and 2-sc-bar only
      pRateSHMS = GetRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"HMS",31);
      pRateSHMS = GetRate(kBeamI[j],kBeamE[j],kSHMSAngle[j]*deg,kSHMSP0[j],"SHMS",31);
    }
    if(j==1) {
      //apply W cuts fir Delta peak
      pRateSHMS = GetRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"HMS",2);
      if(j!=1) pRateSHMS = GetRate(kBeamI[j],kBeamE[j],kSHMSAngle[j]*deg,kSHMSP0[j],"SHMS",2);
    }
    if(j>=2) {
      //apply W cuts for DIS
      pRateSHMS = GetRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"HMS",4);
      pRateSHMS = GetRate(kBeamI[j],kBeamE[j],kSHMSAngle[j]*deg,kSHMSP0[j],"SHMS",4);
    }
  }
}

/*
//check E' @ W=1232 peak
double W=1.232;
double M=0.9383;
double sin2th = pow(sin(11.7/57.3/2.),2.);
double E0=2.10;
double nu = (W*W-M*M+4*E0*E0*sin2th)/(2*M+4*E0*sin2th);
cout<<" E0="<<E0<<",  nu="<<nu<<"  ==> E'="<<E0-nu<<"\n";

double W=1.232;
double M=0.9383;
double sin2th = pow(sin(8.5/57.3/2.),2.);
double E0=2.10;
double nu = (W*W-M*M+4*E0*E0*sin2th)/(2*M+4*E0*sin2th);
cout<<" E0="<<E0<<",  nu="<<nu<<"  ==> E'="<<E0-nu<<"\n";


double Ep = 1.682;
double Q2 = 4*E0*Ep*sin2th;
double W2 =-Q2+M*M+2*M*(E0-Ep);
cout<<" E0="<<E0<<",  E'="<<Ep<<"  ==> W="<<sqrt(W2)<<"\n";
*/
