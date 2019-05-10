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

using namespace std;

#define DEBUG 3

static const double deg = acos(0.0)/90.0;

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

////////////////////////////////////////////////////////////////////////////
//get integrated xs for given element at given vertex z,
//applying x and q2 cuts if their upper limits are larger than zero
double GetInteXS(double pBeamE, double pAngle, double pMomentum, double Z, int N, double VZ, string Name, int ElasOnly=0,
                 double pXmin=-1.0, double pXmax=-1.0, double pQ2min=-1.0, double pQ2max=-1.0)
{
  //for HMS
  double pTheta_tg_max=0.07;
  double pTheta_tg_min=-0.06;
  double pPhi_tg_max=0.024;
  double pPhi_tg_min=-0.024;
  double pEprime_max=pMomentum*1.15;
  double pEprime_min=pMomentum*0.88;
  //HMS entrance information, from https://github.com/JeffersonLab/mc-single-arm/blob/master/src/hms/mc_hms.f
  double pPivot2Entrance=126.2;      //cm
  double pH_entr_min=-4.575+0.000;   //cm
  double pH_entr_max=4.575+0.000;    //cm
  double pV_entr_min=-12.144+0.028;  //cm
  double pV_entr_max=12.144+0.028;   //cm

  if(Name=="SHMS") {
    pTheta_tg_max=0.05;
    pTheta_tg_min=-0.05;
    pPhi_tg_max=0.04;
    pPhi_tg_min=-0.04;
    pEprime_max=pMomentum*1.39;
    pEprime_min=pMomentum*0.79;
    //SHMS entrance, https://github.com/JeffersonLab/mc-single-arm/blob/master/src/shms/mc_shms.f
    pPivot2Entrance=258.0;    //cm
    pH_entr_min=-8.5+0.000;   //cm
    pH_entr_max=8.5+0.000;    //cm
    pV_entr_min=-12.5+0.000;  //cm
    pV_entr_max=12.5+0.000;   //cm
  }
  if(pEprime_max>pBeamE) pEprime_max=pBeamE;
    
  const double AMU = 0.9314941;
  double pMtg=(Z+N)*AMU;
  double pP_elas=0.0;

  double pTheta_tg,pPhi_tg=0;
  double pTheta=-999.,pPhi=-999.,pEprime=-999.;
  double pXs_elas=0.0,pXs=0.0,pInteXs=0.0;

  double pTheta_min = pAngle - 20*deg;
  double pTheta_max = pAngle + 20*deg;
  if(pTheta_min < 5*deg) pTheta_min = 5*deg;
  
  double pPhi_min = -80 *deg ;
  double pPhi_max =  80 *deg ;

  double dTheta = 0.2*deg;
  double dPhi = 0.2*deg;
  double dEprime = 0.01;
  
  double x_entr,y_entr,z_entr=0;

  pTheta = pTheta_min - 0.5*dTheta;
  while (pTheta < pTheta_max) {
    
    pTheta += dTheta;
    //get elastic Eprime 
    pP_elas = GetElasEprime(pBeamE, pTheta, pMtg);
    double dOmega = (cos(pTheta-0.5*dTheta)-cos(pTheta+0.5*dTheta))*dPhi;

    pPhi = pPhi_min - 0.5*dPhi;
    while (pPhi < pPhi_max) {
      pPhi += dPhi;

      //convert Lab frame to Tranportation frame
      //void  P_HCS2TCS(double Theta_hall, double Phi_hall, double EndPlaneTheta_hall, double &Theta_tr, double &Phi_tr)
      Transform::P_HCS2TCS(pTheta, pPhi, pAngle, pTheta_tg, pPhi_tg);
      if(pTheta_tg>pTheta_tg_max || pTheta_tg<pTheta_tg_min)  continue;
      if(pPhi_tg>pPhi_tg_max || pPhi_tg<pPhi_tg_min)  continue;
      
      //project to entrance plane
      //SHMS has a honrizontal bending field, but I assume this field do not affect out results        
      //void Project(double x,double y,double z,double z_drift,double theta,double phi,double &x_out, double &y_out, double &z_out)
      x_entr=y_entr=z_entr=0;
      
      //void  X_HCS2TCS(double x, double y, double z,double EndPlaneTheta_hall, double &x_tr, double &y_tr, double &z_tr)
      double x_tr=0, y_tr=0, z_tr=0;
      Transform::X_HCS2TCS(0,0,VZ,pAngle,x_tr,y_tr,z_tr);
      Transform::Project(x_tr, y_tr, z_tr, pPivot2Entrance-z_tr,pTheta_tg,pPhi_tg,x_entr,y_entr,z_entr);
      if(DEBUG>=5) cout<<" x_tr="<<x_tr<<"  y_tr="<<y_tr<<"  z_tr="<<z_tr<<"  theta_tr="<<pTheta_tg<<"  phi_tr="<<pPhi_tg<<endl;
      if(x_entr<pV_entr_min || x_entr>pV_entr_max) continue;
      if(Name=="HMS" || Name=="SHMS") {
        if(y_entr<pH_entr_min || y_entr>pH_entr_max) continue;
      }
        
      pEprime = pEprime_min - 0.5*dEprime;
      double pEprime_up = pEprime + 0.5*dEprime;
      while (pEprime < pEprime_max) {
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
        
        double pQ2=2.0*pBeamE*pEprime*(1.0-cos(pTheta));
        double pXbj=pQ2/(2.0*0.9383*(pBeamE-pEprime));
        //applying X cuts
        if(pXmax>0.0 && pXmax>pXmin) {
          if(pXbj<pXmin || pXbj>pXmax) continue;
        }
        //applying Q2 cuts
        if(pQ2max>0.0 && pQ2max>pQ2min) {
          if(pQ2<pQ2min || pQ2>pQ2max) continue;
        }
        if(DEBUG>=4) cout<<" pTheta="<<pTheta/deg<<"  pPhi="<<pPhi/deg<<"  pEprime="<<pEprime<<endl;

        if(!ElasOnly) {
          //please note that PBosted xs is for per nucleus
          //double PBosted::GetXS(int Z, int N, double Ei, double Ef, double theta, double Tb=-0.001, double Ta=-0.001);
          pXs=0.0;
          //note that PBosted::GetXS() will not work for Q2>11, I give up these points
          if(pQ2 < 11.0) pXs = PBosted::GetXS(Z, N, pBeamE, pEprime, pTheta, 0.000, 0.000);
          if(isnanf(pXs)) pXs=0.0;
          pXs *= 1000.0;   //turn from ub/MeV/Sr into ub/GeV/Sr
          if(DEBUG>=4) cout<<" Xbj = "<< pXbj<<"  Q2 = "<<pQ2<<"  PBosted::GetXS() = "<<pXs*1.0E3<<" (nb/GeV/Sr)"<<endl;
          pInteXs += pXs*deltaEprime*dOmega;
        }
        
        //check if elastic events are accepted or not
        //please note that elas xs is for per nucleus
        if(fabs(pP_elas-pEprime)<0.5*deltaEprime) {
          //double ElasModel::GetXS(int Z, int N, double Ei, double Theta, double Mtg_GeV=0.0, int iUseLarry=0);
          pXs_elas=0.0;
          pXs_elas = ElasModel::GetXS(Z,N,pBeamE,pTheta,pMtg);
          if(isnanf(pXs_elas)) pXs_elas=0.0;
          pXs_elas *= 1.0;   //turn into ub/Sr
          if(DEBUG>=4) cout<<" pP_elas="<<pP_elas<<" deltaEprime="<<deltaEprime<<",  ElasModel::GetXS() = "<<pXs_elas*1.0E6<<" (pb/Sr)"<<endl;
          pInteXs += pXs_elas*dOmega;
        }
      }
    }    
  }
  
  return pInteXs;
}

////////////////////////////////////////////////////////////////////////////
//Beam Current in uA, All energies are in GeV unit. All angles are in radian unit.
double GetRate(double pBeamCurrent, double pBeamE, double pDetectorAngle, double pDetectorMomentum, string pDetectorName, int pElasOnly=0)
{
  if(pElasOnly)  cout<<"\n===================pure elastic======================\n";
  
  int pre = cout.precision();   //get the original precision of cout
  //define material 

  //const double amu = 0.93149403;//molar mass constant in GeV/c2

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

//density=2.77g/cm3 ,   rad-thick = 19.4246 g/cm2, rad-length=7.038cm
//The average moler mass of a mixture is given by
// 1/M = sum(W_i/M_i)
// where W_i is the mass fration and M_i is the moler mass of the component i 

  double MolMass_pr=1.00794;    //in g/mol
  double MolMass_he3=3.01603;   //in g/mol
  double MolMass_c12=12.01078;  //in g/mol
  //double MolMass_ta=180.94788;  //in g/mol
  
  double Dpr=0.07085;  // g/cm3;  //LH2 density 70.85 g/L
  double Dc12=2.267;   // g/cm3;
  double Dhe3=0.00178; // g/cm3;  //Density at 21.1°C (70°F) : 0.1650 kg/m3
  //double Dta=16.69;   // g/cm3;

  double Dhe3_294p25K = 0.165E-3;  //in unit of g/cm3
  double Dhe3_274p15K = Dhe3_294p25K * 294.25/273.15;
  double Dhe3_10amg = Dhe3_274p15K * 10.0; //density at 10 amagats, 
  Dhe3 = Dhe3_10amg;  //does not depends on temperature any more since the volumn is sealed up
  if(DEBUG>=1) cout<<"  density_he3_@_10amg="<<Dhe3<<" g/cm^3"<<endl;
  
  double MolMass_SiO2=(28.0856*1+15.9994*2);  //in g/mol   
  double MolMass_BaO =(137.327*1+15.9994*1);  //in g/mol
  double MolMass_Al2O3=(26.982*2+15.9994*3);  //in g/mol
  double MolMass_CaO =(40.078*1+15.9994*1);   //in g/mol
  double MolMass_SrO =(87.62*1+15.9994*1);    //in g/mol

  double pMolMass_aver_ge180_rev = 0.603/MolMass_SiO2;
  pMolMass_aver_ge180_rev += 0.182/MolMass_BaO;
  pMolMass_aver_ge180_rev += 0.143/MolMass_Al2O3;
  pMolMass_aver_ge180_rev += 0.065/MolMass_CaO;
  pMolMass_aver_ge180_rev += 0.0025/MolMass_SrO;
  pMolMass_aver_ge180_rev += 0.0045/MolMass_pr;
  
  double pMolMass_aver_ge180 = 1.0 / pMolMass_aver_ge180_rev;
  double pZ_aver_ge180 = pMolMass_aver_ge180 * 0.4829;
  int pN_aver_ge180 = int(pZ_aver_ge180+0.5);  //round up to 1 if > .5
  
  double Dge180 = 2.77;  // g/cm3;
  if(DEBUG>=1) cout<<"  pZ_aver_ge180="<<pZ_aver_ge180<<"    pMolMass_aver_ge180="<<pMolMass_aver_ge180<<endl;

  double pI_na=pBeamCurrent*1000.0;
  double pThickXDens;// in g/cm2
  double pLumi=0, pInteXs, pInteRate;

  cout<<"\nBeam current = "<<pBeamCurrent<<" uA,  Detector = "<<pDetectorName<<",  Angle = "<<pDetectorAngle/deg<<" deg"<<endl;

  //Loop1, 
  const char *L1Name[]={"C12","He3","GE-180","GE-180","LH2"};
  double L1Z[]={6., 2., double(pZ_aver_ge180), double(pZ_aver_ge180), 1.0};  //number of protons in one atom
  double L1N[]={6., 1., double(pN_aver_ge180), double(pN_aver_ge180), 0.0};  //number of neutrons in one atom
  int L1Nuclei_per_molecule[]={1, 1, 1, 1, 2};  //number of atoms per molecule, use 1 for all compounds, n for Hn, Dn
  double L1Dens[]={Dc12, Dhe3, Dge180, Dge180, Dpr};  //g/cm3
  double L1MolMass[]={MolMass_c12, MolMass_he3, pMolMass_aver_ge180, pMolMass_aver_ge180, MolMass_pr*2};  //g
  double L1Thick[]={0.254, 40.0, 0.014, 0.014, 5.0}; //cm
  double L1VZ[]={0.0, 0.0, -20.0, 20.0, 0.0};  // vertex location in cm

  //x an Q2 binning
  double L1Xbj[]={0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775, 0.825};  // x boundary
  double L1Q2[]={2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0};  //Q2 boundary


  cout.setf(ios_base::fixed); 
  cout.precision(4);
  cout<<setw(4)<<"Det"<<setw(8)<<"Target"
      <<setw(10)<<"BeamE"<<setw(10)<<"Current"<<setw(10)<<"DetP0"
      <<setw(10)<<"DetAngle"<<setw(10)<<"VZ(cm)"<<setw(10)<<" Thick(cm)"
      <<setw(12)<<"Lumi(10^33)"<<setw(12)<<"InteXS(pb)"<<setw(12)<<"Rate(Hz)"
      <<endl;

  for(int i=0;i<4;i++)
  {
    pThickXDens=L1Dens[i]*L1Thick[i]; // (g/cm2);
    pLumi=GetLumi10pow30(pI_na,L1Nuclei_per_molecule[i],L1MolMass[i],pThickXDens);

    //I do not want to do z slice for E=2.2 (elastic setting)
    if(i==1 && pBeamE>8.0) {
      //He3, calculate xs for each x-z bin
      for(int ix=0;ix<12;ix++) {
        
        int nTry = 20;
        double pXS_he3 = 0, pXS_tot_he3 = 0;
        //for He3 target, calculate rate for each x bin
        for(int iz=0;iz<nTry;iz++) {
          double pVZ = -0.5*L1Thick[i] + (iz+0.5)* (L1Thick[i]/nTry);
          pXS_he3=GetInteXS(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], pVZ, pDetectorName, pElasOnly, L1Xbj[ix],L1Xbj[ix+1]);
          if(DEBUG>=3) cout<<"  x = "<<(L1Xbj[ix]+L1Xbj[ix+1])/2.0<<"  VZ = "<<setw(6)<<pVZ<<"  InteXS_he3 = "<<setw(8)<<pXS_he3*1.0E6<<" (pb)"<<endl;
          pXS_tot_he3 += pXS_he3;
        }
        pInteXs=pXS_tot_he3/nTry;   //Need to get average, or change the luminosity
        if(DEBUG>=3) cout<<"  Averaged  InteXS_he3 = "<<pInteXs*1.0E6<<" (pb)"<<endl;
        
        double pInteXs_xq=pXS_tot_he3/nTry;
        double pRate_xq=pLumi*pInteXs_xq;
        if(DEBUG>=2) cout<<"  x="<<(L1Xbj[ix]+L1Xbj[ix+1])/2.0<<"  pInteXs_xq = "<<setw(8)<<pInteXs_xq*1.0E6<<" (pb)  Rate="<<pRate_xq<<" (Hz)"<<endl;
  
      }
    } else {
      pInteXs=GetInteXS(pBeamE, pDetectorAngle, pDetectorMomentum, L1Z[i], L1N[i], L1VZ[i], pDetectorName, pElasOnly);
    }
    pInteRate=pLumi*pInteXs;
    
    cout<<setw(4)<<pDetectorName<<setw(8)<<L1Name[i]
        <<setprecision(4)<<setw(10)<<pBeamE
        <<setw(10)<<setprecision(2)<<pBeamCurrent<<setw(10)<<pDetectorMomentum
        <<setprecision(2)<<setw(10)<<pDetectorAngle/deg<<setw(10)<<L1VZ[i]
        <<setw(10)<<setprecision(3)<<L1Thick[i]
        <<setprecision(2)<<setw(12)<<pLumi/1000.<<setw(12)<<pInteXs*1.0E6<<setw(12)<<pInteRate
        <<endl;
  }   
  cout<<endl;
  cout.precision(pre);    //recover cout default precision
  cout.unsetf(std::ios::floatfield);  //recover cout floatfield properties
  
  return 0;
}


////////////////////////////////////////////////////////////////////////////
void A1NRates()
{
  const double kBeamI[]={1.0,30.0,30.0,30.0,30.0};  //in uA
  const double kBeamE[]={2.2,10.5,10.5,10.5,10.5};
  const double kHMSAngle[]={12.5,12.5,12.5,30.0,30.0};
  const double kHMSP0[]={2.15,5.7,6.8,2.9,3.5};
  const double kSHMSAngle[]={12.5,12.5,12.5,30.0,30.0};
  const double kSHMSP0[]={2.15,5.8,7.5,2.4,3.4};
  
  double pRateHMS=0.0,pRateSHMS=0.0;
  
  for(int j=0;j<5;j++)
  {
    //double GetRate(double pBeamE, double pDetectorAngle, double pDetectorMomentum, string pDetectorName)
    pRateHMS = GetRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"HMS",0);  
    pRateSHMS = GetRate(kBeamI[j],kBeamE[j],kSHMSAngle[j]*deg,kSHMSP0[j],"SHMS",0);   
    if(j==0) {
      pRateSHMS = GetRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"DHMS",1);
      pRateSHMS = GetRate(kBeamI[j],kBeamE[j],kHMSAngle[j]*deg,kHMSP0[j],"SHMS",1);
    }
  }
}

