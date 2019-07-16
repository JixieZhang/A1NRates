#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <string>


using namespace std;

extern double D2NRates();
extern double A1NRates();
extern double GetRate(double pBeamCurrent, double pBeamE, double pDetectorAngle, double pDetectorMomentum, string pDetectorName, int pElasOnly=0);


int getrate1_main(int argc, char** argv)
{
  if(argc<5) {
  cout<<" Error: you need to provide at least 5 arguments!\n"
      <<" Calculate rates using method 1.\n"
      <<" Usage: "<<argv[0]<<" <BeamCurrent_uA> <Beam_GeV> <DetectorAngle_deg> <DetectorMomentum_GeV> <DetectorName=HMS|SHMS> [ElasOnly=0]\n"
      <<"        All energies are in GeV unit. All angles are in degree unit.\n"
      <<"        ElasOnly=-1: pure inelastic for full acceptance.\n"
      <<"        ElasOnly=0:  inelastic + elastic for full acceptance.\n"
      <<"        ElasOnly=1:  pure elastic for full acceptance.\n"
      <<"        ElasOnly=2:  inelastic + elastic for full acceptance, with cut of 1.10<W<1.35.\n"
      <<"        ElasOnly=4:  inelastic + elastic for full acceptance. with cut of 2.00<W<100.0.\n"
      <<"        ElasOnly=-30: pure inelastic for 2-SC-Bar acceptance.\n"
      <<"        ElasOnly=30: inelastic + elastic for 2-SC-Bar acceptance.\n"
      <<"        ElasOnly=31: pure elastic for 2-SC-Bar acceptance.\n"
      <<endl;
    exit(-1);
  }
  const double degree = asin(1.0)/90.0;
  double pBeamCurrent, pBeamE, pDetectorAngle, pDetectorMomentum;
  string pDetectorName;

  pBeamCurrent = atof(argv[1]); 
  pBeamE = atof(argv[2]); 
  pDetectorAngle = atof(argv[3])*degree; 
  pDetectorMomentum = atof(argv[4]); 
  pDetectorName = argv[5];
  int pElasOnly = 0;
  if(argc>=6) pElasOnly = atol(argv[6]); 
  
  cout<<"  Beam="<<pBeamE<<"  DetAngle="<<pDetectorAngle/degree<<"  pDetectorMomentum="<<pDetectorMomentum<<endl;

  GetRate(pBeamCurrent, pBeamE, pDetectorAngle, pDetectorMomentum, pDetectorName, pElasOnly);
  return 0;
}


int main(int argc, char** argv)
{
  if(argc<2) {
    cout<<"\n Calcualte rates for A1N or D2N. There are 2 modes to run this program:\n\n";
    cout<<" Usage 1: "<<argv[0]<<" <task=-1|-2> \n"
        <<"          This will calculate rates for all kinematic points for A1N(task=-1) or d2n(task=-2)\n\n";
    cout<<" Usage 2: "<<argv[0]<<" <BeamCurrent_uA> <Beam_GeV> <DetectorAngle_deg> <DetectorMomentum_GeV> <DetectorName=HMS|SHMS> [ElasOnly=0]\n"
        <<"        All energies are in GeV unit. All angles are in degree unit.\n"
        <<"        ElasOnly=-1: pure inelastic for full acceptance.\n"
        <<"        ElasOnly=0:  inelastic + elastic for full acceptance.\n"
        <<"        ElasOnly=1:  pure elastic for full acceptance.\n"
        <<"        ElasOnly=2:  inelastic + elastic for full acceptance, with cut of 1.10<W<1.35.\n"
        <<"        ElasOnly=4:  inelastic + elastic for full acceptance. with cut of 2.00<W<100.0.\n"
        <<"        ElasOnly=-30: pure inelastic for 2-SC-Bar acceptance.\n"
        <<"        ElasOnly=30: inelastic + elastic for 2-SC-Bar acceptance.\n"
        <<"        ElasOnly=31: pure elastic for 2-SC-Bar acceptance.\n"
        <<endl;
    return -1;
  }
  else if(argc==2) {
    //do calculation for all points
    int task = atol(argv[1]);
    if(task==-1) A1NRates();
    else if(task==-2) D2NRates();
    else cout<<" Error: this task is not supported yet ...\n";
  } else {
    // do calculation for only one point
    getrate1_main(argc,argv);
  }
  return 0;
}

/*
./a1nrates >! result.txt;


#elastic (PbPt point)
./a1nrates 1.0 2.1 11.7 2.068 HMS  30 >! result_elas.txt;
./a1nrates 1.0 2.1 11.7 2.068 HMS  31 >> result_elas.txt;
./a1nrates 1.0 2.1 8.5 2.083 SHMS  30 >> result_elas.txt;
./a1nrates 1.0 2.1 8.5 2.083 SHMS  31 >> result_elas.txt;

#Delta
./a1nrates 1.0 2.1 11.7 1.682 HMS  0 >! result_delta.txt;
./a1nrates 1.0 2.1 11.7 1.682 HMS  2 >> result_delta.txt;

#Point A|B|C|D 
./a1nrates 30.0 10.5 12.5 5.700  HMS 0 >! result_ABCD.txt;
./a1nrates 30.0 10.5 12.5 5.800 SHMS 0 >> result_ABCD.txt;
./a1nrates 30.0 10.5 12.5 6.800  HMS 0 >> result_ABCD.txt;
./a1nrates 30.0 10.5 12.5 7.500 SHMS 0 >> result_ABCD.txt;
./a1nrates 30.0 10.5 30.0 2.900  HMS 0 >> result_ABCD.txt;
./a1nrates 30.0 10.5 30.0 2.400 SHMS 0 >> result_ABCD.txt;
./a1nrates 30.0 10.5 30.0 3.500  HMS 0 >> result_ABCD.txt;
./a1nrates 30.0 10.5 30.0 3.400 SHMS 0 >> result_ABCD.txt;

#Point A|B|C|D with W>2.0 cut
./a1nrates 30.0 10.5 12.5 5.700  HMS 4 >! result_ABCD_Wgt2.txt;
./a1nrates 30.0 10.5 12.5 5.800 SHMS 4 >> result_ABCD_Wgt2.txt;
./a1nrates 30.0 10.5 12.5 6.800  HMS 4 >> result_ABCD_Wgt2.txt;
./a1nrates 30.0 10.5 12.5 7.500 SHMS 4 >> result_ABCD_Wgt2.txt;
./a1nrates 30.0 10.5 30.0 2.900  HMS 4 >> result_ABCD_Wgt2.txt;
./a1nrates 30.0 10.5 30.0 2.400 SHMS 4 >> result_ABCD_Wgt2.txt;
./a1nrates 30.0 10.5 30.0 3.500  HMS 4 >> result_ABCD_Wgt2.txt;
./a1nrates 30.0 10.5 30.0 3.400 SHMS 4 >> result_ABCD_Wgt2.txt;

*/
