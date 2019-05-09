#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <time.h>

//#define RATES_DEBUG 1

using namespace std;

extern double A1NRates();
extern double GetRate(double pBeamCurrent, double pBeamE, double pDetectorAngle, double pDetectorMomentum, string pDetectorName);

void usage(int argc, char** argv)
{
  cout<<" Error: you need to provide 5 arguments!\n"
      <<" Usage: "<<argv[0]<<" <pBeamCurrent_uA> <pBeam_GeV> <pDetectorAngle_deg> <pDetectorMomentum_GeV> <pDetectorName=HMS|SHMS>\n"
      <<"        All energies are in GeV unit. All angles are in degree unit."
      <<endl;
}

int _main(int argc, char** argv)
{
  if(argc<2) {
    usage(argc,argv);
    exit(-1);
  }
  const double degree = acos(0.0)/180.0;
  double pBeamCurrent, pBeamE, pDetectorAngle, pDetectorMomentum;
  string pDetectorName;

  pBeamCurrent = atol(argv[1]); 
  pBeamE = atol(argv[2]); 
  pDetectorAngle = atol(argv[3])*degree; 
  pDetectorMomentum = atol(argv[4]); 
  pDetectorName = argv[5];

  GetRate(pBeamCurrent, pBeamE, pDetectorAngle, pDetectorMomentum, pDetectorName);
  return 0;
}

int main(int argc, char** argv)
{
  A1NRates();
  return 0;
}
