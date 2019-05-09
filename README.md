Program to calculate rate using XSModel.
XSModel return DiffXS in ub/MeV-Sr for inelas or ub for elas.
 
To be improved:
1) HMS and SHMS acceptance.
   One need to run single-arm-simulation program to get the angle acceptance. 
   This angle acceptance should have a z dependence.  I am not sure if the
   single-arm-simulation program has already taken z dependence into account.
   If not, one need to develop a transportation packing for single-arm-simulation
   program first.

2) Make A1NRates.cc as a class such that it can be used in root cint.

First release version v1.0.0.

1. Hard-coded HMS and SHMS collimator size and angle acceptance, need to be improved.
2. Hard-coded target material, it will be better to separate them.

