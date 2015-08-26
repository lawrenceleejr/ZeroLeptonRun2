#ifndef BOSONTAGGING_h
#define BOSONTAGGING_h

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cctype>
#include <algorithm>
#include <stdexcept>
#include <exception>
#include <TROOT.h>
#include <vector>
#include <sstream>
#include <TString.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TArrayD.h>
#include <TVectorD.h>
using namespace std;

class BosonTagging{

 public:
  bool ReturnTag(int mycase, float jetpT, float jetM, float jetD2);
  
 private:

};
#endif
