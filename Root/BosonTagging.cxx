#include <iostream>
#include <stdexcept>
#include <algorithm>
#include "ZeroLeptonRun2/BosonTagging.h"


bool BosonTagging::ReturnTag(int mycase, float jetpT, float jetM, float jetD2){ 
  

  //std::cout << "MASS AND D2 " << good_jets_recl[j0].M() << "  " << vD2.at(j0) << "  " << good_jets_recl[j0].Pt() << std::endl;                                             
  // MASS CUT                                                                                                                                                                
  // http://acode-browser.usatlas.bnl.gov/lxr/source/atlas/Reconstruction/Jet/JetSubStructureUtils/Root/BosonTag.cxx                                                         
  // mass_param 1-2 ; mass_window ; d2_params 1-5; dcut direction                                                                                                            
  // W-tagging : https://github.com/mileswu/JetSubstructureTools/blob/master/JetSubStructureUtils/data/config_13TeV_20150528_Wtagging.dat                                    
  // medium  AK10LCTRIMF5R20         73.0577  0.0131854 15.0  1.1036    0.000185103 2.7553e-07    0   0   RIGHT                                                              
  // tight  AK10LCTRIMF5R20         73.0577  0.0131854 15.0  0.67947   0.000592321 -8.09529e-08  0   0   RIGHT                                                               
  // Z-tagging : http://acode-browser.usatlas.bnl.gov/lxr/source/atlas/Reconstruction/Jet/JetSubStructureUtils/data/config_13TeV_20150528_Ztagging.dat                       
  // medium   AK10LCTRIMF5R20                 84.2200  0.0116900 15.0  1.00952   0.000322370 1.68332e-07  0   0   RIGHT                                                      
  // tight    AK10LCTRIMF5R20                 84.2200  0.0116900 15.0  0.514401  0.000851661 -2.0976e-07  0   0   RIGHT                                                      

  // Mass window
  float massWindow = 15.0;
  // Mass params
  std::vector<float> mass_params;
  // W boson
  if(mycase==1 || mycase==2){
    mass_params.push_back(73.0577);
    mass_params.push_back(0.0131854);
  }
  // Z boson
  if(mycase==3 || mycase==4){
    mass_params.push_back(84.2200);
    mass_params.push_back(0.0116900);
  }
  std::vector<float> d2_params;
  
  // W D2 param - MEDIUM
  if(mycase==1){
    d2_params.push_back(1.1036);
    d2_params.push_back(0.000185103);
    d2_params.push_back(2.7553e-07);
    d2_params.push_back(0);
    d2_params.push_back(0);
  }
  // W D2 param - TIGHT
  if(mycase==2){
    d2_params.push_back(0.67947);
    d2_params.push_back(0.000592321);
    d2_params.push_back(-8.09529e-08);
    d2_params.push_back(0);
    d2_params.push_back(0);
  }
  // Z D2 param - MEDIUM
  if(mycase==3){
    d2_params.push_back(1.00952);
    d2_params.push_back(0.000322370);
    d2_params.push_back(1.68332e-07);
    d2_params.push_back(0);
    d2_params.push_back(0);
  }
  // Z D2 param - TIGHT 
  if(mycase==4){
    d2_params.push_back(0.514401);
    d2_params.push_back(0.000851661);
    d2_params.push_back(-2.0976e-07);
    d2_params.push_back(0);
    d2_params.push_back(0);
  }
  // direction
  std::string direction = "RIGHT";

  // http://acode-browser.usatlas.bnl.gov/lxr/source/atlas/Reconstruction/Jet/JetSubStructureUtils/Root/BosonTag.cxx L.417 and later
  float meanMassVal(0.0);
  for(int i=0;i<2;i++){
    meanMassVal += mass_params.at(i)*pow(jetpT/1.e3,i);
  }
  bool passM = (meanMassVal - massWindow < jetM/1.e3) && (jetM/1.e3 < meanMassVal + massWindow);
  
  // D2 CUT                                                                                                                                                                                        
  float d2CutVal= d2_params.at(0)
    + d2_params.at(1) * jetpT/1.e3
    + d2_params.at(2) * pow(jetpT/1.e3, 2)
    + d2_params.at(3) * pow(jetpT/1.e3, 3)
    + d2_params.at(4) * pow(jetpT/1.e3, 4);
  bool passD2= (jetD2 < d2CutVal && direction=="RIGHT") || (d2CutVal < jetD2 && direction=="LEFT");

  return passM&&passD2;
  
}


