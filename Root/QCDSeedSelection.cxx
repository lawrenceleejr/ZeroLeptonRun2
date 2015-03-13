#include "ZeroLeptonRun2/QCDSeedSelection.h"
#include <cmath>

const double QCDSeedSelection::GeV = 1.e+3;

bool QCDSeedSelection::passTrigger_2011Data(int RunNumber, bool EF_j180_a4_EFFS, bool EF_j135_a4_EFFS, bool EF_j100_a4_EFFS, bool EF_j75_a4_EFFS, bool EF_j55_a4_EFFS, bool EF_j240_a4tc_EFFS,  bool EF_j180_a4tc_EFFS,  bool EF_j135_a4tc_EFFS,  bool EF_j100_a4tc_EFFS,  bool EF_j75_a4tc_EFFS,  bool EF_j55_a4tc_EFFS, int nJets, float leadingJetPt, float& trigger_weight){

  bool passTrigger = false;

  if (RunNumber < 178110){
    if (EF_j180_a4_EFFS && nJets && leadingJetPt > 260*GeV){
      passTrigger = true;
      trigger_weight = 1.0;
    }
    if (EF_j135_a4_EFFS && nJets && leadingJetPt > 200*GeV && leadingJetPt < 260*GeV ){
      passTrigger = true;
      trigger_weight = 1.0;
    }
    else if (EF_j100_a4_EFFS && nJets && leadingJetPt > 160*GeV && leadingJetPt < 200*GeV){
      passTrigger = true;
      trigger_weight = 28.2;
    }
    else if (EF_j75_a4_EFFS && nJets && leadingJetPt > 130*GeV && leadingJetPt < 160*GeV){
      passTrigger = true;
      trigger_weight = 102.6;
    }
    else if (EF_j55_a4_EFFS && nJets && leadingJetPt > 100*GeV && leadingJetPt < 130*GeV){
      passTrigger = true;
      trigger_weight = 380.0;
    }
  }//very early 2011 runs
  else{
    // Periods: E F H I J M 
    if (EF_j240_a4tc_EFFS && nJets && leadingJetPt > 335*GeV){
      passTrigger = true;
      
      if (RunNumber >= 179710 && RunNumber <= 180481 ) // Period D
	trigger_weight = 1.000;  
      
      if (RunNumber >= 180614 && RunNumber <= 180776 ) // Period E
	trigger_weight = 1.000;
      
      if (RunNumber >= 182013 && RunNumber <= 182519 ) // Period F
	trigger_weight = 1.000;
      
      if (RunNumber >= 182726 && RunNumber <= 183462 ) // Period G
	trigger_weight = 1.000;
      
      if (RunNumber >= 183544 && RunNumber <= 184169 ) // Period H
	trigger_weight = 1.000; 
      
      if (RunNumber >= 183544 && RunNumber <= 184169 ) // Period H
	trigger_weight = 1.000;
      
      if (RunNumber >= 185353 && RunNumber <= 186493 ) // Period I
	trigger_weight = 1.000;
      
      if (RunNumber >= 186516 && RunNumber <= 186755 ) // Period J
	trigger_weight = 1.000;
      
      if (RunNumber >= 186873 && RunNumber <= 187815 ) // Period K
	trigger_weight = 1.000;

      if (RunNumber >= 188902 && RunNumber <= 190343 ) // Period L
	trigger_weight = 1.002;

      if (RunNumber >= 190503 && RunNumber <= 191933 ) // Period M
	trigger_weight = 1.001;
    }
    if (EF_j180_a4tc_EFFS && nJets && leadingJetPt > 260*GeV && leadingJetPt <= 335*GeV){
      passTrigger = true;
      
      if (RunNumber >= 179710 && RunNumber <= 180481 ) // Period D
	trigger_weight = 1.000;
      
      if (RunNumber >= 180614 && RunNumber <= 180776 ) // Period E
	trigger_weight = 1.000;
      
      if (RunNumber >= 182013 && RunNumber <= 182519 ) // Period F
	trigger_weight = 1.000;
      
      if (RunNumber >= 182726 && RunNumber <= 183462 ) // Period G
	trigger_weight = 1.041;
      
      if (RunNumber >= 183544 && RunNumber <= 184169 ) //Period H
	trigger_weight = 1.022;
      
      if (RunNumber >= 185353 && RunNumber <= 186493 ) // Period I
	trigger_weight = 14.348;
      
      if (RunNumber >= 186516 && RunNumber <= 186755 ) // Period J
	trigger_weight = 18.611;
      
      if (RunNumber >= 186873 && RunNumber <= 187815 ) // Period K
	trigger_weight = 18.951;
      
      if (RunNumber >= 188902 && RunNumber <= 190343 ) // Period L
	trigger_weight = 27.428;
      
      if (RunNumber >= 190503 && RunNumber <= 191933 ) // Period M
	trigger_weight = 28.502;
    }
    if (EF_j135_a4tc_EFFS && nJets && leadingJetPt > 200*GeV && leadingJetPt <= 260*GeV){
      passTrigger = true;
      
      if (RunNumber >= 179710 && RunNumber <= 180481 ) // Period D
	trigger_weight = 1.937;
      
      if (RunNumber >= 180614 && RunNumber <= 180776 ) // Period E
	trigger_weight = 2.461;
      
      if (RunNumber >= 182013 && RunNumber <= 182519 ) // Period F
	trigger_weight = 37.0;
      
      if (RunNumber >= 182726 && RunNumber <= 183462 ) // Period G
	trigger_weight = 45.484;
      
      if (RunNumber >= 183544 && RunNumber <= 184169 ) //Period H
	trigger_weight = 43.978;
      
      if (RunNumber >= 185353 && RunNumber <= 186493 ) // Period I
	trigger_weight = 58.486;
      
      if (RunNumber >= 186516 && RunNumber <= 186755 ) // Period J
	trigger_weight = 76.808;
      
      if (RunNumber >= 186873 && RunNumber <= 187815 ) // Period K
	trigger_weight = 78.443;
      
      if (RunNumber >= 188902 && RunNumber <= 190343 ) // Period L
	trigger_weight = 117.428;
      
      if (RunNumber >= 190503 && RunNumber <= 191933 ) // Period M
	trigger_weight = 122.267;
    }
    if (EF_j100_a4tc_EFFS && nJets && leadingJetPt > 160*GeV && leadingJetPt <= 200*GeV){
      passTrigger = true;
      
      if (RunNumber >= 179710 && RunNumber <= 180481 ) // Period D
	trigger_weight = 46.915;
      
      if (RunNumber >= 180614 && RunNumber <= 180776 ) // Period E
	trigger_weight = 125.299;
      
      if (RunNumber >= 182013 && RunNumber <= 182519 ) // Period F
	trigger_weight = 141.977;
      
      if (RunNumber >= 182726 && RunNumber <= 183462 ) // Period G
	trigger_weight =  173.904;
      
      if (RunNumber >= 183544 && RunNumber <= 184169 ) //Period H
	trigger_weight = 168.562 ; 
      
      if (RunNumber >= 185353 && RunNumber <= 186493 ) // Period I
	trigger_weight = 221.145;
      
      if (RunNumber >= 186516 && RunNumber <= 186755 ) // Period J
	trigger_weight = 293.661;
      
      if (RunNumber >= 186873 && RunNumber <= 187815 ) // Period K
	trigger_weight = 300.146;
      
      if (RunNumber >= 188902 && RunNumber <= 190343 ) // Period L
	trigger_weight = 448.747;
      
      if (RunNumber >= 190503 && RunNumber <= 191933 ) // Period M
	trigger_weight = 475.2;
    }
    if (EF_j75_a4tc_EFFS && nJets && leadingJetPt > 130*GeV && leadingJetPt <= 160*GeV){
      passTrigger = true;
      
      if (RunNumber >= 179710 && RunNumber <= 180481 ) // Period D
	trigger_weight = 159.0;
      
      if (RunNumber >= 180614 && RunNumber <= 180776 ) // Period E
	trigger_weight = 435.045;
      
      if (RunNumber >= 182013 && RunNumber <= 182519 ) // Period F
	trigger_weight = 495.649;

      if (RunNumber >= 182726 && RunNumber <= 183462 ) // Period G
	trigger_weight = 592.157 ;

      if (RunNumber >= 183544 && RunNumber <= 184169 ) //Period H
	trigger_weight =  571.418;  
      
      if (RunNumber >= 185353 && RunNumber <= 186493 ) // Period I
	trigger_weight = 747.185;
      
      if (RunNumber >= 186516 && RunNumber <= 186755 ) // Period J
	trigger_weight = 996.803;
      
      if (RunNumber >= 186873 && RunNumber <= 187815 ) // Period K
	trigger_weight = 1018.601 ;
      
      if (RunNumber >= 188902 && RunNumber <= 190343 ) // Period L
	trigger_weight = 1594.434;
      
      if (RunNumber >= 190503 && RunNumber <= 191933 ) // Period M
	trigger_weight = 1730.306;
    }
    if (EF_j55_a4tc_EFFS && nJets && leadingJetPt > 100*GeV && leadingJetPt <= 130*GeV){
      passTrigger = true;
      
      if (RunNumber >= 179710 && RunNumber <= 180481 ) // Period D
	trigger_weight = 591.791;
      
      if (RunNumber >= 180614 && RunNumber <= 180776 ) // Period E
	trigger_weight = 1550.372;
      
      if (RunNumber >= 182013 && RunNumber <= 182519 ) // Period F
	trigger_weight = 1745.522;
      
      if (RunNumber >= 182726 && RunNumber <= 183462 ) // Period G
	trigger_weight = 2106.321;
      
      if (RunNumber >= 183544 && RunNumber <= 184169 ) //Period H
	trigger_weight = 2030.037;  
      
      if (RunNumber >= 185353 && RunNumber <= 186493 ) // Period I
	trigger_weight = 2639.199;
      
      if (RunNumber >= 186516 && RunNumber <= 186755 ) // Period J
	trigger_weight = 3541.597;
      
      if (RunNumber >= 186873 && RunNumber <= 187815 ) // Period K
	trigger_weight = 3617.199;
      
      if (RunNumber >= 188902 && RunNumber <= 190343 ) // Period L
	trigger_weight = 6089.262;
      if (RunNumber >= 190503 && RunNumber <= 191933 ) // Period M
	trigger_weight = 6445.025;

    }
  }

  return passTrigger;
}

bool QCDSeedSelection::passTrigger_2012Data(int RunNumber,bool EF_j460_a4tchad, bool EF_j360_a4tchad, bool EF_j280_a4tchad, bool EF_j220_a4tchad, bool EF_j180_a4tchad, bool EF_j145_a4tchad,bool EF_j110_a4tchad,bool EF_j80_a4tchad, bool EF_j55_a4tchad,int nJets, float leadingJetPt, float& trigger_weight){
  
  bool passTrigger = false;

  //The prescale weights are derived using data12_8TeV.periodAllYear_DetStatus-v47-pro13-01_CoolRunQuery-00-04-08_Susy.xml
  //use +80 GeV to define Plateau - 25 April 2012 - BUT need to check with Trigger guys
  //I bet this is wrong, but may be ok as an initial guess to start doing stuff.  
  
  if (EF_j460_a4tchad && nJets && leadingJetPt > 510*GeV){
    passTrigger = true;    
    //Prescales for trigger EF_j460_a4tchad
    //------------------------------------------------------------
    if (RunNumber >= 200804 && RunNumber <= 200913 ) trigger_weight = 1.000 ; // for A1-3
    if (RunNumber >= 200926 && RunNumber <= 201191 ) trigger_weight = 1.000 ; // for A4
    if (RunNumber >= 201257 && RunNumber <= 201556 ) trigger_weight = 1.000 ; // for A5-8
    if (RunNumber >= 202660 && RunNumber <= 203195 ) trigger_weight = 1.000 ; // for B1-3
    if (RunNumber >= 203228 && RunNumber <= 203524 ) trigger_weight = 1.000 ; // for B4
    if (RunNumber >= 203602 && RunNumber <= 203792 ) trigger_weight = 1.000 ; // for B5-6
    if (RunNumber >= 203875 && RunNumber <= 204158 ) trigger_weight = 1.000 ; // for B7-9
    if (RunNumber >= 204240 && RunNumber <= 204668 ) trigger_weight = 1.000 ; // for B10-11
    if (RunNumber >= 204707 && RunNumber <= 205017 ) trigger_weight = 1.003 ; // for B12
    if (RunNumber >= 205055 && RunNumber <= 205113 ) trigger_weight = 1.000 ; // for B13-14
    if (RunNumber >= 206248 && RunNumber <= 206614 ) trigger_weight = 1.000 ; // for C1-3
    if (RunNumber >= 206717 && RunNumber <= 207046 ) trigger_weight = 1.000 ; // for C4-C6
    if (RunNumber >= 207113 && RunNumber <= 207397 ) trigger_weight = 1.000 ; // for C7-9
    if (RunNumber >= 207447 && RunNumber <= 208126 ) trigger_weight = 1.000 ; // for D1-3
    if (RunNumber >= 208179 && RunNumber <= 208485 ) trigger_weight = 1.000 ; // for D4-6
    if (RunNumber >= 208631 && RunNumber <= 209025 ) trigger_weight = 1.000 ; // for D7-8
    if (RunNumber >= 209074 && RunNumber <= 209899 ) trigger_weight = 1.000 ; // for E1-3
    if (RunNumber >= 209909 && RunNumber <= 210308 ) trigger_weight = 1.000 ; // for E4-5
    if (RunNumber >= 211522 && RunNumber <= 211620 ) trigger_weight = 1.000 ; // for G1-G2
    if (RunNumber >= 211670 && RunNumber <= 212142 ) trigger_weight = 1.001 ; // for G3-4
    if (RunNumber >= 212144 && RunNumber <= 212272 ) trigger_weight = 1.000 ; // for G5
    if (RunNumber >= 212619 && RunNumber <= 212858 ) trigger_weight = 1.000 ; // for H1-2
    if (RunNumber >= 212967 && RunNumber <= 213250 ) trigger_weight = 1.000 ; // for H3-4
    if (RunNumber >= 213264 && RunNumber <= 213359 ) trigger_weight = 1.000 ; // for H5-6
    if (RunNumber >= 213431 && RunNumber <= 213819 ) trigger_weight = 1.000 ; // for I1-3
    if (RunNumber >= 213900 && RunNumber <= 214216 ) trigger_weight = 1.000 ; // for J1-2
    if (RunNumber >= 214388 && RunNumber <= 214777 ) trigger_weight = 1.000 ; // for J3-4
    if (RunNumber >= 214984 && RunNumber <= 215061 ) trigger_weight = 1.000 ; // for J5-6
    if (RunNumber >= 215063 && RunNumber <= 215091 ) trigger_weight = 1.000 ; // for J7-8
    if (RunNumber >= 215414 && RunNumber <= 215541 ) trigger_weight = 1.000 ; // for L1
    if (RunNumber >= 215559 && RunNumber <= 215643 ) trigger_weight = 1.000 ; // for L2-3
    //------------------------------------------------------------

    
  }

  if (EF_j360_a4tchad && nJets && leadingJetPt > 410*GeV && leadingJetPt <= 510*GeV){
    passTrigger = true;    
    //Prescales for trigger EF_j460_a4tchad
    //------------------------------------------------------------
    if (RunNumber >= 200804 && RunNumber <= 200913 ) trigger_weight = 1.000 ; // for A1-3
    if (RunNumber >= 200926 && RunNumber <= 201191 ) trigger_weight = 1.000 ; // for A4
    if (RunNumber >= 201257 && RunNumber <= 201556 ) trigger_weight = 1.000 ; // for A5-8
    if (RunNumber >= 202660 && RunNumber <= 203195 ) trigger_weight = 1.000 ; // for B1-3
    if (RunNumber >= 203228 && RunNumber <= 203524 ) trigger_weight = 1.000 ; // for B4
    if (RunNumber >= 203602 && RunNumber <= 203792 ) trigger_weight = 1.000 ; // for B5-6
    if (RunNumber >= 203875 && RunNumber <= 204158 ) trigger_weight = 1.000 ; // for B7-9
    if (RunNumber >= 204240 && RunNumber <= 204668 ) trigger_weight = 1.000 ; // for B10-11
    if (RunNumber >= 204707 && RunNumber <= 205017 ) trigger_weight = 1.003 ; // for B12
    if (RunNumber >= 205055 && RunNumber <= 205113 ) trigger_weight = 1.000 ; // for B13-14
    if (RunNumber >= 206248 && RunNumber <= 206614 ) trigger_weight = 1.000 ; // for C1-3
    if (RunNumber >= 206717 && RunNumber <= 207046 ) trigger_weight = 1.000 ; // for C4-C6
    if (RunNumber >= 207113 && RunNumber <= 207397 ) trigger_weight = 1.000 ; // for C7-9
    if (RunNumber >= 207447 && RunNumber <= 208126 ) trigger_weight = 1.000 ; // for D1-3
    if (RunNumber >= 208179 && RunNumber <= 208485 ) trigger_weight = 1.000 ; // for D4-6
    if (RunNumber >= 208631 && RunNumber <= 209025 ) trigger_weight = 1.000 ; // for D7-8
    if (RunNumber >= 209074 && RunNumber <= 209899 ) trigger_weight = 1.000 ; // for E1-3
    if (RunNumber >= 209909 && RunNumber <= 210308 ) trigger_weight = 1.000 ; // for E4-5
    if (RunNumber >= 211522 && RunNumber <= 211620 ) trigger_weight = 1.000 ; // for G1-G2
    if (RunNumber >= 211670 && RunNumber <= 212142 ) trigger_weight = 1.001 ; // for G3-4
    if (RunNumber >= 212144 && RunNumber <= 212272 ) trigger_weight = 1.000 ; // for G5
    if (RunNumber >= 212619 && RunNumber <= 212858 ) trigger_weight = 1.000 ; // for H1-2
    if (RunNumber >= 212967 && RunNumber <= 213250 ) trigger_weight = 1.000 ; // for H3-4
    if (RunNumber >= 213264 && RunNumber <= 213359 ) trigger_weight = 1.000 ; // for H5-6
    if (RunNumber >= 213431 && RunNumber <= 213819 ) trigger_weight = 1.000 ; // for I1-3
    if (RunNumber >= 213900 && RunNumber <= 214216 ) trigger_weight = 1.000 ; // for J1-2
    if (RunNumber >= 214388 && RunNumber <= 214777 ) trigger_weight = 1.000 ; // for J3-4
    if (RunNumber >= 214984 && RunNumber <= 215061 ) trigger_weight = 1.000 ; // for J5-6
    if (RunNumber >= 215063 && RunNumber <= 215091 ) trigger_weight = 1.000 ; // for J7-8
    if (RunNumber >= 215414 && RunNumber <= 215541 ) trigger_weight = 1.000 ; // for L1
    if (RunNumber >= 215559 && RunNumber <= 215643 ) trigger_weight = 1.000 ; // for L2-3
    //------------------------------------------------------------
  }
  
  if (EF_j280_a4tchad && nJets && leadingJetPt > 330*GeV && leadingJetPt <= 410*GeV){
    passTrigger = true;
    //Prescales for trigger EF_j280_a4tchad
    //------------------------------------------------------------
    if (RunNumber >= 200804 && RunNumber <= 200913 ) trigger_weight = 1.000 ; // for A1-3
    if (RunNumber >= 200926 && RunNumber <= 201191 ) trigger_weight = 1.071 ; // for A4
    if (RunNumber >= 201257 && RunNumber <= 201556 ) trigger_weight = 2.660 ; // for A5-8
    if (RunNumber >= 202660 && RunNumber <= 203195 ) trigger_weight = 17.223 ; // for B1-3
    if (RunNumber >= 203228 && RunNumber <= 203524 ) trigger_weight = 23.417 ; // for B4
    if (RunNumber >= 203602 && RunNumber <= 203792 ) trigger_weight = 24.931 ; // for B5-6
    if (RunNumber >= 203875 && RunNumber <= 204158 ) trigger_weight = 33.630 ; // for B7-9
    if (RunNumber >= 204240 && RunNumber <= 204668 ) trigger_weight = 26.396 ; // for B10-11
    if (RunNumber >= 204707 && RunNumber <= 205017 ) trigger_weight = 30.336 ; // for B12
    if (RunNumber >= 205055 && RunNumber <= 205113 ) trigger_weight = 29.326 ; // for B13-14
    if (RunNumber >= 206248 && RunNumber <= 206614 ) trigger_weight = 23.463 ; // for C1-3
    if (RunNumber >= 206717 && RunNumber <= 207046 ) trigger_weight = 29.532 ; // for C4-C6
    if (RunNumber >= 207113 && RunNumber <= 207397 ) trigger_weight = 26.072 ; // for C7-9
    if (RunNumber >= 207447 && RunNumber <= 208126 ) trigger_weight = 27.015 ; // for D1-3
    if (RunNumber >= 208179 && RunNumber <= 208485 ) trigger_weight = 30.706 ; // for D4-6
    if (RunNumber >= 208631 && RunNumber <= 209025 ) trigger_weight = 27.057 ; // for D7-8
    if (RunNumber >= 209074 && RunNumber <= 209899 ) trigger_weight = 30.835 ; // for E1-3
    if (RunNumber >= 209909 && RunNumber <= 210308 ) trigger_weight = 26.965 ; // for E4-5
    if (RunNumber >= 211522 && RunNumber <= 211620 ) trigger_weight = 21.838 ; // for G1-G2
    if (RunNumber >= 211670 && RunNumber <= 212142 ) trigger_weight = 29.641 ; // for G3-4
    if (RunNumber >= 212144 && RunNumber <= 212272 ) trigger_weight = 28.427 ; // for G5
    if (RunNumber >= 212619 && RunNumber <= 212858 ) trigger_weight = 30.620 ; // for H1-2
    if (RunNumber >= 212967 && RunNumber <= 213250 ) trigger_weight = 31.298 ; // for H3-4
    if (RunNumber >= 213264 && RunNumber <= 213359 ) trigger_weight = 43.595 ; // for H5-6
    if (RunNumber >= 213431 && RunNumber <= 213819 ) trigger_weight = 32.029 ; // for I1-3
    if (RunNumber >= 213900 && RunNumber <= 214216 ) trigger_weight = 29.519 ; // for J1-2
    if (RunNumber >= 214388 && RunNumber <= 214777 ) trigger_weight = 28.286 ; // for J3-4
    if (RunNumber >= 214984 && RunNumber <= 215061 ) trigger_weight = 25.373 ; // for J5-6
    if (RunNumber >= 215063 && RunNumber <= 215091 ) trigger_weight = 27.960 ; // for J7-8
    if (RunNumber >= 215414 && RunNumber <= 215541 ) trigger_weight = 27.449 ; // for L1
    if (RunNumber >= 215559 && RunNumber <= 215643 ) trigger_weight = 34.854 ; // for L2-3
    //------------------------------------------------------------
  }
  
  // two trigger end point
  if (EF_j220_a4tchad && nJets && leadingJetPt > 270*GeV && leadingJetPt <= 330*GeV){
    passTrigger = true;
    //Prescales for trigger EF_j220_a4tchad
    //------------------------------------------------------------
    if (RunNumber >= 200804 && RunNumber <= 200913 ) trigger_weight = 1.000 ; // for A1-3
    if (RunNumber >= 200926 && RunNumber <= 201191 ) trigger_weight = 6.186 ; // for A4
    if (RunNumber >= 201257 && RunNumber <= 201556 ) trigger_weight = 26.356 ; // for A5-8
    if (RunNumber >= 202660 && RunNumber <= 203195 ) trigger_weight = 59.186 ; // for B1-3
    if (RunNumber >= 203228 && RunNumber <= 203524 ) trigger_weight = 86.725 ; // for B4
    if (RunNumber >= 203602 && RunNumber <= 203792 ) trigger_weight = 92.083 ; // for B5-6
    if (RunNumber >= 203875 && RunNumber <= 204158 ) trigger_weight = 123.825 ; // for B7-9
    if (RunNumber >= 204240 && RunNumber <= 204668 ) trigger_weight = 97.207 ; // for B10-11
    if (RunNumber >= 204707 && RunNumber <= 205017 ) trigger_weight = 111.688 ; // for B12
    if (RunNumber >= 205055 && RunNumber <= 205113 ) trigger_weight = 107.994 ; // for B13-14
    if (RunNumber >= 206248 && RunNumber <= 206614 ) trigger_weight = 86.390 ; // for C1-3
    if (RunNumber >= 206717 && RunNumber <= 207046 ) trigger_weight = 108.723 ; // for C4-C6
    if (RunNumber >= 207113 && RunNumber <= 207397 ) trigger_weight = 96.007 ; // for C7-9
    if (RunNumber >= 207447 && RunNumber <= 208126 ) trigger_weight = 99.480 ; // for D1-3
    if (RunNumber >= 208179 && RunNumber <= 208485 ) trigger_weight = 113.072 ; // for D4-6
    if (RunNumber >= 208631 && RunNumber <= 209025 ) trigger_weight = 99.623 ; // for D7-8
    if (RunNumber >= 209074 && RunNumber <= 209899 ) trigger_weight = 113.506 ; // for E1-3
    if (RunNumber >= 209909 && RunNumber <= 210308 ) trigger_weight = 99.280 ; // for E4-5
    if (RunNumber >= 211522 && RunNumber <= 211620 ) trigger_weight = 80.359 ; // for G1-G2
    if (RunNumber >= 211670 && RunNumber <= 212142 ) trigger_weight = 109.099 ; // for G3-4
    if (RunNumber >= 212144 && RunNumber <= 212272 ) trigger_weight = 104.649 ; // for G5
    if (RunNumber >= 212619 && RunNumber <= 212858 ) trigger_weight = 112.230 ; // for H1-2
    if (RunNumber >= 212967 && RunNumber <= 213250 ) trigger_weight = 114.354 ; // for H3-4
    if (RunNumber >= 213264 && RunNumber <= 213359 ) trigger_weight = 159.317 ; // for H5-6
    if (RunNumber >= 213431 && RunNumber <= 213819 ) trigger_weight = 117.052 ; // for I1-3
    if (RunNumber >= 213900 && RunNumber <= 214216 ) trigger_weight = 108.264 ; // for J1-2
    if (RunNumber >= 214388 && RunNumber <= 214777 ) trigger_weight = 105.725 ; // for J3-4
    if (RunNumber >= 214984 && RunNumber <= 215061 ) trigger_weight = 94.889 ; // for J5-6
    if (RunNumber >= 215063 && RunNumber <= 215091 ) trigger_weight = 104.521 ; // for J7-8
    if (RunNumber >= 215414 && RunNumber <= 215541 ) trigger_weight = 102.582 ; // for L1
    if (RunNumber >= 215559 && RunNumber <= 215643 ) trigger_weight = 130.212 ; // for L2-3
    //------------------------------------------------------------
  }
  /// one trigger end
  ///EF 180 trigger edit
  if (EF_j180_a4tchad && nJets && leadingJetPt > 230*GeV && leadingJetPt <= 270*GeV){
    passTrigger = true;
    //Prescales for trigger EF_j180_a4tchad
    //------------------------------------------------------------
    if (RunNumber >= 200804 && RunNumber <= 200913 ) trigger_weight = 1.000 ; // for A1-3
    if (RunNumber >= 200926 && RunNumber <= 201191 ) trigger_weight = 83.349 ; // for A4
    if (RunNumber >= 201257 && RunNumber <= 201556 ) trigger_weight = 233.287 ; // for A5-8
    if (RunNumber >= 202660 && RunNumber <= 203195 ) trigger_weight = 168.645 ; // for B1-3
    if (RunNumber >= 203228 && RunNumber <= 203524 ) trigger_weight = 240.955 ; // for B4
    if (RunNumber >= 203602 && RunNumber <= 203792 ) trigger_weight = 256.119 ; // for B5-6
    if (RunNumber >= 203875 && RunNumber <= 204158 ) trigger_weight = 344.715 ; // for B7-9
    if (RunNumber >= 204240 && RunNumber <= 204668 ) trigger_weight = 270.725 ; // for B10-11
    if (RunNumber >= 204707 && RunNumber <= 205017 ) trigger_weight = 311.043 ; // for B12
    if (RunNumber >= 205055 && RunNumber <= 205113 ) trigger_weight = 300.761 ; // for B13-14
    if (RunNumber >= 206248 && RunNumber <= 206614 ) trigger_weight = 240.655 ; // for C1-3
    if (RunNumber >= 206717 && RunNumber <= 207046 ) trigger_weight = 302.806 ; // for C4-C6
    if (RunNumber >= 207113 && RunNumber <= 207397 ) trigger_weight = 267.444 ; // for C7-9
    if (RunNumber >= 207447 && RunNumber <= 208126 ) trigger_weight = 277.090 ; // for D1-3
    if (RunNumber >= 208179 && RunNumber <= 208485 ) trigger_weight = 314.843 ; // for D4-6
    if (RunNumber >= 208631 && RunNumber <= 209025 ) trigger_weight = 277.475 ; // for D7-8
    if (RunNumber >= 209074 && RunNumber <= 209899 ) trigger_weight = 316.165 ; // for E1-3
    if (RunNumber >= 209909 && RunNumber <= 210308 ) trigger_weight = 276.583 ; // for E4-5
    if (RunNumber >= 211522 && RunNumber <= 211620 ) trigger_weight = 223.980 ; // for G1-G2
    if (RunNumber >= 211670 && RunNumber <= 212142 ) trigger_weight = 303.957 ; // for G3-4
    if (RunNumber >= 212144 && RunNumber <= 212272 ) trigger_weight = 291.515 ; // for G5
    if (RunNumber >= 212619 && RunNumber <= 212858 ) trigger_weight = 312.102 ; // for H1-2
    if (RunNumber >= 212967 && RunNumber <= 213250 ) trigger_weight = 317.660 ; // for H3-4
    if (RunNumber >= 213264 && RunNumber <= 213359 ) trigger_weight = 443.254 ; // for H5-6
    if (RunNumber >= 213431 && RunNumber <= 213819 ) trigger_weight = 325.099 ; // for I1-3
    if (RunNumber >= 213900 && RunNumber <= 214216 ) trigger_weight = 300.384 ; // for J1-2
    if (RunNumber >= 214388 && RunNumber <= 214777 ) trigger_weight = 292.109 ; // for J3-4
    if (RunNumber >= 214984 && RunNumber <= 215061 ) trigger_weight = 262.138 ; // for J5-6
    if (RunNumber >= 215063 && RunNumber <= 215091 ) trigger_weight = 288.758 ; // for J7-8
    if (RunNumber >= 215414 && RunNumber <= 215541 ) trigger_weight = 283.433 ; // for L1
    if (RunNumber >= 215559 && RunNumber <= 215643 ) trigger_weight = 359.636 ; // for L2-3
    //------------------------------------------------------------
  }
 
  //// One trigger code block two comment
 
  if (EF_j110_a4tchad && nJets && leadingJetPt > 160*GeV && leadingJetPt <= 230*GeV){
    passTrigger = true;
    //Prescales for trigger EF_j110_a4tchad
    //------------------------------------------------------------
    if (RunNumber >= 200804 && RunNumber <= 200913 ) trigger_weight = 8.415 ; // for A1-3
    if (RunNumber >= 200926 && RunNumber <= 201191 ) trigger_weight = 217.137 ; // for A4
    if (RunNumber >= 201257 && RunNumber <= 201556 ) trigger_weight = 701.063 ; // for A5-8
    if (RunNumber >= 202660 && RunNumber <= 203195 ) trigger_weight = 1628.988 ; // for B1-3
    if (RunNumber >= 203228 && RunNumber <= 203524 ) trigger_weight = 2302.571 ; // for B4
    if (RunNumber >= 203602 && RunNumber <= 203792 ) trigger_weight = 2443.241 ; // for B5-6
    if (RunNumber >= 203875 && RunNumber <= 204158 ) trigger_weight = 3301.764 ; // for B7-9
    if (RunNumber >= 204240 && RunNumber <= 204668 ) trigger_weight = 2576.138 ; // for B10-11
    if (RunNumber >= 204707 && RunNumber <= 205017 ) trigger_weight = 2961.772 ; // for B12
    if (RunNumber >= 205055 && RunNumber <= 205113 ) trigger_weight = 2859.528 ; // for B13-14
    if (RunNumber >= 206248 && RunNumber <= 206614 ) trigger_weight = 2292.202 ; // for C1-3
    if (RunNumber >= 206717 && RunNumber <= 207046 ) trigger_weight = 2883.170 ; // for C4-C6
    if (RunNumber >= 207113 && RunNumber <= 207397 ) trigger_weight = 2546.448 ; // for C7-9
    if (RunNumber >= 207447 && RunNumber <= 208126 ) trigger_weight = 2640.031 ; // for D1-3
    if (RunNumber >= 208179 && RunNumber <= 208485 ) trigger_weight = 2997.937 ; // for D4-6
    if (RunNumber >= 208631 && RunNumber <= 209025 ) trigger_weight = 2644.615 ; // for D7-8
    if (RunNumber >= 209074 && RunNumber <= 209899 ) trigger_weight = 3006.812 ; // for E1-3
    if (RunNumber >= 209909 && RunNumber <= 210308 ) trigger_weight = 2638.303 ; // for E4-5
    if (RunNumber >= 211522 && RunNumber <= 211620 ) trigger_weight = 2121.480 ; // for G1-G2
    if (RunNumber >= 211670 && RunNumber <= 212142 ) trigger_weight = 2879.332 ; // for G3-4
    if (RunNumber >= 212144 && RunNumber <= 212272 ) trigger_weight = 2761.621 ; // for G5
    if (RunNumber >= 212619 && RunNumber <= 212858 ) trigger_weight = 2954.885 ; // for H1-2
    if (RunNumber >= 212967 && RunNumber <= 213250 ) trigger_weight = 3000.150 ; // for H3-4
    if (RunNumber >= 213264 && RunNumber <= 213359 ) trigger_weight = 4183.325 ; // for H5-6
    if (RunNumber >= 213431 && RunNumber <= 213819 ) trigger_weight = 3073.137 ; // for I1-3
    if (RunNumber >= 213900 && RunNumber <= 214216 ) trigger_weight = 2871.040 ; // for J1-2
    if (RunNumber >= 214388 && RunNumber <= 214777 ) trigger_weight = 2931.223 ; // for J3-4
    if (RunNumber >= 214984 && RunNumber <= 215061 ) trigger_weight = 2627.258 ; // for J5-6
    if (RunNumber >= 215063 && RunNumber <= 215091 ) trigger_weight = 2902.027 ; // for J7-8
    if (RunNumber >= 215414 && RunNumber <= 215541 ) trigger_weight = 2853.598 ; // for L1
    if (RunNumber >= 215559 && RunNumber <= 215643 ) trigger_weight = 3595.281 ; // for L2-3
    //------------------------------------------------------------
  }

  if (EF_j80_a4tchad && nJets && leadingJetPt > 130*GeV && leadingJetPt <= 160*GeV){
    passTrigger = true;
    //Prescales for trigger EF_j80_a4tchad
    //------------------------------------------------------------
    if (RunNumber >= 200804 && RunNumber <= 200913 ) trigger_weight = 36.617 ; // for A1-3
    if (RunNumber >= 200926 && RunNumber <= 201191 ) trigger_weight = 955.981 ; // for A4
    if (RunNumber >= 201257 && RunNumber <= 201556 ) trigger_weight = 3130.252 ; // for A5-8
    if (RunNumber >= 202660 && RunNumber <= 203195 ) trigger_weight = 7108.769 ; // for B1-3
    if (RunNumber >= 203228 && RunNumber <= 203524 ) trigger_weight = 9573.527 ; // for B4
    if (RunNumber >= 203602 && RunNumber <= 203792 ) trigger_weight = 10154.120 ; // for B5-6
    if (RunNumber >= 203875 && RunNumber <= 204158 ) trigger_weight = 13696.955 ; // for B7-9
    if (RunNumber >= 204240 && RunNumber <= 204668 ) trigger_weight = 10709.599 ; // for B10-11
    if (RunNumber >= 204707 && RunNumber <= 205017 ) trigger_weight = 12242.402 ; // for B12
    if (RunNumber >= 205055 && RunNumber <= 205113 ) trigger_weight = 11894.048 ; // for B13-14
    if (RunNumber >= 206248 && RunNumber <= 206614 ) trigger_weight = 9520.698 ; // for C1-3
    if (RunNumber >= 206717 && RunNumber <= 207046 ) trigger_weight = 11988.900 ; // for C4-C6
    if (RunNumber >= 207113 && RunNumber <= 207397 ) trigger_weight = 10582.584 ; // for C7-9
    if (RunNumber >= 207447 && RunNumber <= 208126 ) trigger_weight = 10973.889 ; // for D1-3
    if (RunNumber >= 208179 && RunNumber <= 208485 ) trigger_weight = 12463.907 ; // for D4-6
    if (RunNumber >= 208631 && RunNumber <= 209025 ) trigger_weight = 10987.328 ; // for D7-8
    if (RunNumber >= 209074 && RunNumber <= 209899 ) trigger_weight = 12490.944 ; // for E1-3
    if (RunNumber >= 209909 && RunNumber <= 210308 ) trigger_weight = 10969.140 ; // for E4-5
    if (RunNumber >= 211522 && RunNumber <= 211620 ) trigger_weight = 8878.174 ; // for G1-G2
    if (RunNumber >= 211670 && RunNumber <= 212142 ) trigger_weight = 12040.525 ; // for G3-4
    if (RunNumber >= 212144 && RunNumber <= 212272 ) trigger_weight = 11552.765 ; // for G5
    if (RunNumber >= 212619 && RunNumber <= 212858 ) trigger_weight = 12332.077 ; // for H1-2
    if (RunNumber >= 212967 && RunNumber <= 213250 ) trigger_weight = 12518.988 ; // for H3-4
    if (RunNumber >= 213264 && RunNumber <= 213359 ) trigger_weight = 17468.191 ; // for H5-6
    if (RunNumber >= 213431 && RunNumber <= 213819 ) trigger_weight = 12853.131 ; // for I1-3
    if (RunNumber >= 213900 && RunNumber <= 214216 ) trigger_weight = 12007.911 ; // for J1-2
    if (RunNumber >= 214388 && RunNumber <= 214777 ) trigger_weight = 12360.152 ; // for J3-4
    if (RunNumber >= 214984 && RunNumber <= 215061 ) trigger_weight = 11091.373 ; // for J5-6
    if (RunNumber >= 215063 && RunNumber <= 215091 ) trigger_weight = 12220.752 ; // for J7-8
    if (RunNumber >= 215414 && RunNumber <= 215541 ) trigger_weight = 11997.862 ; // for L1
    if (RunNumber >= 215559 && RunNumber <= 215643 ) trigger_weight = 15224.875 ; // for L2-3
    //------------------------------------------------------------
  }

  if (EF_j55_a4tchad && nJets && leadingJetPt > 90*GeV && leadingJetPt <= 130*GeV){
    passTrigger = true;
    //Prescales for trigger EF_j55_a4tchad
    //------------------------------------------------------------
    if (RunNumber >= 200804 && RunNumber <= 200913 ) trigger_weight = 195.325 ; // for A1-3
    if (RunNumber >= 200926 && RunNumber <= 201191 ) trigger_weight = 4966.323 ; // for A4
    if (RunNumber >= 201257 && RunNumber <= 201556 ) trigger_weight = 16167.480 ; // for A5-8
    if (RunNumber >= 202660 && RunNumber <= 203195 ) trigger_weight = 37758.218 ; // for B1-3
    if (RunNumber >= 203228 && RunNumber <= 203524 ) trigger_weight = 50458.688 ; // for B4
    if (RunNumber >= 203602 && RunNumber <= 203792 ) trigger_weight = 53489.679 ; // for B5-6
    if (RunNumber >= 203875 && RunNumber <= 204158 ) trigger_weight = 72078.432 ; // for B7-9
    if (RunNumber >= 204240 && RunNumber <= 204668 ) trigger_weight = 56358.448 ; // for B10-11
    if (RunNumber >= 204707 && RunNumber <= 205017 ) trigger_weight = 64584.200 ; // for B12
    if (RunNumber >= 205055 && RunNumber <= 205113 ) trigger_weight = 62610.418 ; // for B13-14
    if (RunNumber >= 206248 && RunNumber <= 206614 ) trigger_weight = 50149.534 ; // for C1-3
    if (RunNumber >= 206717 && RunNumber <= 207046 ) trigger_weight = 63110.710 ; // for C4-C6
    if (RunNumber >= 207113 && RunNumber <= 207397 ) trigger_weight = 55723.209 ; // for C7-9
    if (RunNumber >= 207447 && RunNumber <= 208126 ) trigger_weight = 57743.354 ; // for D1-3
    if (RunNumber >= 208179 && RunNumber <= 208485 ) trigger_weight = 65589.369 ; // for D4-6
    if (RunNumber >= 208631 && RunNumber <= 209025 ) trigger_weight = 57837.998 ; // for D7-8
    if (RunNumber >= 209074 && RunNumber <= 209899 ) trigger_weight = 65752.318 ; // for E1-3
    if (RunNumber >= 209909 && RunNumber <= 210308 ) trigger_weight = 57688.794 ; // for E4-5
    if (RunNumber >= 211522 && RunNumber <= 211620 ) trigger_weight = 46649.935 ; // for G1-G2
    if (RunNumber >= 211670 && RunNumber <= 212142 ) trigger_weight = 63066.979 ; // for G3-4
    if (RunNumber >= 212144 && RunNumber <= 212272 ) trigger_weight = 60538.954 ; // for G5
    if (RunNumber >= 212619 && RunNumber <= 212858 ) trigger_weight = 64976.136 ; // for H1-2
    if (RunNumber >= 212967 && RunNumber <= 213250 ) trigger_weight = 66170.805 ; // for H3-4
    if (RunNumber >= 213264 && RunNumber <= 213359 ) trigger_weight = 92350.130 ; // for H5-6
    if (RunNumber >= 213431 && RunNumber <= 213819 ) trigger_weight = 67785.950 ; // for I1-3
    if (RunNumber >= 213900 && RunNumber <= 214216 ) trigger_weight = 63110.536 ; // for J1-2
    if (RunNumber >= 214388 && RunNumber <= 214777 ) trigger_weight = 63562.289 ; // for J3-4
    if (RunNumber >= 214984 && RunNumber <= 215061 ) trigger_weight = 57039.299 ; // for J5-6
    if (RunNumber >= 215063 && RunNumber <= 215091 ) trigger_weight = 63017.587 ; // for J7-8
    if (RunNumber >= 215414 && RunNumber <= 215541 ) trigger_weight = 61965.647 ; // for L1
    if (RunNumber >= 215559 && RunNumber <= 215643 ) trigger_weight = 78668.614 ; // for L2-3
    //------------------------------------------------------------
  }
  
  return passTrigger;

}

bool QCDSeedSelection::passBJetTrigger_2012Data(int RunNumber,bool EF_b360, bool EF_b280, bool EF_b220, bool EF_b180,bool EF_b110, bool EF_b80, bool EF_b55 ,int nJets, float leadingJetPt, float& trigger_weight){

  bool passTrigger = false;

  //Bjet one trigger pt spectrum //  
  if (EF_b360 && nJets && leadingJetPt > 540*GeV){
    passTrigger = true;
    //Prescales for trigger EF_b360_loose_j360_a4tchad_L2j140
    //------------------------------------------------------------
    if (RunNumber >= 200804 && RunNumber <= 200913 ) trigger_weight = 1.000 ; // for A1-3
    if (RunNumber >= 200926 && RunNumber <= 201191 ) trigger_weight = 1.000 ; // for A4
    if (RunNumber >= 201257 && RunNumber <= 201556 ) trigger_weight = 1.000 ; // for A5-8
    if (RunNumber >= 202660 && RunNumber <= 203195 ) trigger_weight = 1.000 ; // for B1-3
    if (RunNumber >= 203228 && RunNumber <= 203524 ) trigger_weight = 1.000 ; // for B4
    if (RunNumber >= 203602 && RunNumber <= 203792 ) trigger_weight = 1.000 ; // for B5-6
    if (RunNumber >= 203875 && RunNumber <= 204158 ) trigger_weight = 1.000 ; // for B7-9
    if (RunNumber >= 204240 && RunNumber <= 204668 ) trigger_weight = 1.000 ; // for B10-11
    if (RunNumber >= 204707 && RunNumber <= 205017 ) trigger_weight = 1.003 ; // for B12
    if (RunNumber >= 205055 && RunNumber <= 205113 ) trigger_weight = 1.000 ; // for B13-14
    if (RunNumber >= 206248 && RunNumber <= 206614 ) trigger_weight = 1.000 ; // for C1-3
    if (RunNumber >= 206717 && RunNumber <= 207046 ) trigger_weight = 1.000 ; // for C4-C6
    if (RunNumber >= 207113 && RunNumber <= 207397 ) trigger_weight = 1.000 ; // for C7-9
    if (RunNumber >= 207447 && RunNumber <= 208126 ) trigger_weight = 1.000 ; // for D1-3
    if (RunNumber >= 208179 && RunNumber <= 208485 ) trigger_weight = 1.000 ; // for D4-6
    if (RunNumber >= 208631 && RunNumber <= 209025 ) trigger_weight = 1.000 ; // for D7-8
    if (RunNumber >= 209074 && RunNumber <= 209899 ) trigger_weight = 1.000 ; // for E1-3
    if (RunNumber >= 209909 && RunNumber <= 210308 ) trigger_weight = 1.000 ; // for E4-5
    if (RunNumber >= 211522 && RunNumber <= 211620 ) trigger_weight = 1.000 ; // for G1-G2
    if (RunNumber >= 211670 && RunNumber <= 212142 ) trigger_weight = 1.001 ; // for G3-4
    if (RunNumber >= 212144 && RunNumber <= 212272 ) trigger_weight = 1.000 ; // for G5
    if (RunNumber >= 212619 && RunNumber <= 212858 ) trigger_weight = 1.000 ; // for H1-2
    if (RunNumber >= 212967 && RunNumber <= 213250 ) trigger_weight = 1.000 ; // for H3-4
    if (RunNumber >= 213264 && RunNumber <= 213359 ) trigger_weight = 1.000 ; // for H5-6
    if (RunNumber >= 213431 && RunNumber <= 213819 ) trigger_weight = 1.000 ; // for I1-3
    if (RunNumber >= 213900 && RunNumber <= 214216 ) trigger_weight = 1.000 ; // for J1-2
    if (RunNumber >= 214388 && RunNumber <= 214777 ) trigger_weight = 1.000 ; // for J3-4
    if (RunNumber >= 214984 && RunNumber <= 215061 ) trigger_weight = 1.000 ; // for J5-6
    if (RunNumber >= 215063 && RunNumber <= 215091 ) trigger_weight = 1.000 ; // for J7-8
    if (RunNumber >= 215414 && RunNumber <= 215541 ) trigger_weight = 1.000 ; // for L1
    if (RunNumber >= 215559 && RunNumber <= 215643 ) trigger_weight = 1.000 ; // for L2-3
    //------------------------------------------------------------
  }

  if (EF_b280 && nJets && leadingJetPt > 420*GeV && leadingJetPt <= 540*GeV){
    passTrigger = true;

    //Prescales for trigger EF_b280_loose_j280_a4tchad_L2j140
    //------------------------------------------------------------
    if (RunNumber >= 200804 && RunNumber <= 200913 ) trigger_weight = 1.000 ; // for A1-3
    if (RunNumber >= 200926 && RunNumber <= 201191 ) trigger_weight = 3.847 ; // for A4
    if (RunNumber >= 201257 && RunNumber <= 201556 ) trigger_weight = 10.688 ; // for A5-8
    if (RunNumber >= 202660 && RunNumber <= 203195 ) trigger_weight = 8.579 ; // for B1-3
    if (RunNumber >= 203228 && RunNumber <= 203524 ) trigger_weight = 10.423 ; // for B4
    if (RunNumber >= 203602 && RunNumber <= 203792 ) trigger_weight = 11.526 ; // for B5-6
    if (RunNumber >= 203875 && RunNumber <= 204158 ) trigger_weight = 16.222 ; // for B7-9
    if (RunNumber >= 204240 && RunNumber <= 204668 ) trigger_weight = 12.756 ; // for B10-11
    if (RunNumber >= 204707 && RunNumber <= 205017 ) trigger_weight = 14.642 ; // for B12
    if (RunNumber >= 205055 && RunNumber <= 205113 ) trigger_weight = 14.161 ; // for B13-14
    if (RunNumber >= 206248 && RunNumber <= 206614 ) trigger_weight = 11.347 ; // for C1-3
    if (RunNumber >= 206717 && RunNumber <= 207046 ) trigger_weight = 14.257 ; // for C4-C6
    if (RunNumber >= 207113 && RunNumber <= 207397 ) trigger_weight = 12.597 ; // for C7-9
    if (RunNumber >= 207447 && RunNumber <= 208126 ) trigger_weight = 13.049 ; // for D1-3
    if (RunNumber >= 208179 && RunNumber <= 208485 ) trigger_weight = 14.820 ; // for D4-6
    if (RunNumber >= 208631 && RunNumber <= 209025 ) trigger_weight = 13.073 ; // for D7-8
    if (RunNumber >= 209074 && RunNumber <= 209899 ) trigger_weight = 14.883 ; // for E1-3
    if (RunNumber >= 209909 && RunNumber <= 210308 ) trigger_weight = 13.026 ; // for E4-5
    if (RunNumber >= 211522 && RunNumber <= 211620 ) trigger_weight = 10.561 ; // for G1-G2
    if (RunNumber >= 211670 && RunNumber <= 212142 ) trigger_weight = 14.308 ; // for G3-4
    if (RunNumber >= 212144 && RunNumber <= 212272 ) trigger_weight = 13.727 ; // for G5
    if (RunNumber >= 212619 && RunNumber <= 212858 ) trigger_weight = 14.894 ; // for H1-2
    if (RunNumber >= 212967 && RunNumber <= 213250 ) trigger_weight = 15.334 ; // for H3-4
    if (RunNumber >= 213264 && RunNumber <= 213359 ) trigger_weight = 21.395 ; // for H5-6
    if (RunNumber >= 213431 && RunNumber <= 213819 ) trigger_weight = 15.688 ; // for I1-3
    if (RunNumber >= 213900 && RunNumber <= 214216 ) trigger_weight = 15.328 ; // for J1-2
    if (RunNumber >= 214388 && RunNumber <= 214777 ) trigger_weight = 20.399 ; // for J3-4
    if (RunNumber >= 214984 && RunNumber <= 215061 ) trigger_weight = 18.317 ; // for J5-6
    if (RunNumber >= 215063 && RunNumber <= 215091 ) trigger_weight = 20.171 ; // for J7-8
    if (RunNumber >= 215414 && RunNumber <= 215541 ) trigger_weight = 19.802 ; // for L1
    if (RunNumber >= 215559 && RunNumber <= 215643 ) trigger_weight = 25.106 ; // for L2-3
    //------------------------------------------------------------
  }
  
  // two trigger end point
  if (EF_b220 && nJets && leadingJetPt > 330*GeV  && leadingJetPt <= 420*GeV){
    passTrigger = true;
    
    //Prescales for trigger EF_b220_loose_j220_a4tchad_L2j140
    //------------------------------------------------------------
    if (RunNumber >= 200804 && RunNumber <= 200913 ) trigger_weight = 1.185 ; // for A1-3
    if (RunNumber >= 200926 && RunNumber <= 201191 ) trigger_weight = 11.831 ; // for A4
    if (RunNumber >= 201257 && RunNumber <= 201556 ) trigger_weight = 29.064 ; // for A5-8
    if (RunNumber >= 202660 && RunNumber <= 203195 ) trigger_weight = 29.720 ; // for B1-3
    if (RunNumber >= 203228 && RunNumber <= 203524 ) trigger_weight = 49.664 ; // for B4
    if (RunNumber >= 203602 && RunNumber <= 203792 ) trigger_weight = 53.162 ; // for B5-6
    if (RunNumber >= 203875 && RunNumber <= 204158 ) trigger_weight = 72.120 ; // for B7-9
    if (RunNumber >= 204240 && RunNumber <= 204668 ) trigger_weight = 56.629 ; // for B10-11
    if (RunNumber >= 204707 && RunNumber <= 205017 ) trigger_weight = 65.060 ; // for B12
    if (RunNumber >= 205055 && RunNumber <= 205113 ) trigger_weight = 62.896 ; // for B13-14
    if (RunNumber >= 206248 && RunNumber <= 206614 ) trigger_weight = 50.335 ; // for C1-3
    if (RunNumber >= 206717 && RunNumber <= 207046 ) trigger_weight = 63.341 ; // for C4-C6
    if (RunNumber >= 207113 && RunNumber <= 207397 ) trigger_weight = 55.930 ; // for C7-9
    if (RunNumber >= 207447 && RunNumber <= 208126 ) trigger_weight = 57.960 ; // for D1-3
    if (RunNumber >= 208179 && RunNumber <= 208485 ) trigger_weight = 65.868 ; // for D4-6
    if (RunNumber >= 208631 && RunNumber <= 209025 ) trigger_weight = 58.036 ; // for D7-8
    if (RunNumber >= 209074 && RunNumber <= 209899 ) trigger_weight = 66.136 ; // for E1-3
    if (RunNumber >= 209909 && RunNumber <= 210308 ) trigger_weight = 57.847 ; // for E4-5
    if (RunNumber >= 211522 && RunNumber <= 211620 ) trigger_weight = 45.052 ; // for G1-G2
    if (RunNumber >= 211670 && RunNumber <= 212142 ) trigger_weight = 61.179 ; // for G3-4
    if (RunNumber >= 212144 && RunNumber <= 212272 ) trigger_weight = 58.672 ; // for G5
    if (RunNumber >= 212619 && RunNumber <= 212858 ) trigger_weight = 61.913 ; // for H1-2
    if (RunNumber >= 212967 && RunNumber <= 213250 ) trigger_weight = 62.218 ; // for H3-4
    if (RunNumber >= 213264 && RunNumber <= 213359 ) trigger_weight = 86.733 ; // for H5-6
    if (RunNumber >= 213431 && RunNumber <= 213819 ) trigger_weight = 63.668 ; // for I1-3
    if (RunNumber >= 213900 && RunNumber <= 214216 ) trigger_weight = 60.617 ; // for J1-2
    if (RunNumber >= 214388 && RunNumber <= 214777 ) trigger_weight = 68.635 ; // for J3-4
    if (RunNumber >= 214984 && RunNumber <= 215061 ) trigger_weight = 61.595 ; // for J5-6
    if (RunNumber >= 215063 && RunNumber <= 215091 ) trigger_weight = 67.858 ; // for J7-8
    if (RunNumber >= 215414 && RunNumber <= 215541 ) trigger_weight = 66.611 ; // for L1
    if (RunNumber >= 215559 && RunNumber <= 215643 ) trigger_weight = 84.562 ; // for L2-3
    //------------------------------------------------------------
  }
  /// End b-jet leading trigger spectrum

  if (EF_b180 && nJets && leadingJetPt > 270*GeV && leadingJetPt <= 330*GeV){
    passTrigger = true;
    //Prescales for trigger EF_b180_loose_j180_a4tchad_L2j140
    //------------------------------------------------------------
    if (RunNumber >= 200804 && RunNumber <= 200913 ) trigger_weight = 2.389 ; // for A1-3
    if (RunNumber >= 200926 && RunNumber <= 201191 ) trigger_weight = 36.926 ; // for A4
    if (RunNumber >= 201257 && RunNumber <= 201556 ) trigger_weight = 89.581 ; // for A5-8
    if (RunNumber >= 202660 && RunNumber <= 203195 ) trigger_weight = 90.058 ; // for B1-3
    if (RunNumber >= 203228 && RunNumber <= 203524 ) trigger_weight = 120.593 ; // for B4
    if (RunNumber >= 203602 && RunNumber <= 203792 ) trigger_weight = 128.565 ; // for B5-6
    if (RunNumber >= 203875 && RunNumber <= 204158 ) trigger_weight = 173.675 ; // for B7-9
    if (RunNumber >= 204240 && RunNumber <= 204668 ) trigger_weight = 136.359 ; // for B10-11
    if (RunNumber >= 204707 && RunNumber <= 205017 ) trigger_weight = 156.677 ; // for B12
    if (RunNumber >= 205055 && RunNumber <= 205113 ) trigger_weight = 151.466 ; // for B13-14
    if (RunNumber >= 206248 && RunNumber <= 206614 ) trigger_weight = 121.213 ; // for C1-3
    if (RunNumber >= 206717 && RunNumber <= 207046 ) trigger_weight = 151.001 ; // for C4-C6
    if (RunNumber >= 207113 && RunNumber <= 207397 ) trigger_weight = 133.293 ; // for C7-9
    if (RunNumber >= 207447 && RunNumber <= 208126 ) trigger_weight = 138.118 ; // for D1-3
    if (RunNumber >= 208179 && RunNumber <= 208485 ) trigger_weight = 156.963 ; // for D4-6
    if (RunNumber >= 208631 && RunNumber <= 209025 ) trigger_weight = 138.356 ; // for D7-8
    if (RunNumber >= 209074 && RunNumber <= 209899 ) trigger_weight = 157.583 ; // for E1-3
    if (RunNumber >= 209909 && RunNumber <= 210308 ) trigger_weight = 137.857 ; // for E4-5
    if (RunNumber >= 211522 && RunNumber <= 211620 ) trigger_weight = 111.797 ; // for G1-G2
    if (RunNumber >= 211670 && RunNumber <= 212142 ) trigger_weight = 151.492 ; // for G3-4
    if (RunNumber >= 212144 && RunNumber <= 212272 ) trigger_weight = 145.320 ; // for G5
    if (RunNumber >= 212619 && RunNumber <= 212858 ) trigger_weight = 154.837 ; // for H1-2
    if (RunNumber >= 212967 && RunNumber <= 213250 ) trigger_weight = 156.757 ; // for H3-4
    if (RunNumber >= 213264 && RunNumber <= 213359 ) trigger_weight = 218.449 ; // for H5-6
    if (RunNumber >= 213431 && RunNumber <= 213819 ) trigger_weight = 160.443 ; // for I1-3
    if (RunNumber >= 213900 && RunNumber <= 214216 ) trigger_weight = 148.818 ; // for J1-2
    if (RunNumber >= 214388 && RunNumber <= 214777 ) trigger_weight = 146.898 ; // for J3-4
    if (RunNumber >= 214984 && RunNumber <= 215061 ) trigger_weight = 131.830 ; // for J5-6
    if (RunNumber >= 215063 && RunNumber <= 215091 ) trigger_weight = 145.220 ; // for J7-8
    if (RunNumber >= 215414 && RunNumber <= 215541 ) trigger_weight = 142.614 ; // for L1
    if (RunNumber >= 215559 && RunNumber <= 215643 ) trigger_weight = 180.812 ; // for L2-3
    //------------------------------------------------------------
  }

  //bjet one trig block 2 
    
  if (EF_b110 && nJets && leadingJetPt > 165*GeV && leadingJetPt <= 270*GeV /* leadingJetPt > 140*GeV && leadingJetPt <= 220*GeV*/){      
    passTrigger = true;
    
    //Prescales for trigger EF_b110_looseEF_j110_a4tchad
    //------------------------------------------------------------
    if (RunNumber >= 200804 && RunNumber <= 200913 ) trigger_weight = 12.622 ; // for A1-3
    if (RunNumber >= 200926 && RunNumber <= 201191 ) trigger_weight = 177.828 ; // for A4
    if (RunNumber >= 201257 && RunNumber <= 201556 ) trigger_weight = 324.282 ; // for A5-8
    if (RunNumber >= 202660 && RunNumber <= 203195 ) trigger_weight = 240.250 ; // for B1-3
    if (RunNumber >= 203228 && RunNumber <= 203524 ) trigger_weight = 359.775 ; // for B4
    if (RunNumber >= 203602 && RunNumber <= 203792 ) trigger_weight = 382.157 ; // for B5-6
    if (RunNumber >= 203875 && RunNumber <= 204158 ) trigger_weight = 516.209 ; // for B7-9
    if (RunNumber >= 204240 && RunNumber <= 204668 ) trigger_weight = 403.799 ; // for B10-11
    if (RunNumber >= 204707 && RunNumber <= 205017 ) trigger_weight = 463.977 ; // for B12
    if (RunNumber >= 205055 && RunNumber <= 205113 ) trigger_weight = 448.591 ; // for B13-14
    if (RunNumber >= 206248 && RunNumber <= 206614 ) trigger_weight = 358.158 ; // for C1-3
    if (RunNumber >= 206717 && RunNumber <= 207046 ) trigger_weight = 450.494 ; // for C4-C6
    if (RunNumber >= 207113 && RunNumber <= 207397 ) trigger_weight = 397.881 ; // for C7-9
    if (RunNumber >= 207447 && RunNumber <= 208126 ) trigger_weight = 412.526 ; // for D1-3
    if (RunNumber >= 208179 && RunNumber <= 208485 ) trigger_weight = 468.430 ; // for D4-6
    if (RunNumber >= 208631 && RunNumber <= 209025 ) trigger_weight = 413.220 ; // for D7-8
    if (RunNumber >= 209074 && RunNumber <= 209899 ) trigger_weight = 469.814 ; // for E1-3
    if (RunNumber >= 209909 && RunNumber <= 210308 ) trigger_weight = 412.235 ; // for E4-5
    if (RunNumber >= 211522 && RunNumber <= 211620 ) trigger_weight = 336.743 ; // for G1-G2
    if (RunNumber >= 211670 && RunNumber <= 212142 ) trigger_weight = 457.035 ; // for G3-4
    if (RunNumber >= 212144 && RunNumber <= 212272 ) trigger_weight = 438.354 ; // for G5
    if (RunNumber >= 212619 && RunNumber <= 212858 ) trigger_weight = 469.031 ; // for H1-2
    if (RunNumber >= 212967 && RunNumber <= 213250 ) trigger_weight = 476.214 ; // for H3-4
    if (RunNumber >= 213264 && RunNumber <= 213359 ) trigger_weight = 664.021 ; // for H5-6
    if (RunNumber >= 213431 && RunNumber <= 213819 ) trigger_weight = 487.801 ; // for I1-3
    if (RunNumber >= 213900 && RunNumber <= 214216 ) trigger_weight = 452.492 ; // for J1-2
    if (RunNumber >= 214388 && RunNumber <= 214777 ) trigger_weight = 446.896 ; // for J3-4
    if (RunNumber >= 214984 && RunNumber <= 215061 ) trigger_weight = 401.119 ; // for J5-6
    if (RunNumber >= 215063 && RunNumber <= 215091 ) trigger_weight = 441.977 ; // for J7-8
    if (RunNumber >= 215414 && RunNumber <= 215541 ) trigger_weight = 433.993 ; // for L1
    if (RunNumber >= 215559 && RunNumber <= 215643 ) trigger_weight = 550.816 ; // for L2-3
    //------------------------------------------------------------
  }
  
  
  if (EF_b80 && nJets && leadingJetPt > 120*GeV && leadingJetPt <= 165*GeV){ 
    passTrigger = true;
    //Prescales for trigger EF_b80_looseEF_j80_a4tchad
    //------------------------------------------------------------
    if (RunNumber >= 200804 && RunNumber <= 200913 ) trigger_weight = 47.601 ; // for A1-3
    if (RunNumber >= 200926 && RunNumber <= 201191 ) trigger_weight = 703.209 ; // for A4
    if (RunNumber >= 201257 && RunNumber <= 201556 ) trigger_weight = 1349.511 ; // for A5-8
    if (RunNumber >= 202660 && RunNumber <= 203195 ) trigger_weight = 997.592 ; // for B1-3
    if (RunNumber >= 203228 && RunNumber <= 203524 ) trigger_weight = 1227.374 ; // for B4
    if (RunNumber >= 203602 && RunNumber <= 203792 ) trigger_weight = 1300.200 ; // for B5-6
    if (RunNumber >= 203875 && RunNumber <= 204158 ) trigger_weight = 1751.187 ; // for B7-9
    if (RunNumber >= 204240 && RunNumber <= 204668 ) trigger_weight = 1369.443 ; // for B10-11
    if (RunNumber >= 204707 && RunNumber <= 205017 ) trigger_weight = 1564.952 ; // for B12
    if (RunNumber >= 205055 && RunNumber <= 205113 ) trigger_weight = 1521.390 ; // for B13-14
    if (RunNumber >= 206248 && RunNumber <= 206614 ) trigger_weight = 1217.087 ; // for C1-3
    if (RunNumber >= 206717 && RunNumber <= 207046 ) trigger_weight = 1531.084 ; // for C4-C6
    if (RunNumber >= 207113 && RunNumber <= 207397 ) trigger_weight = 1352.277 ; // for C7-9
    if (RunNumber >= 207447 && RunNumber <= 208126 ) trigger_weight = 1402.124 ; // for D1-3
    if (RunNumber >= 208179 && RunNumber <= 208485 ) trigger_weight = 1592.454 ; // for D4-6
    if (RunNumber >= 208631 && RunNumber <= 209025 ) trigger_weight = 1404.368 ; // for D7-8
    if (RunNumber >= 209074 && RunNumber <= 209899 ) trigger_weight = 1597.426 ; // for E1-3
    if (RunNumber >= 209909 && RunNumber <= 210308 ) trigger_weight = 1401.056 ; // for E4-5
    if (RunNumber >= 211522 && RunNumber <= 211620 ) trigger_weight = 1117.962 ; // for G1-G2
    if (RunNumber >= 211670 && RunNumber <= 212142 ) trigger_weight = 1514.925 ; // for G3-4
    if (RunNumber >= 212144 && RunNumber <= 212272 ) trigger_weight = 1453.257 ; // for G5
    if (RunNumber >= 212619 && RunNumber <= 212858 ) trigger_weight = 1536.578 ; // for H1-2
    if (RunNumber >= 212967 && RunNumber <= 213250 ) trigger_weight = 1545.231 ; // for H3-4
    if (RunNumber >= 213264 && RunNumber <= 213359 ) trigger_weight = 2156.387 ; // for H5-6
    if (RunNumber >= 213431 && RunNumber <= 213819 ) trigger_weight = 1583.293 ; // for I1-3
    if (RunNumber >= 213900 && RunNumber <= 214216 ) trigger_weight = 1496.035 ; // for J1-2
    if (RunNumber >= 214388 && RunNumber <= 214777 ) trigger_weight = 1626.332 ; // for J3-4
    if (RunNumber >= 214984 && RunNumber <= 215061 ) trigger_weight = 1459.385 ; // for J5-6
    if (RunNumber >= 215063 && RunNumber <= 215091 ) trigger_weight = 1607.994 ; // for J7-8
    if (RunNumber >= 215414 && RunNumber <= 215541 ) trigger_weight = 1578.664 ; // for L1
    if (RunNumber >= 215559 && RunNumber <= 215643 ) trigger_weight = 2003.279 ; // for L2-3
    //------------------------------------------------------------
  }

  if (EF_b55 && nJets && leadingJetPt > 90*GeV && leadingJetPt <= 120*GeV){   
    passTrigger = true;
    //Prescales for trigger EF_b55_looseEF_j55_a4tchad
    //------------------------------------------------------------
    if (RunNumber >= 200804 && RunNumber <= 200913 ) trigger_weight = 195.325 ; // for A1-3
    if (RunNumber >= 200926 && RunNumber <= 201191 ) trigger_weight = 2796.567 ; // for A4
    if (RunNumber >= 201257 && RunNumber <= 201556 ) trigger_weight = 4965.912 ; // for A5-8
    if (RunNumber >= 202660 && RunNumber <= 203195 ) trigger_weight = 3661.927 ; // for B1-3
    if (RunNumber >= 203228 && RunNumber <= 203524 ) trigger_weight = 5203.868 ; // for B4
    if (RunNumber >= 203602 && RunNumber <= 203792 ) trigger_weight = 5514.405 ; // for B5-6
    if (RunNumber >= 203875 && RunNumber <= 204158 ) trigger_weight = 7430.750 ; // for B7-9
    if (RunNumber >= 204240 && RunNumber <= 204668 ) trigger_weight = 5810.141 ; // for B10-11
    if (RunNumber >= 204707 && RunNumber <= 205017 ) trigger_weight = 6658.153 ; // for B12
    if (RunNumber >= 205055 && RunNumber <= 205113 ) trigger_weight = 6454.666 ; // for B13-14
    if (RunNumber >= 206248 && RunNumber <= 206614 ) trigger_weight = 5170.053 ; // for C1-3
    if (RunNumber >= 206717 && RunNumber <= 207046 ) trigger_weight = 6506.264 ; // for C4-C6
    if (RunNumber >= 207113 && RunNumber <= 207397 ) trigger_weight = 5744.997 ; // for C7-9
    if (RunNumber >= 207447 && RunNumber <= 208126 ) trigger_weight = 5953.264 ; // for D1-3
    if (RunNumber >= 208179 && RunNumber <= 208485 ) trigger_weight = 6761.778 ; // for D4-6
    if (RunNumber >= 208631 && RunNumber <= 209025 ) trigger_weight = 5962.907 ; // for D7-8
    if (RunNumber >= 209074 && RunNumber <= 209899 ) trigger_weight = 6779.720 ; // for E1-3
    if (RunNumber >= 209909 && RunNumber <= 210308 ) trigger_weight = 5947.896 ; // for E4-5
    if (RunNumber >= 211522 && RunNumber <= 211620 ) trigger_weight = 4780.590 ; // for G1-G2
    if (RunNumber >= 211670 && RunNumber <= 212142 ) trigger_weight = 6488.883 ; // for G3-4
    if (RunNumber >= 212144 && RunNumber <= 212272 ) trigger_weight = 6223.532 ; // for G5
    if (RunNumber >= 212619 && RunNumber <= 212858 ) trigger_weight = 6437.873 ; // for H1-2
    if (RunNumber >= 212967 && RunNumber <= 213250 ) trigger_weight = 6362.573 ; // for H3-4
    if (RunNumber >= 213264 && RunNumber <= 213359 ) trigger_weight = 8879.827 ; // for H5-6
    if (RunNumber >= 213431 && RunNumber <= 213819 ) trigger_weight = 6517.871 ; // for I1-3
    if (RunNumber >= 213900 && RunNumber <= 214216 ) trigger_weight = 6078.157 ; // for J1-2
    if (RunNumber >= 214388 && RunNumber <= 214777 ) trigger_weight = 6170.317 ; // for J3-4
    if (RunNumber >= 214984 && RunNumber <= 215061 ) trigger_weight = 5537.236 ; // for J5-6
    if (RunNumber >= 215063 && RunNumber <= 215091 ) trigger_weight = 6081.201 ; // for J7-8
    if (RunNumber >= 215414 && RunNumber <= 215541 ) trigger_weight = 5958.221 ; // for L1
    if (RunNumber >= 215559 && RunNumber <= 215643 ) trigger_weight = 7564.308 ; // for L2-3
    //------------------------------------------------------------
  }

 
  return passTrigger;
}

bool QCDSeedSelection::pass2B35JetTrigger_2012Data(int RunNumber, bool EF_2b35_j145_j35, bool EF_2b35_j145_j100 ,bool EF_2b35_j110_2j35 ,bool EF_2b35_3j35, int nJets, float leadingJetPt, float& trigger_weight){
  
  bool passTrigger=false;
  
  if ( EF_2b35_3j35 && leadingJetPt > 90*GeV  && leadingJetPt <= 165*GeV ){
    passTrigger=true;
    trigger_weight=1.0;
  }

  else if ( (EF_2b35_j110_2j35 ) && leadingJetPt > 165*GeV  && leadingJetPt <= 210*GeV ){
    passTrigger=true;
    trigger_weight=1.0;
  }
   
  else if ( (EF_2b35_j145_j35 || EF_2b35_j145_j100 ) && leadingJetPt > 210*GeV ){
    passTrigger=true;
    trigger_weight=1.0;
    
  }
  
  return passTrigger;

}



bool QCDSeedSelection::selectSeedEvent(double MissingET, double sumET, double etMissSigCut){

  const double etMissSig = MissingET / sqrt(sumET);
  if (etMissSig < etMissSigCut) return true;
  else return false;

}





