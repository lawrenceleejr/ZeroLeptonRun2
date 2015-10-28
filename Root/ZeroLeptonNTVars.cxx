#include "ZeroLeptonRun2/ZeroLeptonNTVars.h"
#include "TTree.h"

std::string NTVars::toString()
{ 
  return
 
std::string("RunNumber/i:EventNumber/i:LumiBlockNumber/i:veto/i:eventWeight/F:pileupWeight/F:pileupWeightUp/F:pileupWeightDown/F:genWeight/F:ttbarWeightHT/F:ttbarWeightPt2/F:ttbarAvgPt/F:WZweight/F:nJet/i:met/F:metPhi/F:dPhi/F:dPhiR/F:meffInc/F:hardproc/I:nBJet/I:nCJet/I:bTagWeight/F:bTagWeightBUp/F:bTagWeightBDown/F:bTagWeightCUp/F:bTagWeightCDown/F:bTagWeightLUp/F:bTagWeightLDown/F:cTagWeight/F:cTagWeightBUp/F:cTagWeightBDown/F:cTagWeightCUp/F:cTagWeightCDown/F:cTagWeightLUp/F:cTagWeightLDown/F:normWeight/F:normWeightUp/F:normWeightDown/F:cleaning/i:timing/F:jet1Emf/F:jet2Emf/F:jet1Chf/F:jet2Chf/F:pdfId1/I:pdfId2/I:tauN/i:tauJetBDTLoose/i:tauLooseN/i:tauMt/F:SherpaBugMET/F:dPhiBadTile/F"); 
}

void NTVars::Reset()
{
  RunNumber = 0;
  EventNumber = 0;
  LumiBlockNumber = 0;
  veto = 0;
  weight = 1.f;
  pileupWeight = pileupWeightUp = pileupWeightDown = 1.f;
  genWeight = 1.f;
  ttbarWeightHT = 1.f;
  ttbarWeightPt2 = 1.f;
  ttbarAvgPt = 0.f;
  WZweight = 1.f;
  Njet = 0;
  MET = 0.f;
  METPhi = 0.f;
  deltaPhi = 0.f;
  deltaPhiRemaining = 0.f;
  MeffIncl=0.f;
  hardproc=-1;
  nBJet = 0;
  nCJet = 0;
  bTagWeight      = 1.f;
  bTagWeightBUp   = 1.f;
  bTagWeightBDown = 1.f;
  bTagWeightCUp   = 1.f;
  bTagWeightCDown = 1.f;
  bTagWeightLUp   = 1.f;
  bTagWeightLDown = 1.f;
  cTagWeight      = 1.f;
  cTagWeightBUp   = 1.f;
  cTagWeightBDown = 1.f;
  cTagWeightCUp   = 1.f;
  cTagWeightCDown = 1.f;
  cTagWeightLUp   = 1.f;
  cTagWeightLDown = 1.f;
  normWeight=-1.f;
  normWeightDown=-1.f;
  normWeightUp=-1.f;
  cleaning=0;
  timing=0;
  emfjet0=0;
  emfjet1=0;
  chfjet0=0;
  chfjet1=0;
  pdfId1 = pdfId2 = 0;
  tauN = 0;
  tauJetBDTLoose = 0;
  tauLooseN = 0;
  tauMt = 0.f;
  SherpaBugMET = 0.f;
  dPhiBadTile = 0.f;

  // Clear vectors
  jetPt.clear();
  jetEta.clear();
  jetPhi.clear();
  jetM.clear();
  jetBTag.clear();
  jetFlav.clear();
  jetTagU.clear();
  jetTagB.clear();
  jetTagC.clear();
  jetSmearSystW.clear();
  jetFracSamplingMax.clear();
  jetFracSamplingMaxIndex.clear();
  
  tauPt.clear();
  tauEta.clear();
  tauPhi.clear();
  tauLooseSF.clear();
  tauLooseSFStatUp.clear();
  tauLooseSFStatDown.clear();
  tauLooseSFSystUp.clear();
  tauLooseSFSystDown.clear();

  systWeights.clear();
}

NTVarsRead::NTVarsRead(): ntv()
{
  p_jetPt   = &ntv.jetPt;
  p_jetEta  = &ntv.jetEta;
  p_jetPhi  = &ntv.jetPhi;
  p_jetM    = &ntv.jetM;
  p_jetBTag = &ntv.jetBTag;
  p_jetFlav = &ntv.jetFlav;
  p_jetTagU = &ntv.jetTagU;
  p_jetTagB = &ntv.jetTagB;
  p_jetTagC = &ntv.jetTagC;
  p_jetSmearSystW = &ntv.jetSmearSystW;
  p_jetFracSamplingMax = &ntv.jetFracSamplingMax;
  p_jetFracSamplingMaxIndex = &ntv.jetFracSamplingMaxIndex;
  
  p_tauPt   = &ntv.tauPt;
  p_tauEta  = &ntv.tauEta;
  p_tauPhi  = &ntv.tauPhi;
  p_tauLooseSF         = &ntv.tauLooseSF;
  p_tauLooseSFStatUp   = &ntv.tauLooseSFStatUp;
  p_tauLooseSFStatDown = &ntv.tauLooseSFStatDown;
  p_tauLooseSFSystUp   = &ntv.tauLooseSFSystUp;
  p_tauLooseSFSystDown = &ntv.tauLooseSFSystDown; 

  p_systWeights = &ntv.systWeights;
}


void NTVarsRead::setAddresses(TTree* tree, bool addJetSmearSystW)
{
  tree->GetBranch("NTVars")->SetAddress(&ntv.RunNumber);
  tree->GetBranch("jetPt")->SetAddress(&p_jetPt);
  tree->GetBranch("jetEta")->SetAddress(&p_jetEta);
  tree->GetBranch("jetPhi")->SetAddress(&p_jetPhi);
  tree->GetBranch("jetM")->SetAddress(&p_jetM);
  tree->GetBranch("jetBTag")->SetAddress(&p_jetBTag);
  tree->GetBranch("jetFlav")->SetAddress(&p_jetFlav);
  tree->GetBranch("jetTagU")->SetAddress(&p_jetTagU);
  tree->GetBranch("jetTagB")->SetAddress(&p_jetTagB);
  tree->GetBranch("jetTagC")->SetAddress(&p_jetTagC);
  if ( addJetSmearSystW  ) {
    tree->GetBranch("jetSmearSystW")->SetAddress(&p_jetSmearSystW);
  }
  tree->GetBranch("jetFracSamplingMax")->SetAddress(&p_jetFracSamplingMax);
  tree->GetBranch("jetFracSamplingMaxIndex")->SetAddress(&p_jetFracSamplingMaxIndex);

  tree->GetBranch("tauPt")->SetAddress(&p_tauPt);
  tree->GetBranch("tauEta")->SetAddress(&p_tauEta);
  tree->GetBranch("tauPhi")->SetAddress(&p_tauPhi);
  tree->GetBranch("tauLooseSF")->SetAddress(&p_tauLooseSF);
  tree->GetBranch("tauLooseSFStatUp")->SetAddress(&p_tauLooseSFStatUp);
  tree->GetBranch("tauLooseSFStatDown")->SetAddress(&p_tauLooseSFStatDown);
  tree->GetBranch("tauLooseSFSystUp")->SetAddress(&p_tauLooseSFSystUp);
  tree->GetBranch("tauLooseSFSystDown")->SetAddress(&p_tauLooseSFSystDown);

  tree->GetBranch("systWeights")->SetAddress(&p_systWeights);

}


std::string NTReclusteringVars::toString()
{ 
  return std::string("NWcandidates/i:test/i:nJetsRecl/I"); 
}

void NTReclusteringVars::Reset()
{ 
  NWcandidates = 0;
  test = 0;
  nJetsRecl = 0;
  
  // Clear vectors
  RTjets10SubJetIndeces.clear();
  RTjetPt.clear();
  RTjetEta.clear();
  RTjetPhi.clear();
  RTjetM.clear();
  ReclJetMass.clear();
  ReclJetPt.clear();
  ReclJetEta.clear();
  ReclJetPhi.clear();
  D2.clear();
  isWmedium.clear();
  isWtight.clear();
  isZmedium.clear();
  isZtight.clear();

}

NTReclusteringVarsRead::NTReclusteringVarsRead(): RTntv()
{
  
  //RTjet variables
  p_RTjets10SubJetIndeces =  &RTntv.RTjets10SubJetIndeces;
  p_RTjetPt   = &RTntv.RTjetPt;
  p_RTjetEta  = &RTntv.RTjetEta;
  p_RTjetPhi  = &RTntv.RTjetPhi;
  p_RTjetM    = &RTntv.RTjetM;
  p_ReclJetMass = &RTntv.ReclJetMass;
  p_ReclJetPt   = &RTntv.ReclJetPt;
  p_ReclJetEta  = &RTntv.ReclJetEta;
  p_ReclJetPhi  = &RTntv.ReclJetPhi;
  p_D2          = &RTntv.D2;
  p_isWmedium   = &RTntv.isWmedium;
  p_isWtight    = &RTntv.isWtight;
  p_isZmedium   = &RTntv.isZmedium;
  p_isZtight    = &RTntv.isZtight;

}


void NTReclusteringVarsRead::setAddresses(TTree* tree)
{
  tree->GetBranch("NTReclusteringVars")->SetAddress(&RTntv.NWcandidates);
  
  //RTjet variables  
  tree->GetBranch("RTjets10SubJetIndeces")->SetAddress(&p_RTjets10SubJetIndeces);
  tree->GetBranch("RTjetPt")->SetAddress(&p_RTjetPt);
  tree->GetBranch("RTjetEta")->SetAddress(&p_RTjetEta);
  tree->GetBranch("RTjetPhi")->SetAddress(&p_RTjetPhi);
  tree->GetBranch("RTjetM")->SetAddress(&p_RTjetM);
  tree->GetBranch("ReclJetMass")->SetAddress(&p_ReclJetMass);
  tree->GetBranch("ReclJetPt")->SetAddress(&p_ReclJetPt);
  tree->GetBranch("ReclJetEta")->SetAddress(&p_ReclJetEta);
  tree->GetBranch("ReclJetPhi")->SetAddress(&p_ReclJetPhi);
  tree->GetBranch("D2")->SetAddress(&p_D2);
  tree->GetBranch("isWmedium")->SetAddress(&p_isWmedium);
  tree->GetBranch("isWtight")->SetAddress(&p_isWtight);
  tree->GetBranch("isZmedium")->SetAddress(&p_isZmedium);
  tree->GetBranch("isZtight")->SetAddress(&p_isZtight);

}

std::string NTCRWTVars::toString()
{ 
  return std::string("lep1Pt/F:lep1Eta/F:lep1Phi/F:lep1sign/I:mt/F:Wpt/F:dphilMET/F:Weta/F:lep1ptvarcone20/F:lep1ptvarcone30/F:lep1topoetcone20/F"); 
  //return std::string("lep1Pt/F:lep1Eta/F:lep1Phi/F:lep1sign/I:mt/F:Wpt/F:leptonWeight/F:leptonWeightUp/F:leptonWeightDown/F:triggerWeight/F:triggerWeightUp/F:triggerWeightDown/F:lep1Iso/F:lep1DRjet/F:lep1jetJVF/F"); 
}

void NTCRWTVars::Reset()
{
  lep1Pt = lep1Eta = lep1Phi = 0.f;
  lep1sign = 0;
  mt = 0.f;
  Wpt = 0.f;
  dphilMET = 0.f;
  Weta = 0.f;
  lep1ptvarcone20 = 0.f;
  lep1ptvarcone30 = 0.f;
  lep1topoetcone20 = 0.f;
  //leptonWeight     = 1.f;
  //leptonWeightUp   = 1.f;  
  //leptonWeightDown = 1.f;
  //triggerWeight     = 1.f;
  //triggerWeightUp   = 1.f;
  //triggerWeightDown = 1.f;
  //lep1Iso = 0.f;
  //lep1DRjet = 999.f;
  //lep1jetJVF = -999.f;
}
  



std::string NTCRZVars::toString()
{ 
  return std::string("lep1Pt/F:lep2Pt/F:lep1Eta/F:lep2Eta/F:lep1Phi/F:lep2Phi/F:lep1sign/I:lep2sign/I:mll/F:Zpt/F:leptonWeight/F:leptonWeightUp/F:leptonWeightDown/F:triggerWeight/F:triggerWeightUp/F:triggerWeightDown/F:fakemet/F:fakemetPhi/F:lep1Iso/F:lep2Iso/F:lep1DRjet/F:lep2DRjet/F:lep1jetJVF/F:lep2jetJVF/F"); 
}


void NTCRZVars::Reset()
{
  lep1Pt = lep1Eta = lep1Phi = 0.f;
  lep1sign = 0;
  lep2Pt = lep2Eta = lep2Phi = 0.f;
  lep2sign = 0;
  mll = 0.f;
  Zpt = 0.f;
  leptonWeight     = 1.f;
  leptonWeightUp   = 1.f;  
  leptonWeightDown = 1.f;
  triggerWeight     = 1.f;
  triggerWeightUp   = 1.f;
  triggerWeightDown = 1.f;
  fakemet = fakemetPhi = 0.f;
  lep1Iso = lep2Iso = 0.f;
  lep1DRjet = lep2DRjet = 999.f;
  lep1jetJVF = lep2jetJVF = -999.f;
}


std::string NTCR3LVars::toString()
{
  return std::string("lep1Pt/F:lep2Pt/F:lep3Pt/F:lep1Eta/F:lep2Eta/F:lep3Eta/F:lep1Phi/F:lep2Phi/F:lep3Phi/F:lep1sign/I:lep2sign/I:lep3sign/I:mll/F:Zpt/F:leptonWeight/F:leptonWeightUp/F:leptonWeightDown/F:triggerWeight/F:triggerWeightUp/F:triggerWeightDown/F:fakemet/F:fakemetPhi/F:lep1ptvarcone20/F:lep2ptvarcone20/F:lep3ptvarcone20/F:lep1ptvarcone30/F:lep2ptvarcone30/F:lep3ptvarcone30/F:lep1topoetcone20/F:lep2topoetcone20/F:lep3topoetcone20/F:lep1DRjet/F:lep2DRjet/F:lep3DRjet/F:lep1jetJVF/F:lep2jetJVF/F:lep3jetJVF/F:mt/F:Wpt/F:lepfromW/I:lepptfromW/F");
}

void NTCR3LVars::Reset()
{
  lep1Pt = lep1Eta = lep1Phi = 0.f;
  lep1sign = 0;
  lep2Pt = lep2Eta = lep2Phi = 0.f;
  lep2sign = 0;
  lep3Pt = lep3Eta = lep3Phi = 0.f;
  lep3sign = 0;
  lepfromW = 0;
  lepptfromW = 0 ; 
  mll = 0.f;
  Zpt = 0.f;
  mt = 0.f;
  Wpt = 0.f;
  leptonWeight     = 1.f;
  leptonWeightUp   = 1.f;
  leptonWeightDown = 1.f;
  triggerWeight     = 1.f;
  triggerWeightUp   = 1.f;
  triggerWeightDown = 1.f;
  fakemet = fakemetPhi = 0.f;
  lep1ptvarcone20 = lep2ptvarcone20 = lep3ptvarcone20 = 0.f;
  lep1ptvarcone30 = lep2ptvarcone30 = lep3ptvarcone30 = 0.f;
  lep1topoetcone20 = lep2topoetcone20 = lep3topoetcone20 = 0.f;
  lep1DRjet = lep2DRjet = lep3DRjet = 999.f;
  lep1jetJVF = lep2jetJVF = lep3jetJVF = -999.f;
}


std::string NTExtraVars::toString()
{ 
  return std::string("mettrack/F:mettrack_phi/F:mT2/F:mT2_noISR/F:Ap/F");

}

void NTExtraVars::Reset()
{ 
  mettrack = 0.f;
  mettrack_phi = 0.f;
  mT2 = -1.f;
  mT2_noISR = 0.f;
  mT2=0.f;
  mT2_noISR=0.f;
  Ap =0.f;
}

std::string NTRJigsawVars::toString()
{ 
  return std::string("RJVars_PP_Mass/F:RJVars_PP_InvGamma/F:RJVars_PP_dPhiBetaR/F:RJVars_PP_dPhiVis/F:RJVars_PP_CosTheta/F:RJVars_PP_dPhiDecayAngle/F:RJVars_PP_VisShape/F:RJVars_PP_MDeltaR/F:RJVars_P1_Mass/F:RJVars_P1_CosTheta/F:RJVars_P2_Mass/F:RJVars_P2_CosTheta/F:RJVars_I1_Depth/F:RJVars_I2_Depth/F:RJVars_V1_N/F:RJVars_V2_N/F:RJVars_MG/F:RJVars_DeltaBetaGG/F:RJVars_dphiVG/F:RJVars_P_0_CosTheta/F:RJVars_C_0_CosTheta/F:RJVars_P_0_dPhiGC/F:RJVars_P_0_MassRatioGC/F:RJVars_P_0_Jet1_pT/F:RJVars_P_0_Jet2_pT/F:RJVars_P_0_PInvHS/F:RJVars_P_1_CosTheta/F:RJVars_C_1_CosTheta/F:RJVars_P_1_dPhiGC/F:RJVars_P_1_MassRatioGC/F:RJVars_P_1_Jet1_pT/F:RJVars_P_1_Jet2_pT/F:RJVars_P_1_PInvHS/F:RJVars_QCD_dPhiR/F:RJVars_QCD_Rpt/F:RJVars_QCD_Rmsib/F:RJVars_QCD_Rpsib/F:RJVars_QCD_Delta1/F:RJVars_QCD_Delta2/F");
}

void NTRJigsawVars::Reset()
{ 

  RJVars_PP_Mass           =0.f; 
  RJVars_PP_InvGamma       =0.f; 
  RJVars_PP_dPhiBetaR      =0.f; 
  RJVars_PP_dPhiVis        =0.f; 
  RJVars_PP_CosTheta       =0.f; 
  RJVars_PP_dPhiDecayAngle =0.f; 
  RJVars_PP_VisShape       =0.f; 
  RJVars_PP_MDeltaR        =0.f; 
  RJVars_P1_Mass           =0.f; 
  RJVars_P1_CosTheta       =0.f; 
  RJVars_P2_Mass           =0.f; 
  RJVars_P2_CosTheta       =0.f; 
  RJVars_I1_Depth          =0.f; 
  RJVars_I2_Depth          =0.f; 
  RJVars_V1_N              =0.f; 
  RJVars_V2_N              =0.f;     
  RJVars_MG                =0.f;       
  RJVars_DeltaBetaGG       =0.f;       
  RJVars_dphiVG            =0.f;       
  RJVars_P_0_CosTheta      =0.f;       
  RJVars_C_0_CosTheta      =0.f;       
  RJVars_P_0_dPhiGC        =0.f;     
  RJVars_P_0_MassRatioGC   =0.f;   
  RJVars_P_0_Jet1_pT       =0.f; 
  RJVars_P_0_Jet2_pT       =0.f; 
  RJVars_P_0_PInvHS        =0.f;       
  RJVars_P_1_CosTheta      =0.f;       
  RJVars_C_1_CosTheta      =0.f;       
  RJVars_P_1_dPhiGC        =0.f;     
  RJVars_P_1_MassRatioGC   =0.f;      
  RJVars_P_1_Jet1_pT       =0.f; 
  RJVars_P_1_Jet2_pT       =0.f; 
  RJVars_P_1_PInvHS        =0.f;       
  RJVars_QCD_dPhiR         =0.f;  
  RJVars_QCD_Rpt           =0.f;  
  RJVars_QCD_Rmsib         =0.f;  
  RJVars_QCD_Rpsib         =0.f;  
  RJVars_QCD_Delta1        =0.f;  
  RJVars_QCD_Delta2        =0.f;  
}

std::string NTTheoryVars::toString()
{
  return std::string("mu1ScaleWeightUp/F:mu1ScaleWeightDown/F:mu2ScaleWeightUp/F:mu2ScaleWeightDown/F:matchScaleWeightUp/F:matchScaleWeightDown/F:HFWeight/F:nPartonsWeight/F:nTruthJet/I:nParton/I:nTau/I"); 
// scaleFormWeightUp/F:scaleFormWeightDown/F:
} 


void NTTheoryVars::Reset()
{ 
  mu1ScaleWeightUp = mu1ScaleWeightDown = mu2ScaleWeightUp = mu2ScaleWeightDown = matchScaleWeightUp = matchScaleWeightDown = HFWeight = nPartonsWeight = 1.f; // scaleFormWeightUp = scaleFormWeightDown = 
  nTruthJet = nParton = nTau =0;
}


void NTTheoryVars::Copy(NTTheoryVars& other)
{
  mu1ScaleWeightUp = other.mu1ScaleWeightUp; 
  mu1ScaleWeightDown = other.mu1ScaleWeightDown; 
  mu2ScaleWeightUp = other.mu2ScaleWeightUp; 
  mu2ScaleWeightDown = other.mu2ScaleWeightDown; 
  matchScaleWeightUp = other.matchScaleWeightUp; 
  matchScaleWeightDown = other.matchScaleWeightDown; 
  HFWeight = other.HFWeight;
  nPartonsWeight = other.nPartonsWeight;
  nTruthJet = other.nTruthJet;
  nParton = other.nParton;
  nTau = other.nTau;
}


std::string NTISRVars::toString()
{ 
  return std::string("nJetISR/I:isrjetIndex/I:isrjetPt/F:isrjetEta/F:isrjetPhi/F:jet1Alpha/F:jet2Alpha/F:jet3Alpha/F:jet4Alpha/F:jet5Alpha/F:jet1minPtDistinction/F:jet2minPtDistinction/F:jet3minPtDistinction/F:jet4minPtDistinction/F:jet5minPtDistinction/F:jet1minDeltaDistinction/F:jet2minDeltaDistinction/F:jet3minDeltaDistinction/F:jet4minDeltaDistinction/F:jet5minDeltaDistinction/F:jet1minEtaGap/F:jet2minEtaGap/F:jet3minEtaGap/F:jet4minEtaGap/F:jet5minEtaGap/F:jet1maxEtaOtherJets/F:jet2maxEtaOtherJets/F:jet3maxEtaOtherJets/F:jet4maxEtaOtherJets/F:jet5maxEtaOtherJets/F:jet1DPhiMET/F:jet2DPhiMET/F:jet3DPhiMET/F:jet4DPhiMET/F:jet5DPhiMET/F"); 
}

void NTISRVars::Reset()
{ 
  nISRJets = ISRjet_index = 0; 
  isrjetPt = 0.f;
  
  isrjetEta = isrjetPhi = jet1Alpha = jet2Alpha = jet3Alpha = jet4Alpha = jet5Alpha = jet1minPtDistinction = jet2minPtDistinction = jet3minPtDistinction = jet4minPtDistinction = jet5minPtDistinction = jet1minDeltaDistinction = jet2minDeltaDistinction = jet3minDeltaDistinction = jet4minDeltaDistinction = jet5minDeltaDistinction = jet1minEtaGap = jet2minEtaGap = jet3minEtaGap = jet4minEtaGap = jet5minEtaGap = jet1maxEtaOtherJets = jet2maxEtaOtherJets = jet3maxEtaOtherJets = jet4maxEtaOtherJets = jet5maxEtaOtherJets = jet1DPhiMET = jet2DPhiMET = jet3DPhiMET = jet4DPhiMET = jet5DPhiMET = 999.f;
}



void bookNTVars(TTree* tree, NTVars& ntv, bool addJetSmearSystW)
{
  tree->Branch("NTVars",&ntv,NTVars::toString().c_str());
  tree->Branch("jetPt",&(ntv.jetPt));
  tree->Branch("jetEta",&(ntv.jetEta));
  tree->Branch("jetPhi",&(ntv.jetPhi));
  tree->Branch("jetM",&(ntv.jetM));
  tree->Branch("jetBTag",&(ntv.jetBTag));
  tree->Branch("jetFlav",&(ntv.jetFlav));
  tree->Branch("jetTagU",&(ntv.jetTagU));
  tree->Branch("jetTagB",&(ntv.jetTagB));
  tree->Branch("jetTagC",&(ntv.jetTagC));                
  tree->Branch("jetFracSamplingMax",&(ntv.jetFracSamplingMax));
  tree->Branch("jetFracSamplingMaxIndex",&(ntv.jetFracSamplingMaxIndex));
  if ( addJetSmearSystW ) {
    tree->Branch("jetSmearSystW",&(ntv.jetSmearSystW)); 
  }

  tree->Branch("tauPt",&(ntv.tauPt));
  tree->Branch("tauEta",&(ntv.tauEta));
  tree->Branch("tauPhi",&(ntv.tauPhi));
  tree->Branch("tauLooseSF",&(ntv.tauLooseSF));
  tree->Branch("tauLooseSFStatUp",&(ntv.tauLooseSFStatUp));
  tree->Branch("tauLooseSFStatDown",&(ntv.tauLooseSFStatDown));
  tree->Branch("tauLooseSFSystUp",&(ntv.tauLooseSFSystUp));
  tree->Branch("tauLooseSFSystDown",&(ntv.tauLooseSFSystDown));

  tree->Branch("systWeights",&(ntv.systWeights));
  treePolicies(tree);
}

void bookNTReclusteringVars(TTree* tree, NTReclusteringVars& RTntv)
{
  tree->Branch("NTReclusteringVars",&RTntv,NTReclusteringVars::toString().c_str());
  
  //RT jets
  tree->Branch("RTjets10SubJetIndeces",&(RTntv.RTjets10SubJetIndeces));  
  tree->Branch("RTjetPt",&(RTntv.RTjetPt));
  tree->Branch("RTjetEta",&(RTntv.RTjetEta));
  tree->Branch("RTjetPhi",&(RTntv.RTjetPhi));
  tree->Branch("RTjetM",&(RTntv.RTjetM));
  tree->Branch("ReclJetMass",&(RTntv.ReclJetMass));
  tree->Branch("ReclJetPt",&(RTntv.ReclJetPt));
  tree->Branch("ReclJetEta",&(RTntv.ReclJetEta));
  tree->Branch("ReclJetPhi",&(RTntv.ReclJetPhi));
  tree->Branch("D2",&(RTntv.D2));
  tree->Branch("isWmedium",&(RTntv.isWmedium));
  tree->Branch("isWtight",&(RTntv.isWtight));
  tree->Branch("isZmedium",&(RTntv.isZmedium));
  tree->Branch("isZtight",&(RTntv.isZtight));
  treePolicies(tree);
}

std::string NTCRYVars::toString()
{
  return std::string("phPt/F:phEta/F:phPhi/F:origmet/F:origmetPhi/F:phTopoetcone20/F:phPtvarcone20/F:phPtcone20/F:phTopoetcone40/F:phPtvarcone40/F:phPtcone40/F:phLoose/I:phTight/I:phTruthType/I:phTruthOrigin/I:phisEMvalue/I:phSignal/I");
}

void NTCRYVars::Reset()
{
  phPt= 0.f;
  phEta= 0.f;
  phPhi= 0.f;
  origmet = 0.f;
  origmetPhi = 0.f;
  phTopoetcone20= 0.f;
  phPtvarcone20= 0.f;
  phPtcone20= 0.f;
  phTopoetcone40= 0.f;
  phPtvarcone40= 0.f;
  phPtcone40= 0.f;
  phLoose= 0;
  phTight= 0;
  phTruthType= 0;
  phTruthOrigin= 0;
  phisEMvalue= 0;
  phSignal= 0;
}


//NTCRYVarsRead::NTCRYVarsRead(): ntv()
//{
//  p_phPt = &ntv.phPt;
//  p_phEta = &ntv.phEta;
//  p_phPhi = &ntv.phPhi;
//  p_phSignal = &ntv.phSignal;
//  p_phTopoetcone20 = &ntv.phTopoetcone20;
//  p_phPtvarcone20 = &ntv.phPtvarcone20;
//  p_phPtcone20 = &ntv.phPtcone20;
//  p_phTopoetcone40 = &ntv.phTopoetcone40;
//  p_phPtvarcone40 = &ntv.phPtvarcone40;
//  p_phPtcone40 = &ntv.phPtcone40;
//  //p_phisEMTight  = &ntv.phisEMTight;
//  p_phLoose = &ntv.phLoose;
//  p_phTight = &ntv.phTight;
//  p_phTruthType = &ntv.phTruthType;
//  p_phTruthOrigin = &ntv.phTruthOrigin;
//  p_phisEMvalue = &ntv.phisEMvalue;
//}

//void bookNTCRYVars(TTree* tree, NTCRYVars& cryntv)
//{
//  tree->Branch("NTCRYVars",&cryntv,NTCRYVars::toString().c_str());
//  tree->Branch("phPt",&(cryntv.phPt));
//  tree->Branch("phEta",&(cryntv.phEta));
//  tree->Branch("phPhi",&(cryntv.phPhi));
//  tree->Branch("phSignal",&(cryntv.phSignal));
//  tree->Branch("phTopoetcone20",&(cryntv.phTopoetcone20));
//  tree->Branch("phPtvarcone20",&(cryntv.phPtvarcone20));
//  tree->Branch("phPtcone20",&(cryntv.phPtcone20));
//  tree->Branch("phTopoetcone40",&(cryntv.phTopoetcone40));
//  tree->Branch("phPtvarcone40",&(cryntv.phPtvarcone40));
//  tree->Branch("phPtcone40",&(cryntv.phPtcone40));
//  //tree->Branch("phisEMTight",&(cryntv.phisEMTight));
//  tree->Branch("phLoose",&(cryntv.phLoose));
//  tree->Branch("phTight",&(cryntv.phTight));
//  tree->Branch("phTruthType",&(cryntv.phTruthType));
//  tree->Branch("phTruthOrigin",&(cryntv.phTruthOrigin));
//  tree->Branch("phisEMvalue",&(cryntv.phisEMvalue));
//  treePolicies(tree);
//}
//
//void NTCRYVarsRead::setAddresses(TTree* tree)
//{
//  tree->GetBranch("NTCRYVars")->SetAddress(&ntv.origmet);
//  tree->GetBranch("phPt")->SetAddress(&p_phPt);
//  tree->GetBranch("phEta")->SetAddress(&p_phEta);
//  tree->GetBranch("phPhi")->SetAddress(&p_phPhi);
//  tree->GetBranch("phSignal")->SetAddress(&p_phSignal);
//  tree->GetBranch("phTopoetcone20")->SetAddress(&p_phTopoetcone20);
//  tree->GetBranch("phPtvarcone20")->SetAddress(&p_phPtvarcone20);
//  tree->GetBranch("phPtcone20")->SetAddress(&p_phPtcone20);
//  tree->GetBranch("phTopoetcone40")->SetAddress(&p_phTopoetcone40);
//  tree->GetBranch("phPtvarcone40")->SetAddress(&p_phPtvarcone40);
//  tree->GetBranch("phPtcone40")->SetAddress(&p_phPtcone40);
//  //tree->GetBranch("phisEMTight")->SetAddress(&p_phisEMTight);
//  tree->GetBranch("phLoose")->SetAddress(&p_phLoose);
//  tree->GetBranch("phTight")->SetAddress(&p_phTight);
//  tree->GetBranch("phTruthType")->SetAddress(&p_phTruthType);
//  tree->GetBranch("phTruthOrigin")->SetAddress(&p_phTruthOrigin);
//  tree->GetBranch("phisEMvalue")->SetAddress(&p_phisEMvalue);
//}

