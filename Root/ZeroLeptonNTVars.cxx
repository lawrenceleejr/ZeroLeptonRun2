#include "ZeroLeptonRun2/ZeroLeptonNTVars.h"
#include "TTree.h"

std::string NTVars::toString()
{ 
  return
 
std::string("RunNumber/i:EventNumber/i:veto/i:eventWeight/F:pileupWeight/F:pileupWeightUp/F:pileupWeightDown/F:genWeight/F:ttbarWeightHT/F:ttbarWeightPt2/F:ttbarAvgPt/F:WZweight/F:nJet/i:met/F:metPhi/F:dPhi/F:dPhiR/F:meffInc/F:hardproc/I:nBJet/I:nCJet/I:bTagWeight/F:bTagWeightBUp/F:bTagWeightBDown/F:bTagWeightCUp/F:bTagWeightCDown/F:bTagWeightLUp/F:bTagWeightLDown/F:cTagWeight/F:cTagWeightBUp/F:cTagWeightBDown/F:cTagWeightCUp/F:cTagWeightCDown/F:cTagWeightLUp/F:cTagWeightLDown/F:normWeight/F:normWeightUp/F:normWeightDown/F:cleaning/i:timing/F:jet1Emf/F:jet2Emf/F:jet1Chf/F:jet2Chf/F:pdfId1/I:pdfId2/I:tauN/i:tauJetBDTLoose/i:tauMt/F:SherpaBugMET/F:metNOCHCORRCELL/F:metLHTOPO/F:metLHTOPONOCHCORRCELL/F"); 
}

void NTVars::Reset()
{
  RunNumber = 0;
  EventNumber= 0;
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
  tauMt = 0.f;
  SherpaBugMET = 0.f;
  metNOCHCORRCELL = 0.f;
  metLHTOPO = 0.f;
  metLHTOPONOCHCORRCELL = 0.f;

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
  jetBCH_CORR_CELL.clear();
  jetFracSamplingMax.clear();
  jetFracSamplingMaxIndex.clear();
  
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
  p_jetBCH_CORR_CELL = &ntv.jetBCH_CORR_CELL;
  p_jetFracSamplingMax = &ntv.jetFracSamplingMax;
  p_jetFracSamplingMaxIndex = &ntv.jetFracSamplingMaxIndex;
  
 
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
  tree->GetBranch("jetBCH_CORR_CELL")->SetAddress(&p_jetBCH_CORR_CELL);
  tree->GetBranch("jetFracSamplingMax")->SetAddress(&p_jetFracSamplingMax);
  tree->GetBranch("jetFracSamplingMaxIndex")->SetAddress(&p_jetFracSamplingMaxIndex);

}


std::string NTReclusteringVars::toString()
{ 
  return std::string("NWcandidates/i:test/i"); 
}

void NTReclusteringVars::Reset()
{ 
  NWcandidates = 0;
  test = 0;
  
  // Clear vectors
  RTjets10SubJetIndeces.clear();
  RTjetPt.clear();
  RTjetEta.clear();
  RTjetPhi.clear();
  RTjetM.clear();
  
}

NTReclusteringVarsRead::NTReclusteringVarsRead(): RTntv()
{
  
  //RTjet variables
  p_RTjets10SubJetIndeces =  &RTntv.RTjets10SubJetIndeces;
  p_RTjetPt   = &RTntv.RTjetPt;
  p_RTjetEta  = &RTntv.RTjetEta;
  p_RTjetPhi  = &RTntv.RTjetPhi;
  p_RTjetM    = &RTntv.RTjetM;

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
}

std::string NTCRWTVars::toString()
{ 
  return std::string("lep1Pt/F:lep1Eta/F:lep1Phi/F:lep1sign/I:mt/F:Wpt/F"); 
  //return std::string("lep1Pt/F:lep1Eta/F:lep1Phi/F:lep1sign/I:mt/F:Wpt/F:leptonWeight/F:leptonWeightUp/F:leptonWeightDown/F:triggerWeight/F:triggerWeightUp/F:triggerWeightDown/F:lep1Iso/F:lep1DRjet/F:lep1jetJVF/F"); 
}

void NTCRWTVars::Reset()
{
  lep1Pt = lep1Eta = lep1Phi = 0.f;
  lep1sign = 0;
  mt = 0.f;
  Wpt = 0.f;
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



std::string NTExtraVars::toString()
{ 
  return std::string("mettrack/F:mettrack_phi/F:mT2/F:mT2_noISR/F:gaminvRp1/F:shatR/F:mdeltaR/F:cosptR/F:gamma_R/F:dphi_BETA_R/F:dphi_leg1_leg2/F:costhetaR/F:dphi_BETA_Rp1_BETA_R/F:gamma_Rp1/F:costhetaRp1/F:Ap/F");

}

void NTExtraVars::Reset()
{ 
  mettrack = 0.f;
  mettrack_phi = 0.f;
  mT2 = -1.f;
  mT2_noISR = 0.f;
  mT2=0.f;
  mT2_noISR=0.f;
  gaminvRp1 =0.f;
  shatR =0.f;
  mdeltaR =0.f;
  cosptR =0.f;
  gamma_R=0.f;
  dphi_BETA_R =0.f; 
  dphi_leg1_leg2 =0.f; 
  costhetaR =0.f;
  dphi_BETA_Rp1_BETA_R=0.f;
  gamma_Rp1=0.f; 
  costhetaRp1=0.f; 
  Ap =0.f;
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
  tree->Branch("jetBCH_CORR_CELL",&(ntv.jetBCH_CORR_CELL));
  tree->Branch("jetFracSamplingMax",&(ntv.jetFracSamplingMax));
  tree->Branch("jetFracSamplingMaxIndex",&(ntv.jetFracSamplingMaxIndex));
  if ( addJetSmearSystW ) {
    tree->Branch("jetSmearSystW",&(ntv.jetSmearSystW)); 
  }
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

}


void NTCRYVars::Reset()
{
  phPt.clear();
  phEta.clear();
  phPhi.clear();
  phSignal.clear();
  origmet = 0.f;
  origmetPhi = 0.f;
}


NTCRYVarsRead::NTCRYVarsRead(): ntv()
{
  p_phPt = &ntv.phPt;
  p_phEta = &ntv.phEta;
  p_phPhi = &ntv.phPhi;
  p_phSignal = &ntv.phSignal;
}

void bookNTCRYVars(TTree* tree, NTCRYVars& cryntv)
{
  tree->Branch("NTCRYVars",&cryntv,NTCRYVars::toString().c_str());
  tree->Branch("phPt",&(cryntv.phPt));
  tree->Branch("phEta",&(cryntv.phEta));
  tree->Branch("phPhi",&(cryntv.phPhi));
  tree->Branch("phSignal",&(cryntv.phSignal));
}

void NTCRYVarsRead::setAddresses(TTree* tree)
{
  tree->GetBranch("NTCRYVars")->SetAddress(&ntv.origmet);
  tree->GetBranch("phPt")->SetAddress(&p_phPt);
  tree->GetBranch("phEta")->SetAddress(&p_phEta);
  tree->GetBranch("phPhi")->SetAddress(&p_phPhi);
  tree->GetBranch("phSignal")->SetAddress(&p_phSignal);
}

