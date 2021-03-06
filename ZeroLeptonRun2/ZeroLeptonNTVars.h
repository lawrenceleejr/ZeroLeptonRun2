#ifndef ZeroLeptonNTVars_h
#define ZeroLeptonNTVars_h

#include "TTree.h"

#include<vector>
#include<string>
#include<bitset>
#define NMAXPDF 60

class NTVars
{
 public:
  NTVars() { Reset(); }

  static std::string toString();
  void Reset();

  unsigned int RunNumber, EventNumber, LumiBlockNumber, veto;
  float weight, pileupWeight, pileupWeightUp, pileupWeightDown, genWeight;
  float ttbarWeightHT, ttbarWeightPt2, ttbarAvgPt, WZweight;
  int Njet;
  float MET, METPhi, deltaPhi, deltaPhiRemaining, MeffIncl;
  //unsigned int region;
  int hardproc;
  int nBJet;
  int nCJet;
  float bTagWeight, bTagWeightBUp, bTagWeightBDown, bTagWeightCUp, bTagWeightCDown, bTagWeightLUp, bTagWeightLDown;
  float cTagWeight, cTagWeightBUp, cTagWeightBDown, cTagWeightCUp, cTagWeightCDown, cTagWeightLUp, cTagWeightLDown;
  float normWeight,normWeightUp,normWeightDown;
  unsigned int cleaning;
  float timing;
  float emfjet0,emfjet1;
  float chfjet0,chfjet1;
  int pdfId1, pdfId2;
  int tauN, tauJetBDTLoose;
  int tauLooseN;
  float tauMt;
  float SherpaBugMET;
  float dPhiBadTile;

  // WARNING: if you add another vector you need to update NTVarsRead below

  // STL vectors, store 4-momenta, b-tagged weight and truth flavor of
  // every jet with pT above a threshold (40 GeV)
  std::vector< float > jetPt;
  std::vector< float > jetEta;
  std::vector< float > jetPhi;
  std::vector< float > jetM;
  std::vector< float > jetBTag;
  std::vector< int > jetFlav;
  std::vector< float > jetTagU;
  std::vector< float > jetTagB;
  std::vector< float > jetTagC;
  std::vector< float > jetSmearSystW;
  std::vector< float > jetFracSamplingMax;
  std::vector< float > jetFracSamplingMaxIndex;

  std::vector< float > tauPt;
  std::vector< float > tauEta;
  std::vector< float > tauPhi;
  std::vector< float > tauLooseSF;
  std::vector< float > tauLooseSFStatUp;
  std::vector< float > tauLooseSFStatDown;
  std::vector< float > tauLooseSFSystUp;
  std::vector< float > tauLooseSFSystDown;

  std::vector< float > systWeights;
  std::vector< float > triggerBits;
};

// this is a helper class to read an ntuple with an NTVars, as it need the
// address of pointer to the vector<> like:
/*
    NTVarsRead inVars;
    TFile* f = TFile::Open("myfile.root");
    TTree* tree = dynamic_cast<TTree*>(f->Get("SRAllNT"));
    inVars.setAddresses(tree);
    tree.GetEntry(1);
*/

class NTVarsRead
{
public:
  NTVarsRead();
  void setAddresses(TTree* tree, bool addJetSmearSystW = false);

  NTVars ntv;
private:
  std::vector< float >* p_jetPt;
  std::vector< float >* p_jetEta;
  std::vector< float >* p_jetPhi;
  std::vector< float >* p_jetM;
  std::vector< float >* p_jetBTag;
  std::vector< int >* p_jetFlav;
  std::vector< float >* p_jetTagU;
  std::vector< float >* p_jetTagB;
  std::vector< float >* p_jetTagC;
  std::vector< float >* p_jetSmearSystW;
  std::vector< float >* p_jetFracSamplingMax;
  std::vector< float >* p_jetFracSamplingMaxIndex;
  std::vector< float >* p_tauPt;
  std::vector< float >* p_tauEta;
  std::vector< float >* p_tauPhi;
  std::vector< float >* p_tauLooseSF;
  std::vector< float >* p_tauLooseSFStatUp;
  std::vector< float >* p_tauLooseSFStatDown;
  std::vector< float >* p_tauLooseSFSystUp;
  std::vector< float >* p_tauLooseSFSystDown;

  std::vector< float >* p_systWeights;
  std::vector< float >* p_triggerBits;
};


class NTReclusteringVars
{
 public:
  NTReclusteringVars() { Reset(); }

  static std::string toString();
  void Reset();

  unsigned int NWcandidates;
  int test;
  int nJetsRecl;

  // WARNING: if you add another vector you need to update NTVarsRead below

  // STL vectors of RT jets
  std::vector< std::vector< int > > RTjets10SubJetIndeces;
  std::vector< float > RTjetM;
  std::vector< float > RTjetPt;
  std::vector< float > RTjetEta;
  std::vector< float > RTjetPhi;
  std::vector<float> ReclJetMass;
  std::vector<float> ReclJetPt;
  std::vector<float> ReclJetEta;
  std::vector<float> ReclJetPhi;
  std::vector<float> D2;
  std::vector<bool> isWmedium;
  std::vector<bool> isWtight;
  std::vector<bool> isZmedium;
  std::vector<bool> isZtight;
};


class NTReclusteringVarsRead
{
public:
  NTReclusteringVarsRead();
  void setAddresses(TTree* tree);

  NTReclusteringVars RTntv;
private:
  std::vector< std::vector< int > >* p_RTjets10SubJetIndeces;
  std::vector< float >* p_RTjetM;
  std::vector< float >* p_RTjetPt;
  std::vector< float >* p_RTjetEta;
  std::vector< float >* p_RTjetPhi;
  std::vector< float >* p_ReclJetMass;
  std::vector< float >* p_ReclJetPt;
  std::vector< float >* p_ReclJetEta;
  std::vector< float >* p_ReclJetPhi;
  std::vector< float >* p_D2;
  std::vector< bool >* p_isWmedium;
  std::vector< bool >* p_isWtight;
  std::vector< bool >* p_isZmedium;
  std::vector< bool >* p_isZtight;

};

class NTCRZVars
{
 public:
  NTCRZVars() { Reset(); }

  static std::string toString();

  void Reset();

  float lep1Pt, lep2Pt;
  float lep1Eta, lep2Eta;
  float lep1Phi, lep2Phi;
  int lep1sign, lep2sign;
  float mll, Zpt;
  float leptonWeight, leptonWeightUp, leptonWeightDown;
  float triggerWeight, triggerWeightUp, triggerWeightDown;
  float fakemet,fakemetPhi;
  float lep1Iso, lep2Iso;
  float lep1DRjet, lep2DRjet;
  float lep1jetJVF, lep2jetJVF;
};

class NTCRWTVars
{
 public:
  NTCRWTVars() { Reset(); }

  static std::string toString();
  void Reset();
  float lep1Pt, lep1Eta, lep1Phi;
  int lep1sign;
  float mt, Wpt;
  float dphilMET;
  float Weta;
  float lep1ptvarcone20, lep1ptvarcone30, lep1topoetcone20;
  //float leptonWeight, leptonWeightUp, leptonWeightDown;
  //float triggerWeight, triggerWeightUp, triggerWeightDown;
  //float lep1Iso;
  //float lep1DRjet;
  //float lep1jetJVF;
};


class NTCR3LVars
{
 public:
  NTCR3LVars() { Reset(); }

  static std::string toString();

  void Reset();

  float lep1Pt, lep2Pt, lep3Pt;
  float lep1Eta, lep2Eta, lep3Eta;
  float lep1Phi, lep2Phi, lep3Phi;
  int lep1sign, lep2sign, lep3sign;
  float mll, Zpt;
  float leptonWeight, leptonWeightUp, leptonWeightDown;
  float triggerWeight, triggerWeightUp, triggerWeightDown;
  float fakemet,fakemetPhi;
  float lep1ptvarcone20, lep2ptvarcone20, lep3ptvarcone20;
  float lep1ptvarcone30, lep2ptvarcone30, lep3ptvarcone30;
  float lep1topoetcone20, lep2topoetcone20, lep3topoetcone20;
  float lep1DRjet, lep2DRjet, lep3DRjet;
  float lep1jetJVF, lep2jetJVF, lep3jetJVF;
  float mt, Wpt;
  int lepfromW;
  float lepptfromW;
};


class NTCRYVars
{
 public:
  NTCRYVars() { Reset(); }

  static std::string toString();// {return std::string("origmet/F:origmetPhi/F");}
  void Reset();

  // WARNING: if you add another vector you need to update NTCRYVarsRead below
  float phPt;
  float phEta;
  float phPhi;
  float origmet;
  float origmetPhi;
  float phTopoetcone20;
  float phPtvarcone20;
  float phPtcone20;
  float phTopoetcone40;
  float phPtvarcone40;
  float phPtcone40;
  int   phLoose;
  int   phTight;
  int   phTruthType;
  int   phTruthOrigin;
  int   phisEMvalue;
  int   phSignal;

};

//class NTCRYVarsRead
//{
//public:
//  NTCRYVarsRead();
//  void setAddresses(TTree* tree);
//
//  NTCRYVars ntv;
//private:
//  std::vector< float >* p_phPt;
//  std::vector< float >* p_phEta;
//  std::vector< float >* p_phPhi;
//  std::vector< bool >*  p_phSignal;
//  std::vector< float >* p_phTopoetcone20;
//  std::vector< float >* p_phPtvarcone20;
//  std::vector< float >* p_phPtcone20;
//  std::vector< float >* p_phTopoetcone40;
//  std::vector< float >* p_phPtvarcone40;
//  std::vector< float >* p_phPtcone40;
//  //std::vector< int >*   p_phisEMTight;
//  std::vector< bool >*  p_phLoose;
//  std::vector< bool >*  p_phTight;
//  std::vector< int  >*  p_phTruthType;
//  std::vector< int  >*  p_phTruthOrigin;
//  std::vector< int  >*  p_phisEMvalue;
//};

class NTPdfVars {
public:
  NTPdfVars() {
    Reset(0);
  }
  void Reset(int n)
  {
    NpdfWeight = n; // no need to initialize pdfWeight[] as it is dynamic aray
  }
  int NpdfWeight;
  float pdfWeight[NMAXPDF];
};

class NTExtraVars {
public:
  NTExtraVars() {
    Reset();
  }

  static std::string toString();

  void Reset();

  float mettrack;
  float mettrack_phi;
  float mT2;
  float mT2_noISR;
  float Ap;
};

class NTRJigsawVars {
public:

  NTRJigsawVars() {
    Reset();
  }


  static std::string toString();

  void Print();
  void Reset();

  float HT1CM;
  float PIoHT1CM;
  float cosS;
  float NVS;
  float RPT_HT1CM;
  float MS;
  float ddphiP;
  float sdphiP;
  float pPP_Ia;
  float pPP_Ib;
  float pT_jet1a;
  float pT_jet1b;
  float pT_jet2a;
  float pT_jet2b;
  float pTPP_jet1;
  float pTPP_jet2;
  float pTPP_jet1a;
  float pTPP_jet1b;
  float pTPP_jet2a;
  float pTPP_jet2b;
  float pTPP_jet3a;
  float pTPP_jet3b;
  float pPP_jet1a;
  float pPP_jet1b;
  float pPP_jet2a;
  float pPP_jet2b;
  float pPP_jet3a;
  float pPP_jet3b;
  float R_H2PP_H3PP;
  float R_pTj2_HT3PP;
  float R_HT5PP_H5PP;
  float R_H2PP_H5PP;
  float minR_pTj2i_HT3PPi;
  float maxR_H1PPi_H2PPi;
  float R_HT9PP_H9PP;
  float R_H2PP_H9PP;
  float RPZ_HT3PP;
  float RPZ_HT5PP;
  float RPZ_HT9PP;
  float RPT_HT3PP;
  float RPT_HT5PP;
  float RPT_HT9PP;
  float PP_InvGamma;
  float PP_dPhiBetaR;
  float PP_dPhiVis;
  float PP_CosTheta;
  float PP_dPhiDecayAngle;
  float PP_VisShape;
  float PP_MDeltaR;
  float P1_Mass;
  float P1_CosTheta;
  float P2_Mass;
  float P2_CosTheta;
  float I1_Depth;
  float I2_Depth;
  float dphiPV1a;
  float cosV1a;
  float dphiCV2a;
  float cosV2a;
  float dphiPV1b;
  float cosV1b;
  float dphiCV2b;
  float cosV2b;
  float NJa;
  float NJb;
  float dphiVG;
  float QCD_dPhiR;
  float QCD_Rsib;
  float QCD_Delta1;
  float H2PP;
  float H3PP;
  float H4PP;
  float H6PP;
  float H10PP;
  float HT10PP;
  float H2Pa;
  float H2Pb;
  float H3Pa;
  float H3Pb;
  float H4Pa;
  float H4Pb;
  float H5Pa;
  float H5Pb;
  float H2Ca;
  float H2Cb;
  float H3Ca;
  float H3Cb;
  float HT4PP;
  float HT6PP;
  float sangle;
  float dangle;
};

class NTTheoryVars {
public:
  NTTheoryVars() {
    Reset();
  }

  static std::string toString();

  void Reset();

  void Copy(NTTheoryVars& other);

  float mu1ScaleWeightUp, mu1ScaleWeightDown, mu2ScaleWeightUp, mu2ScaleWeightDown, matchScaleWeightUp, matchScaleWeightDown, HFWeight, nPartonsWeight; // , scaleFormWeightUp, scaleFormWeightDown
  int nTruthJet, nParton,nTau;
};

class NTISRVars {
public:
  NTISRVars() {
    Reset();
  }
  static std::string toString();

  void Reset();

  int nISRJets, ISRjet_index;
  float isrjetPt, isrjetEta, isrjetPhi;
  float jet1Alpha,jet2Alpha,jet3Alpha,jet4Alpha,jet5Alpha,jet1minPtDistinction,jet2minPtDistinction,jet3minPtDistinction,jet4minPtDistinction,jet5minPtDistinction,jet1minDeltaDistinction,jet2minDeltaDistinction,jet3minDeltaDistinction,jet4minDeltaDistinction,jet5minDeltaDistinction,jet1minEtaGap,jet2minEtaGap,jet3minEtaGap,jet4minEtaGap,jet5minEtaGap,jet1maxEtaOtherJets,jet2maxEtaOtherJets,jet3maxEtaOtherJets,jet4maxEtaOtherJets,jet5maxEtaOtherJets,jet1DPhiMET,jet2DPhiMET,jet3DPhiMET,jet4DPhiMET,jet5DPhiMET;
};

inline void treePolicies(TTree* tree)
{
  tree->SetAutoSave(-10000000);
  tree->SetAutoFlush(-10000000);
  tree->SetCacheSize(0);
  //tree->SetCacheSize(5000000);
}

void bookNTVars(TTree* tree, NTVars& ntv, bool addJetSmearSystW);

void bookNTReclusteringVars(TTree* tree, NTReclusteringVars& RTntv);

inline void bookNTExtraVars(TTree* tree, NTExtraVars& extrantv)
{
  tree->Branch("NTExtraVars",&extrantv,NTExtraVars::toString().c_str());
  treePolicies(tree);
}


inline void bookNTRJigsawVars(TTree* tree, NTRJigsawVars& rjigsawntv)
{
  tree->Branch("NTRJigsawVars",&rjigsawntv,NTRJigsawVars::toString().c_str());
  treePolicies(tree);
}

inline void bookNTTheoryVars(TTree* tree, NTTheoryVars& theoryntv)
{
  tree->Branch("NTTheoryVars",&theoryntv,NTTheoryVars::toString().c_str());
  treePolicies(tree);
}

inline void bookNTISRVars(TTree* tree, NTISRVars& isrntv)
{
  tree->Branch("NTISRVars",&isrntv,NTISRVars::toString().c_str());
  treePolicies(tree);
}

inline void bookNTCRWTVars(TTree* tree, NTCRWTVars& crwtntv)
{
  tree->Branch("NTCRWTVars",&crwtntv,NTCRWTVars::toString().c_str());
  treePolicies(tree);
}

inline void bookNTCRZVars(TTree* tree, NTCRZVars& crzntv)
{
  tree->Branch("NTCRZVars",&crzntv,NTCRZVars::toString().c_str());
  treePolicies(tree);
}

inline void bookNTCR3LVars(TTree* tree, NTCR3LVars& cr3lntv)
{
  tree->Branch("NTCR3LVars",&cr3lntv,NTCR3LVars::toString().c_str());
  treePolicies(tree);
}

inline void bookNTCRYVars(TTree* tree, NTCRYVars& cryntv)
{
  tree->Branch("NTCRYVars",&cryntv,NTCRYVars::toString().c_str());
  treePolicies(tree);
}


#endif
