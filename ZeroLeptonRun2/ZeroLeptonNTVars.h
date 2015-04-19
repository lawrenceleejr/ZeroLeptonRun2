#ifndef ZeroLeptonNTVars_h
#define ZeroLeptonNTVars_h

#include "TTree.h"

#include<vector>
#include<string>

#define NMAXPDF 60

class NTVars
{
 public:
  NTVars() { Reset(); }

  static std::string toString();
  
  void Reset();

  unsigned int RunNumber, EventNumber, veto;
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
  float tauMt;
  float SherpaBugMET;
  float metNOCHCORRCELL;
  float metLHTOPO;
  float metLHTOPONOCHCORRCELL;

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
  std::vector< float > jetBCH_CORR_CELL;
  std::vector< float > jetFracSamplingMax;
  std::vector< float > jetFracSamplingMaxIndex;

  
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
  std::vector< float >* p_jetBCH_CORR_CELL;
  std::vector< float >* p_jetFracSamplingMax;
  std::vector< float >* p_jetFracSamplingMaxIndex;
 
};


class NTReclusteringVars
{
 public:
  NTReclusteringVars() { Reset(); }

  static std::string toString();
  
  void Reset();

  unsigned int NWcandidates;
  int test;
  
  // WARNING: if you add another vector you need to update NTVarsRead below

  // STL vectors of RT jets  
  std::vector< std::vector< int > > RTjets10SubJetIndeces;
  std::vector< float > RTjetM;  
  std::vector< float > RTjetPt;
  std::vector< float > RTjetEta;
  std::vector< float > RTjetPhi;

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
  //float leptonWeight, leptonWeightUp, leptonWeightDown;
  //float triggerWeight, triggerWeightUp, triggerWeightDown;
  //float lep1Iso;
  //float lep1DRjet;
  //float lep1jetJVF;
};

class NTCRYVars
{
 public:
  NTCRYVars() { Reset(); }
  
  static std::string toString();


  void Reset();
  
  int nPh;
  unsigned int phQuality;
  float phPt, phPhi, phEta, phIso;
  int nLep;
  float lep1Pt, lep2Pt;
  float lep1Eta, lep2Eta;
  float lep1Phi, lep2Phi;
  int lep1sign, lep2sign;
  float photonWeight, photonWeightUp, photonWeightDown;
  float leptonWeight, leptonWeightUp, leptonWeightDown;
  float triggerWeight, triggerWeightUp, triggerWeightDown;
  float mass;
  float fakemet,fakemetPhi;
  float dRPhLep;
  float lep1Iso, lep2Iso;
  float lep1DRjet, lep2DRjet;
  float lep1jetJVF, lep2jetJVF;
};

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
  float gaminvRp1 ;
  float shatR ;
  float mdeltaR ;
  float cosptR ;
  float gamma_R;
  float dphi_BETA_R ; 
  float dphi_leg1_leg2 ; 
  float costhetaR ;
  float dphi_BETA_Rp1_BETA_R;
  float gamma_Rp1; 
  float costhetaRp1;
  float Ap;
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

void bookNTVars(TTree* tree, NTVars& ntv, bool addJetSmearSystW);

void bookNTReclusteringVars(TTree* tree, NTReclusteringVars& RTntv);

inline void bookNTExtraVars(TTree* tree, NTExtraVars& extrantv)
{
  tree->Branch("NTExtraVars",&extrantv,NTExtraVars::toString().c_str());
}

inline void bookNTTheoryVars(TTree* tree, NTTheoryVars& theoryntv)
{
  tree->Branch("NTTheoryVars",&theoryntv,NTTheoryVars::toString().c_str());
}

inline void bookNTISRVars(TTree* tree, NTISRVars& isrntv)
{
  tree->Branch("NTISRVars",&isrntv,NTISRVars::toString().c_str());
}

inline void bookNTCRWTVars(TTree* tree, NTCRWTVars& crwtntv)
{
  tree->Branch("NTCRWTVars",&crwtntv,NTCRWTVars::toString().c_str());
}

inline void bookNTCRZVars(TTree* tree, NTCRZVars& crzntv)
{
  tree->Branch("NTCRZVars",&crzntv,NTCRZVars::toString().c_str());
}


inline void bookNTCRYVars(TTree* tree, NTCRYVars& cryntv)
{
  tree->Branch("NTCRYVars",&cryntv,NTCRYVars::toString().c_str());
}



#endif
