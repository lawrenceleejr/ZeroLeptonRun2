

// cut values
#ifndef ZeroLeptonCutVal_h
#define ZeroLeptonCutVal_h


class ZeroLeptonCutVal
{
 public:

  ZeroLeptonCutVal();
  void ReadCutValues(std::string);

  double m_cutEtMiss;
  double m_cutEtMiss1Jet;
  double m_cutJetPt0;
  double m_cutJetPt1;
  double m_cutJetPt2;
  double m_cutJetPt3;
  double m_cutJetPt4;
  double m_cutJetPt5;
  double m_cutJetPtWRes;
  
  float m_cutMETsig[8][4];
  float m_cutMEToverMeff[8][4];
  float m_cutMeff[8][4];
  float m_CRWT_cutMETsig[8][4];
  float m_CRWT_cutMEToverMeff[8][4];
  float m_CRWT_cutMeff[8][4];
  bool  m_cutMap[8][4];
  
  double m_cutDeltaPhi;
  double m_cutDeltaPhi2;
  double m_cutDeltaPhiQCD;

  double m_cutPhotonPtCRY;
  double m_cutPhotonIsolCRY;
  double m_cutPhotonQualityCRY;

  double m_cutQCDSmearPt;
  double m_cutQCDSmearEta;
};

#endif
