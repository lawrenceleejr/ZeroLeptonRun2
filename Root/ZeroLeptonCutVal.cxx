
#include "ZeroLeptonRun2/ZeroLeptonCxxUtils.h"
#include "ZeroLeptonRun2/ZeroLeptonParamFile.h"
#include "ZeroLeptonRun2/ZeroLeptonCutVal.h"

ZeroLeptonCutVal::ZeroLeptonCutVal()
{
  m_cutEtMiss = -1.;
  m_cutEtMissTruthTest = -1.;
  m_cutEtMiss1Jet = -1.;
  m_cutJetPt0 = -1.;
  m_cutJetPt1 = -1.;
  m_cutJetPt2 = -1.;
  m_cutJetPt3 = -1.;
  m_cutJetPt4 = -1.;
  m_cutJetPt5 = -1.;

  m_cutRJigsawJetPt = -1.;
  
  for ( size_t i = 0; i < 8; i++ )
  {
    for ( size_t j = 0; j < 4; j++ )
    {
      m_cutMETsig[i][j] = -1.;
      m_cutMEToverMeff[i][j] = -1.;
      m_cutMeff[i][j] = -1.;
      m_CRWT_cutMETsig[i][j] = -1.;
      m_CRWT_cutMEToverMeff[i][j] = -1.;
      m_CRWT_cutMeff[i][j] = -1.;
      m_cutMap[i][j] = false;
    }
  }
  
  m_cutDeltaPhi = -1.;
  m_cutDeltaPhi2 = -1.;
  m_cutDeltaPhiQCD = -1.;

  m_cutPhotonPtCRY = -1.;
  m_cutPhotonIsolCRY = -1.;
  m_cutPhotonQualityCRY = -1.;

  m_cutQCDSmearPt = -1.;
  m_cutQCDSmearEta = -1.;
}

void ZeroLeptonCutVal::ReadCutValues(std::string paramfilename)
{
  
  ZeroLeptonParamFile params(paramfilename);
  m_cutEtMiss=params.find<double>("cutEtMiss",160000);
  m_cutEtMissTruthTest=params.find<double>("cutEtMissTruthTest",160000);
  m_cutEtMiss1Jet=params.find<double>("cutEtMiss1Jet",300000);
  m_cutJetPt0=params.find<double>("cutJetPt0",130000);
  m_cutJetPt1=params.find<double>("cutJetPt1",50000);
  m_cutJetPt2=params.find<double>("cutJetPt2",50000);
  m_cutJetPt3=params.find<double>("cutJetPt3",50000);
  m_cutJetPt4=params.find<double>("cutJetPt4",50000);
  m_cutJetPt5=params.find<double>("cutJetPt5",50000);
  m_cutJetPtWRes=params.find<double>("cutJetPtWRes",40000);

  m_cutRJigsawJetPt=params.find<double>("cutRJigsawJetPt",30000);


  m_cutDeltaPhi=params.find<double>("cutDeltaPhi",0.2);
  m_cutDeltaPhi2=params.find<double>("cutDeltaPhi2",0.4);

  m_cutDeltaPhiQCD=params.find<double>("cutDeltaPhiQCD",0.2);

  // MET significance cuts
  m_cutMETsig[0][0]=0.;
  m_cutMETsig[0][1]=0.;
  m_cutMETsig[0][2]=0.;
  m_cutMETsig[0][3]=0.; 
  m_CRWT_cutMETsig[0][0]=0.;
  m_CRWT_cutMETsig[0][1]=0.;
  m_CRWT_cutMETsig[0][2]=0.;
  m_CRWT_cutMETsig[0][3]=0.;
 
  m_cutMETsig[1][0]=params.find<double>("cutMetSigSR2jLoose",0.);
  m_cutMETsig[1][1]=params.find<double>("cutMetSigSR2jMedium",0.);
  m_cutMETsig[1][2]=params.find<double>("cutMetSigSR2jTight");
  m_CRWT_cutMETsig[1][0]=params.find<double>("cutMetSigCRWT2jLoose",0.);
  m_CRWT_cutMETsig[1][1]=params.find<double>("cutMetSigCRWT2jMedium",0.);
  m_CRWT_cutMETsig[1][2]=params.find<double>("cutMetSigCRWT2jTight");

  m_cutMETsig[3][0]=params.find<double>("cutMetSigSR4jVeryLoose",0.);
  m_cutMETsig[3][1]=params.find<double>("cutMetSigSR4jLoose",0.);
  m_CRWT_cutMETsig[3][0]=params.find<double>("cutMetSigCRWT4jVeryLoose",0.);
  m_CRWT_cutMETsig[3][1]=params.find<double>("cutMetSigCRWT4jLoose",0.);

  // MET/meff(Nj) cuts
  m_cutMEToverMeff[0][0]=0.;
  m_cutMEToverMeff[0][1]=0.;
  m_cutMEToverMeff[0][2]=0.;

  m_cutMEToverMeff[2][2]=params.find<double>("cutMetOverMeffSR3jTight",0.);
  m_CRWT_cutMEToverMeff[2][2]=params.find<double>("cutMetOverMeffCRWT3jTight",0.);

  m_cutMEToverMeff[3][2]=params.find<double>("cutMetOverMeffSR4jMedium",0.);
  m_cutMEToverMeff[3][3]=params.find<double>("cutMetOverMeffSR4jTight",0.);
  m_CRWT_cutMEToverMeff[3][2]=params.find<double>("cutMetOverMeffCRWT4jMedium",0.);
  m_CRWT_cutMEToverMeff[3][3]=params.find<double>("cutMetOverMeffCRWT4jTight",0.);

  m_cutMEToverMeff[4][1]=params.find<double>("cutMetOverMeffSR5jMedium",0.);
  m_CRWT_cutMEToverMeff[4][1]=params.find<double>("cutMetOverMeffCRWT5jMedium",0.);

  m_cutMEToverMeff[5][0]=params.find<double>("cutMetOverMeffSR6jLoose",0.);
  m_cutMEToverMeff[5][1]=params.find<double>("cutMetOverMeffSR6jMedium",0.);
  m_cutMEToverMeff[5][2]=params.find<double>("cutMetOverMeffSR6jTight",0.);
  m_cutMEToverMeff[5][3]=params.find<double>("cutMetOverMeffSR6jVeryTight",0.);
  m_CRWT_cutMEToverMeff[5][0]=params.find<double>("cutMetOverMeffCRWT6jLoose",0.);
  m_CRWT_cutMEToverMeff[5][1]=params.find<double>("cutMetOverMeffCRWT6jMedium",0.);
  m_CRWT_cutMEToverMeff[5][2]=params.find<double>("cutMetOverMeffCRWT6jTight",0.);
  m_CRWT_cutMEToverMeff[5][3]=params.find<double>("cutMetOverMeffCRWT6jVeryTight",0.);

  m_cutMEToverMeff[6][0]=params.find<double>("cutMetOverMeffSR2jWuu",0.);
  m_cutMEToverMeff[7][0]=params.find<double>("cutMetOverMeffSR3jWur",0.);
  m_CRWT_cutMEToverMeff[6][0]=params.find<double>("cutMetOverMeffCRWT2jWuu",0.);
  m_CRWT_cutMEToverMeff[7][0]=params.find<double>("cutMetOverMeffCRWT3jWur",0.);

  // Inclusive meff cuts
  m_cutMeff[0][0]=0.;
  m_cutMeff[0][1]=0.;
  m_cutMeff[0][2]=0.;
  m_CRWT_cutMeff[0][0]=0.;
  m_CRWT_cutMeff[0][1]=0.;
  m_CRWT_cutMeff[0][2]=0.;

  m_cutMeff[1][0]=params.find<double>("cutMeffSR2jLoose",0.);
  m_cutMeff[1][1]=params.find<double>("cutMeffSR2jMedium",0.);
  m_cutMeff[1][2]=params.find<double>("cutMeffSR2jTight",0.);
  m_CRWT_cutMeff[1][0]=params.find<double>("cutMeffCRWT2jLoose",0.);
  m_CRWT_cutMeff[1][1]=params.find<double>("cutMeffCRWT2jMedium",0.);
  m_CRWT_cutMeff[1][2]=params.find<double>("cutMeffCRWT2jTight",0.);

  m_cutMeff[2][2]=params.find<double>("cutMeffSR3jTight",0.);
  m_CRWT_cutMeff[2][2]=params.find<double>("cutMeffCRWT3jTight",0.);

  m_cutMeff[3][0]=params.find<double>("cutMeffSR4jVeryLoose",0.);
  m_cutMeff[3][1]=params.find<double>("cutMeffSR4jLoose",0.);
  m_cutMeff[3][2]=params.find<double>("cutMeffSR4jMedium",0.);
  m_cutMeff[3][3]=params.find<double>("cutMeffSR4jTight",0.);
  m_CRWT_cutMeff[3][0]=params.find<double>("cutMeffCRWT4jVeryLoose",0.);
  m_CRWT_cutMeff[3][1]=params.find<double>("cutMeffCRWT4jLoose",0.);
  m_CRWT_cutMeff[3][2]=params.find<double>("cutMeffCRWT4jMedium",0.);
  m_CRWT_cutMeff[3][3]=params.find<double>("cutMeffCRWT4jTight",0.);

  m_cutMeff[4][1]=params.find<double>("cutMeffSR5jMedium",0.);
  m_CRWT_cutMeff[4][1]=params.find<double>("cutMeffCRWT5jMedium",0.);

  m_cutMeff[5][0]=params.find<double>("cutMeffSR6jLoose",0.);
  m_cutMeff[5][1]=params.find<double>("cutMeffSR6jMedium",0.);
  m_cutMeff[5][2]=params.find<double>("cutMeffSR6jTight",0.);
  m_cutMeff[5][3]=params.find<double>("cutMeffSR6jVeryTight",0.);
  m_CRWT_cutMeff[5][0]=params.find<double>("cutMeffCRWT6jLoose",0.);
  m_CRWT_cutMeff[5][1]=params.find<double>("cutMeffCRWT6jMedium",0.);
  m_CRWT_cutMeff[5][2]=params.find<double>("cutMeffCRWT6jTight",0.);
  m_CRWT_cutMeff[5][3]=params.find<double>("cutMeffCRWT6jVeryTight",0.);

  m_cutMeff[6][0]=params.find<double>("cutMeffSR2jWuu",0.);
  m_cutMeff[7][0]=params.find<double>("cutMeffSR3jWur",0.);

  // Cut map flag telling if a region is used
  m_cutMap[0][0]=params.find<bool>("cutMapSR2jLoose",0.);
  m_cutMap[0][1]=params.find<bool>("cutMapSR2jMedium",0.);
  m_cutMap[0][2]=params.find<bool>("cutMapSR2jTight",0.);

  m_cutMap[1][2]=params.find<bool>("cutMapSR3jTight",0.);

  m_cutMap[2][0]=params.find<bool>("cutMapSR4jVeryLoose",0.);
  m_cutMap[2][1]=params.find<bool>("cutMapSR4jLoose",0.);
  m_cutMap[2][2]=params.find<bool>("cutMapSR4jMedium",0.);
  m_cutMap[2][3]=params.find<bool>("cutMapSR4jTight",0.);

  m_cutMap[3][1]=params.find<bool>("cutMapSR5jMedium",0.);

  m_cutMap[4][0]=params.find<bool>("cutMapSR6jLoose",0.);
  m_cutMap[4][1]=params.find<bool>("cutMapSR6jMedium",0.);
  m_cutMap[4][2]=params.find<bool>("cutMapSR6jTight",0.);
  m_cutMap[4][3]=params.find<bool>("cutMapSR6jVeryTight",0.);

  m_cutMap[5][0]=params.find<bool>("cutMapSR2jWuu",0.);
  m_cutMap[6][0]=params.find<bool>("cutMapSR3jWur",0.);

  // Control region cuts 
  m_cutPhotonPtCRY=params.find<double>("cutPhotonPtCRY",85000.);
  m_cutPhotonIsolCRY=params.find<double>("cutPhotonIsolCRY",5000.);
  m_cutPhotonQualityCRY=params.find<int>("cutPhotonQualityCRY",2);

  m_cutQCDSmearPt=params.find<double>("cutQCDSmearPt",20000.);
  m_cutQCDSmearEta=params.find<double>("cutQCDSmearEta",2.8);
}

