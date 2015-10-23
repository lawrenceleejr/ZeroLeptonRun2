#include "ZeroLeptonRun2/PhysObjProxyUtils.h"
#include "ZeroLeptonRun2/PhysObjProxies.h"
#include "ZeroLeptonRun2/ZeroLeptonNTVars.h"
#include "xAODJet/Jet.h"
#include "TLorentzVector.h"
//#include "OxbridgeKinetics/OKTwoParents.h"

#include "ZeroLeptonRun2/Sphericity.h"

#include <cmath>
#include <iostream>

#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"

PhysObjProxyUtils::PhysObjProxyUtils(bool IsData):
  m_IsData(IsData),
     LAB(nullptr),
   PP(nullptr),
   Pa(nullptr),
   Pb(nullptr),
   Ca(nullptr),
   Cb(nullptr),
   SAV1a(nullptr),
   SAV1b(nullptr),
   SAV2a(nullptr),
   SAV2b(nullptr),
   V1a(nullptr),
   V1b(nullptr),
   V2a(nullptr),
   V2b(nullptr),
   Ia(nullptr),
   Ib(nullptr),
   INV(nullptr),
   VIS(nullptr),
   InvMassJigsaw(nullptr),
   InvRapidityJigsaw(nullptr),
   InvCBJigsaw(nullptr),
   CombPPJigsaw(nullptr),
   CombPaJigsaw(nullptr),
   CombPbJigsaw(nullptr),
   LAB_bkg(nullptr),
  S_bkg(nullptr),
  V_bkg(nullptr),
  I_bkg(nullptr),
  INV_bkg(nullptr),
  VIS_bkg(nullptr),
  InvMass_bkg(nullptr),
  InvRapidity_bkg(nullptr)
{
}

PhysObjProxyUtils::~PhysObjProxyUtils()
{
  delete LAB;
  delete PP;
  delete Pa;
  delete Pb;
  delete Ca;
  delete Cb;
  delete SAV1a;
  delete SAV1b;
  delete SAV2a;
  delete SAV2b;
  delete V1a;
  delete V1b;
  delete V2a;
  delete V2b;
  delete Ia;
  delete Ib;
  delete INV;
  delete VIS;
  delete InvMassJigsaw;
  delete InvRapidityJigsaw;
  delete InvCBJigsaw;
  delete CombPPJigsaw;
  delete CombPaJigsaw;
  delete CombPbJigsaw;
  delete LAB_bkg;
  delete S_bkg;
  delete V_bkg;
  delete I_bkg;
  delete INV_bkg;
  delete VIS_bkg;
  delete InvMass_bkg;
  delete InvRapidity_bkg;
}


void PhysObjProxyUtils::EnergyWeightedTime(const std::vector<JetProxy>& jets,
					   std::vector<float>& time) const
{
  time.resize(5,-999.f);

  double denom = 0.;
  double num = 0.;
  std::size_t maxi = jets.size();
  if ( maxi > 6 ) maxi = 6;
  for ( std::size_t i = 0; i < maxi; i++ ) {
    if ( jets[i].jet() ) {  // skip leptons
      denom = denom + jets[i].E();
      float time = -99999.f;
      jets[i].jet()->getAttribute(xAOD::JetAttribute::Timing,time);
      num = num + jets[i].E() * time;
    }
    if ( i == 1 ) time[0] = num/denom; // 2jets
    else if ( i == 2 ) time[1] = num/denom; // 3jets
    else if ( i == 3 ) time[2] = num/denom; // 4jets
    else if ( i == 4 ) time[3] = num/denom; //5jets
    else if ( i == 5 ) time[4] = num/denom; //6jets
  }

  /*
  std::cout << " jet timing ";
  for ( std::size_t i = 0; i <5; ++i ) std::cout << " " << time[i];
  std::cout << std::endl;
  */
}

double PhysObjProxyUtils::SmallestdPhi(const std::vector<JetProxy>& jets, double met_phi) const
{
  if ( jets.size() < 2 ) return 999.;
  double dphi1 = std::acos(std::cos(jets[0].Phi() - met_phi));
  double dphi2 = std::acos(std::cos(jets[1].Phi() - met_phi));
  double dphi3 = 999.;
  if ( jets.size() > 2 && jets[2].Pt() > 40000. ) {
    dphi3= std::acos(std::cos(jets[2].Phi() - met_phi));
  }
  double min1 = std::min(dphi1,dphi2);
  return std::min(min1,dphi3);
}

double PhysObjProxyUtils::SmallestRemainingdPhi(const std::vector<JetProxy>& jets, double met_phi) const
{
  double remainingDPhi=999;
  double dphiMin=999;
  unsigned int jetcount = 0;
  for ( std::vector<JetProxy>::const_iterator itjet = jets.begin();
	  itjet != jets.end(); ++itjet )
    {
      jetcount++;
      if ( jetcount>3 && itjet->Pt()>40000 ) {
	remainingDPhi = std::acos(std::cos(itjet->Phi() - met_phi));
	dphiMin = std::min(remainingDPhi,dphiMin);
      }
    }
  return dphiMin;
}


double PhysObjProxyUtils::Meff(const std::vector<JetProxy>& jets, size_t njets, double met, double jetPtCut, double extraJetPtCut)
{
  double meff=met;
  if ( jets.size() < njets ) njets=jets.size();
  for(size_t i=0; i<njets; i++)
  {
    if ( i<=3 && jets[i].Pt() > jetPtCut ) meff += jets[i].Pt();
    if ( i>3  && jets[i].Pt() > extraJetPtCut ) meff += jets[i].Pt();
  }
  return meff;
}


void PhysObjProxyUtils::ComputeSphericity(const std::vector<JetProxy>& jets, double & Sp, double & ST, double & Ap)
{

  Sp=-1;
  ST=-1;
  Ap=-1;

  vector<TLorentzVector> v_tlv;

  //prepare vector<TLorentzVector> of jets to use
  for(size_t ijet=0; ijet<jets.size(); ijet++)  {
    if ( jets[ijet].Pt() < 40000. ) continue;
    TLorentzVector jet;
    jet.SetPtEtaPhiM(jets[ijet].Pt(),
		     jets[ijet].Eta(),
		       jets[ijet].Phi(),
		     jets[ijet].M());
    v_tlv.push_back(jet);
  }

  int njet = v_tlv.size();
  if(v_tlv.size() < (size_t)njet || v_tlv.size()==0)return ;


  Sphericity sp; //construct
  sp.SetTLV(v_tlv, njet);
  sp.GetSphericity(Sp, ST, Ap);
};


void PhysObjProxyUtils::RJigsawInit(){
  using namespace RestFrames;

  // cleanup previously computed variables
  if( LAB)               delete  LAB                   ;  LAB               = nullptr;
  if( PP)                delete  PP                    ;  PP                = nullptr;
  if( Pa)                delete  Pa                    ;  Pa                = nullptr;
  if( Pb)                delete  Pb                    ;  Pb                = nullptr;
  if( Ca)                delete  Ca                    ;  Ca                = nullptr;
  if( Cb)                delete  Cb                    ;  Cb                = nullptr;
  if( SAV1a)             delete  SAV1a                 ;  SAV1a             = nullptr;
  if( SAV1b)             delete  SAV1b                 ;  SAV1b             = nullptr;
  if( SAV2a)             delete  SAV2a                 ;  SAV2a             = nullptr;
  if( SAV2b)             delete  SAV2b                 ;  SAV2b             = nullptr;
  if( V1a)               delete  V1a                   ;  V1a               = nullptr;
  if( V1b)               delete  V1b                   ;  V1b               = nullptr;
  if( V2a)               delete  V2a                   ;  V2a               = nullptr;
  if( V2b)               delete  V2b                   ;  V2b               = nullptr;
  if( Ia)                delete  Ia                    ;  Ia                = nullptr;
  if( Ib)                delete  Ib                    ;  Ib                = nullptr;
  if( INV)               delete  INV                   ;  INV               = nullptr;
  if( VIS)               delete  VIS                   ;  VIS               = nullptr;
  if( InvMassJigsaw)     delete  InvMassJigsaw         ;  InvMassJigsaw     = nullptr;
  if( InvRapidityJigsaw) delete  InvRapidityJigsaw     ;  InvRapidityJigsaw = nullptr;
  if( InvCBJigsaw)       delete  InvCBJigsaw           ;  InvCBJigsaw       = nullptr;
  if( CombPPJigsaw)      delete  CombPPJigsaw          ;  CombPPJigsaw      = nullptr;
  if( CombPaJigsaw)      delete  CombPaJigsaw          ;  CombPaJigsaw      = nullptr;
  if( CombPbJigsaw)      delete  CombPbJigsaw          ;  CombPbJigsaw      = nullptr;
  if( LAB_bkg)           delete  LAB_bkg               ;  LAB_bkg           = nullptr;
  if( S_bkg)             delete  S_bkg                 ;  S_bkg             = nullptr;
  if( V_bkg)             delete  V_bkg                 ;  V_bkg             = nullptr;
  if( I_bkg)             delete  I_bkg                 ;  I_bkg             = nullptr;
  if( INV_bkg)           delete  INV_bkg               ;  INV_bkg           = nullptr;
  if( VIS_bkg)           delete  VIS_bkg               ;  VIS_bkg           = nullptr;
  if( InvMass_bkg)       delete  InvMass_bkg           ;  InvMass_bkg       = nullptr;
  if( InvRapidity_bkg)   delete  InvRapidity_bkg       ;  InvRapidity_bkg   = nullptr;

  // RestFrames stuff

  LAB = new LabRecoFrame("LAB","lab");
  PP = new DecayRecoFrame("PP","PP");
  Pa = new DecayRecoFrame("Pa","P_{a}");
  Pb = new DecayRecoFrame("Pb","P_{b}");
  Ca = new DecayRecoFrame("Ca","C_{a}");
  Cb = new DecayRecoFrame("Cb","C_{b}");
  SAV1a = new SelfAssemblingRecoFrame("SAV1a","SA_{V1a}");
  SAV1b = new SelfAssemblingRecoFrame("SAV1b","SA_{V1b}");
  SAV2a = new SelfAssemblingRecoFrame("SAV2a","SA_{V2a}");
  SAV2b = new SelfAssemblingRecoFrame("SAV2b","SA_{V2b}");
  V1a = new VisibleRecoFrame("V1a","V_{1a}");
  V1b = new VisibleRecoFrame("V1b","V_{1b}");
  V2a = new VisibleRecoFrame("V2a","V_{2a}");
  V2b = new VisibleRecoFrame("V2b","V_{2b}");
  Ia = new InvisibleRecoFrame("Ia","I_{a}");
  Ib = new InvisibleRecoFrame("Ib","I_{b}");

  LAB->SetChildFrame(*PP);
  PP->AddChildFrame(*Pa);
  PP->AddChildFrame(*Pb);
  Pa->AddChildFrame(*SAV1a);
  Pb->AddChildFrame(*SAV1b);
  Pa->AddChildFrame(*Ca);
  Pb->AddChildFrame(*Cb);

  SAV1a->AddChildFrame(*V1a);
  SAV1b->AddChildFrame(*V1b);

  Ca->AddChildFrame(*SAV2a);
  Cb->AddChildFrame(*SAV2b);
  Ca->AddChildFrame(*Ia);
  Cb->AddChildFrame(*Ib);

  SAV2a->AddChildFrame(*V2a);
  SAV2b->AddChildFrame(*V2b);

  LAB->InitializeTree();

  INV = new InvisibleGroup("INV","Invisible State Jigsaws");
  INV->AddFrame(*Ia);
  INV->AddFrame(*Ib);

  VIS = new CombinatoricGroup("VIS","Visible Object Jigsaws");
  VIS->AddFrame(*V1a);
  VIS->SetNElementsForFrame(*V1a,0,false);
  VIS->AddFrame(*V1b);
  VIS->SetNElementsForFrame(*V1b,0,false);
  VIS->AddFrame(*V2a);
  VIS->SetNElementsForFrame(*V2a,1,false);
  VIS->AddFrame(*V2b);
  VIS->SetNElementsForFrame(*V2b,1,false);

  InvMassJigsaw = new SetMassInvJigsaw("InvMassJigsaw", "Invisible system mass Jigsaw");
  INV->AddJigsaw(*InvMassJigsaw);

  InvRapidityJigsaw = new SetRapidityInvJigsaw("InvRapidityJigsaw", "Invisible system rapidity Jigsaw");
  INV->AddJigsaw(*InvRapidityJigsaw);
  InvRapidityJigsaw->AddVisibleFrames(LAB->GetListVisibleFrames());

  InvCBJigsaw = new ContraBoostInvJigsaw("InvCBJigsaw","Contraboost invariant Jigsaw");
  INV->AddJigsaw(*InvCBJigsaw);
  InvCBJigsaw->AddVisibleFrames(Pa->GetListVisibleFrames(), 0);
  InvCBJigsaw->AddVisibleFrames(Pb->GetListVisibleFrames(), 1);
  InvCBJigsaw->AddInvisibleFrame(*Ia, 0);
  InvCBJigsaw->AddInvisibleFrame(*Ib, 1);

  CombPPJigsaw = new MinMassesCombJigsaw("CombPPJigsaw","Minimize m _{V_{a,b}} Jigsaw");
  VIS->AddJigsaw(*CombPPJigsaw);
  CombPPJigsaw->AddFrame(*V1a,0);
  CombPPJigsaw->AddFrame(*V1b,1);
  CombPPJigsaw->AddFrame(*V2a,0);
  CombPPJigsaw->AddFrame(*V2b,1);

  CombPaJigsaw = new MinMassesCombJigsaw("C1HEM_JIGSAW","Minimize m _{C_{a}} Jigsaw");
  VIS->AddJigsaw(*CombPaJigsaw);
  CombPaJigsaw->AddFrame(*V1a,0);
  CombPaJigsaw->AddFrame(*V2a,1);

  CombPbJigsaw = new MinMassesCombJigsaw("C2HEM_JIGSAW","Minimize m _{C_{b}} Jigsaw");
  VIS->AddJigsaw(*CombPbJigsaw);
  CombPbJigsaw->AddFrame(*V1b,0);
  CombPbJigsaw->AddFrame(*V2b,1);

  LAB->InitializeAnalysis();

  LAB_bkg = new LabRecoFrame("LAB","lab");
  S_bkg   = new SelfAssemblingRecoFrame("CM","CM");
  V_bkg   = new VisibleRecoFrame("V_bkg","Vis");
  I_bkg   = new InvisibleRecoFrame("I_bkg","Inv");

  LAB_bkg->SetChildFrame(*S_bkg);
  S_bkg->AddChildFrame(*V_bkg);
  S_bkg->AddChildFrame(*I_bkg);

  LAB_bkg->InitializeTree();

  INV_bkg = new InvisibleGroup("INV_bkg","Invisible State Jigsaws");
  INV_bkg->AddFrame(*I_bkg);

  VIS_bkg = new CombinatoricGroup("VIS_bkg","Visible Object Jigsaws");
  VIS_bkg->AddFrame(*V_bkg);
  VIS_bkg->SetNElementsForFrame(*V_bkg,1,false);

  InvMass_bkg = new SetMassInvJigsaw("InvMass_bkg", "Invisible system mass Jigsaw");
  INV_bkg->AddJigsaw(*InvMass_bkg);

  InvRapidity_bkg = new SetRapidityInvJigsaw("InvRapidity_bkg", "Invisible system rapidity Jigsaw");
  INV_bkg->AddJigsaw(*InvRapidity_bkg);
  InvRapidity_bkg->AddVisibleFrames(LAB_bkg->GetListVisibleFrames());

  LAB_bkg->InitializeAnalysis();

  return;

}


void PhysObjProxyUtils::CalculateRJigsawVariables(const std::vector<JetProxy>& jets,
						  Double_t metx,
						  Double_t mety,
						  std::map<TString,float>& RJigsawVariables,
						  Double_t jetPtCut,
						  size_t njetsCut//todo reimplement this
						  ){
  using namespace RestFrames;
  TVector3 ETMiss(metx , mety , 0) ;

  vector<TLorentzVector> Jets;//translate to the code from Chris

  for( size_t jj = 0; jj < std::min(njetsCut, jets.size()); ++jj){
  // for( std::vector<JetProxy>::const_iterator ijet = jets.begin();
  //      ijet != jets.end();
  //      ++ijet
  //      ){
    JetProxy const & ijet = jets.at(jj);
    if( (ijet).Pt() > jetPtCut &&
	(ijet).Eta() < 2.8 ) //todo FIXME hardcode
      Jets.push_back( (ijet) );//we might not need to do this but let's be a bit safer
  }

  // need two jets to play
  if(Jets.size() < 2){
    RJigsawVariables = std::map<TString, float>();
    return;
  }

  LAB->ClearEvent();

  vector<RFKey> jetID;
  for(int i = 0; i < int(Jets.size()); i++){
    jetID.push_back(VIS->AddLabFrameFourVector(Jets[i]));
  }
  INV->SetLabFrameThreeVector(ETMiss);
  if(!LAB->AnalyzeEvent()) cout << "Something went wrong..." << endl;

  float const m_NJet = Jets.size();
  float const m_NJ1a = VIS->GetNElementsInFrame(*V1a);
  float const m_NJ1b = VIS->GetNElementsInFrame(*V1b);
  float const m_NJ2a = VIS->GetNElementsInFrame(*V2a);
  float const m_NJ2b = VIS->GetNElementsInFrame(*V2b);
  float const m_NJa = m_NJ1a+m_NJ2a;
  float const m_NJb = m_NJ1b+m_NJ2b;

  //  if(ETMiss.Mag() < 100. || m_NJet < 2)
  //  return;

  LAB_bkg->ClearEvent();
  double HT = 0.;
  vector<RFKey> jetID_bkg;
  for(int i = 0; i < int(Jets.size()); i++){
    Jets[i].SetPtEtaPhiM(Jets[i].Pt(),0.0,Jets[i].Phi(),Jets[i].M());
    jetID_bkg.push_back(VIS_bkg->AddLabFrameFourVector(Jets[i]));
    HT += Jets[i].Pt();
  }
  INV_bkg->SetLabFrameThreeVector(ETMiss);
  if(!LAB_bkg->AnalyzeEvent()) cout << "Something went wrong..." << endl;

  // QCD clean-up
  TLorentzVector Psib = I_bkg->GetSiblingFrame().GetFourVector(*LAB_bkg);
  TLorentzVector Pmet = I_bkg->GetFourVector(*LAB_bkg);

  float m_Rsib = max(0.,Psib.Vect().Dot(Pmet.Vect().Unit()));
  m_Rsib = m_Rsib / (Pmet.Pt() + m_Rsib);

  TVector3 boostQCD = (Pmet+Psib).BoostVector();
  Psib.Boost(-boostQCD);
  double cosQCD = -1.*Psib.Vect().Unit().Dot(boostQCD.Unit());
  cosQCD = (1.-cosQCD)/2.;
  float const m_deltaQCD = (cosQCD-m_Rsib)/(cosQCD+m_Rsib);

  // signal variables
  TLorentzVector vP_Va = Pa->GetVisibleFourVector(*Pa);
  TLorentzVector vP_Vb = Pb->GetVisibleFourVector(*Pb);
  float const m_MP = (vP_Va.M2()-vP_Vb.M2())/(2.*(vP_Va.E()-vP_Vb.E()));

  TLorentzVector vP_V1aPP = V1a->GetFourVector(*PP);
  TLorentzVector vP_V2aPP = V2a->GetFourVector(*PP);
  TLorentzVector vP_V1bPP = V1b->GetFourVector(*PP);
  TLorentzVector vP_V2bPP = V2b->GetFourVector(*PP);
  TLorentzVector vP_IaPP  = Ia->GetFourVector(*PP);
  TLorentzVector vP_IbPP  = Ib->GetFourVector(*PP);

  TLorentzVector vP_V1aPa = V1a->GetFourVector(*Pa);
  TLorentzVector vP_V2aPa = V2a->GetFourVector(*Pa);
  TLorentzVector vP_IaPa  = Ia->GetFourVector(*Pa);
  TLorentzVector vP_V1bPb = V1b->GetFourVector(*Pb);
  TLorentzVector vP_V2bPb = V2b->GetFourVector(*Pb);
  TLorentzVector vP_IbPb  = Ib->GetFourVector(*Pb);

  float const m_H2PP = (vP_V1aPP + vP_V2aPP + vP_V1bPP + vP_V2bPP).P() + (vP_IaPP+vP_IbPP).P();
  float const m_H3PP = (vP_V1aPP + vP_V2aPP).P() + (vP_V1bPP + vP_V2bPP).P() + (vP_IaPP + vP_IbPP).P();
  float const m_H4PP = (vP_V1aPP + vP_V2aPP).P() + (vP_V1bPP + vP_V2bPP).P() + vP_IaPP.P() + vP_IbPP.P();
  float const m_H6PP = vP_V1aPP.P() + vP_V2aPP.P() + vP_V1bPP.P() + vP_V2bPP.P() + vP_IaPP.P() + vP_IbPP.P();

  float const m_H2Pa = (vP_V1aPa + vP_V2aPa).P() + vP_IaPa.P();
  float const m_H2Pb = (vP_V1bPb + vP_V2bPb).P() + vP_IbPb.P();
  float const m_H3Pa = vP_V1aPa.P() + vP_V2aPa.P() + vP_IaPa.P();
  float const m_H3Pb = vP_V1bPb.P() + vP_V2bPb.P() + vP_IbPb.P();

  float  m_H4Pa = 0.;
  float  m_H4Pb = 0.;
  float  m_H5Pa = 0.;
  float  m_H5Pb = 0.;

  if(m_NJ1a > 1){
    m_H4Pa += SAV1a->GetChildFrame(0).GetMomentum(*Pa);
    m_H4Pa += SAV1a->GetChildFrame(1).GetMomentum(*Pa);
    m_H5Pa += m_H4Pa;
  } else {
    m_H4Pa += vP_V1aPa.P();
    m_H5Pa += vP_V1aPa.P();
  }
  if(m_NJ1b > 1){
    m_H4Pb += SAV1b->GetChildFrame(0).GetMomentum(*Pb);
    m_H4Pb += SAV1b->GetChildFrame(1).GetMomentum(*Pb);
    m_H5Pb += m_H4Pb;
  } else {
    m_H4Pb += vP_V1bPb.P();
    m_H5Pb += vP_V1bPb.P();
  }
  m_H4Pa += vP_V2aPa.P();
  m_H4Pb += vP_V2bPb.P();

  if(m_NJ2a > 1){
    m_H5Pa += SAV2a->GetChildFrame(0).GetMomentum(*Pa);
    m_H5Pa += SAV2a->GetChildFrame(1).GetMomentum(*Pa);
  } else {
    m_H5Pa += vP_V2aPa.P();
  }
  if(m_NJ2b > 1){
    m_H5Pb += SAV2b->GetChildFrame(0).GetMomentum(*Pb);
    m_H5Pb += SAV2b->GetChildFrame(1).GetMomentum(*Pb);
  } else {
    m_H5Pb += vP_V2bPb.P();
  }
  m_H4Pa += vP_IaPa.P();
  m_H5Pa += vP_IaPa.P();
  m_H4Pb += vP_IbPb.P();
  m_H5Pb += vP_IbPb.P();

  TLorentzVector vP_IaCa  = Ia->GetFourVector(*Ca);
  TLorentzVector vP_IbCb  = Ib->GetFourVector(*Cb);

  float const m_H2Ca = 2.*vP_IaCa.P();
  float const m_H2Cb = 2.*vP_IbCb.P();
  float m_H3Ca = 0;
  float m_H3Cb = 0;

  if(m_NJ2a > 1)
    m_H3Ca = vP_IaCa.P()+
      SAV2a->GetChildFrame(0).GetMomentum(*Ca)+
      SAV2a->GetChildFrame(1).GetMomentum(*Ca);
  else
    m_H3Ca = m_H2Ca;

  if(m_NJ2b > 1)
    m_H3Cb = vP_IbCb.P()+
      SAV2b->GetChildFrame(0).GetMomentum(*Cb)+
      SAV2b->GetChildFrame(1).GetMomentum(*Cb);
  else
    m_H3Cb = m_H2Cb;

  double P_P = Pa->GetMomentum(*PP);

  double const m_MPP = 2.*sqrt(P_P*P_P + m_MP*m_MP);
  TVector3 vP_PP = PP->GetFourVector(*LAB).Vect();
  double m_pTCM = vP_PP.Pt();
  double m_pZCM = fabs(vP_PP.Pz());
  float const m_RPT = m_pTCM / (m_pTCM + m_MPP/4.);
  float const m_RPZ = m_pZCM;

  float const m_PP_VisShape = PP->GetVisibleShape();

  float const m_gaminvPP = 2.*m_MP/m_MPP;
  float const m_MDR = m_PP_VisShape*PP->GetMass();

  float const m_cosPP = PP->GetCosDecayAngle();
  float const m_dphiVP = PP->GetDeltaPhiDecayVisible();
  float const m_dphiPPV = PP->GetDeltaPhiBoostVisible();
  float const m_cosP = Pa->GetCosDecayAngle(*Ia);

  // gluino hemishpere variables
  float const m_dphiPCa = Pa->GetDeltaPhiDecayPlanes(*Ca);
  float const m_dphiPCb = Pb->GetDeltaPhiDecayPlanes(*Cb);

  // inside gluino hemisphere variables
  float const m_dphiPV1a = Pa->GetDeltaPhiDecayPlanes(*SAV1a);
  float const m_dphiPV1b = Pb->GetDeltaPhiDecayPlanes(*SAV1b);
  float const m_cosV1a = SAV1a->GetCosDecayAngle();
  float const m_cosV1b = SAV1b->GetCosDecayAngle();
  float const m_dphiCV2a = Ca->GetDeltaPhiDecayPlanes(*SAV2a);
  float const m_dphiCV2b = Cb->GetDeltaPhiDecayPlanes(*SAV2b);
  float const m_cosV2a = SAV2a->GetCosDecayAngle();
  float const m_cosV2b = SAV2b->GetCosDecayAngle();

  // float const m_MET = ETMiss.Pt();
  // float const m_Meff = NTVars_meffInc;
  // float const m_Aplan = NTExtraVars_Ap;
  // float const m_dphi = NTVars_dPhi;
  // float const m_dphiR = NTVars_dPhiR;

  float const m_pT_jet1 = Jets[0].Pt();
  float const m_pT_jet2 = Jets[1].Pt();
  // if(m_NJet >= 3)
  //   m_pT_jet3 = Jets[2].Pt();
  // else
  //   m_pT_jet3 = 0.;
  // if(m_NJet >= 3)
  //   m_pT_jet4 = Jets[3].Pt();
  // else
  //   m_pT_jet4 = 0.;
  // if(m_NJet >= 3)
  //   m_pT_jet5 = Jets[4].Pt();
  // else
  //   m_pT_jet5 = 0.;
  // if(m_NJet >= 3)
  //   m_pT_jet6 = Jets[5].Pt();
  // else
  //   m_pT_jet6 = 0.;

  float const m_pTPP_V1a = V1a->GetTransverseMomentum(*PP);
  float const m_pTPP_V2a = V2a->GetTransverseMomentum(*PP);
  float const m_pTPP_V1b = V1b->GetTransverseMomentum(*PP);
  float const m_pTPP_V2b = V2b->GetTransverseMomentum(*PP);
  float const m_pTPP_Ia = Ia->GetTransverseMomentum(*PP);
  float const m_pTPP_Ib = Ib->GetTransverseMomentum(*PP);

  float const m_pPP_V1a = V1a->GetMomentum(*PP);
  float const m_pPP_V2a = V2a->GetMomentum(*PP);
  float const m_pPP_V1b = V1b->GetMomentum(*PP);
  float const m_pPP_V2b = V2b->GetMomentum(*PP);
  float const m_pPP_Ia = Ia->GetMomentum(*PP);
  float const m_pPP_Ib = Ib->GetMomentum(*PP);

  float  m_pT_jet1a = 0.;
  float  m_pT_jet2a = 0.;
  float  m_pT_jet1b = 0.;
  float  m_pT_jet2b = 0.;

  int N = jetID.size();
  for(int j = 0; j < N; j++){
    RestFrame const& frame = VIS->GetFrame(jetID[j]);
    double pT = VIS->GetLabFrameFourVector(jetID[j]).Pt();

    if(V1a->IsSame(frame) || V2a->IsSame(frame)){
      if(pT > m_pT_jet1a){
	m_pT_jet2a = m_pT_jet1a;
	m_pT_jet1a = pT;
      } else {
	if(pT > m_pT_jet2a) m_pT_jet2a = pT;
      }
    }
    if(V1b->IsSame(frame) || V2b->IsSame(frame)){
      if(pT > m_pT_jet1b){
	m_pT_jet2b = m_pT_jet1b;
	m_pT_jet1b = pT;
      } else {
	if(pT > m_pT_jet2b) m_pT_jet2b = pT;
      }
    }
  }

  float  m_pTPP_jet1a = 0.;
  float  m_pTPP_jet2a = 0.;
  float  m_pTPP_jet1b = 0.;
  float  m_pTPP_jet2b = 0.;

  float m_pPP_jet1a = 0;
  float m_pPP_jet2a = 0;
  float m_pPP_jet1b = 0;
  float m_pPP_jet2b = 0;

  for(int j = 0; j < N; j++){
    RestFrame const& frame = VIS->GetFrame(jetID[j]);
    double pT = PP->GetTransverseMomentum(VIS->GetLabFrameFourVector(jetID[j]));
    double p  = PP->GetFourVector(VIS->GetLabFrameFourVector(jetID[j])).P();

    if(V1a->IsSame(frame) || V2a->IsSame(frame)){
      if(pT > m_pTPP_jet1a){
	m_pTPP_jet2a = m_pTPP_jet1a;
	m_pPP_jet2a  = m_pPP_jet1a;
	m_pTPP_jet1a = pT;
	m_pPP_jet1a  = p;
      } else {
	if(pT > m_pTPP_jet2a){
	  m_pTPP_jet2a = pT;
	  m_pPP_jet2a  = p;
	}
      }
    }
    if(V1b->IsSame(frame) || V2b->IsSame(frame)){
      if(pT > m_pTPP_jet1b){
	m_pTPP_jet2b = m_pTPP_jet1b;
	m_pPP_jet2b  = m_pPP_jet1b;
	m_pTPP_jet1b = pT;
	m_pPP_jet1b  = p;
      } else {
	if(pT > m_pTPP_jet2b){
	  m_pTPP_jet2b = pT;
	  m_pPP_jet2b  = p;
	}
      }
    }
  }

  float PG = Pa->GetMomentum(*PP);
  float MGG = 2.*sqrt(PG*PG + m_MP*m_MP);//todo MG is MP right?
  float gaminvGG = 2.*m_MP/m_MPP;
  float gaminv = PP->GetVisibleShape();
  float beta = sqrt(1.- gaminv*gaminv);
  float betaGG = sqrt(1.- gaminvGG*gaminvGG);

  //*** velocity difference between 'massive' and 'mass-less'
  float const DeltaBetaGG = -(betaGG-beta)/(1.-betaGG*beta);

  // m_weight = GetEventWeight();
  // m_veto = NTVars_veto;
  // m_cleaning = NTVars_cleaning;
  // m_timing = NTVars_timing;


  // if(jetID_R.size() < 2){
  //   RJigsawVariables = std::map<TString, float>();
  //   return;
  // }

  ////////////////////////////////////////////////////////////////////////////////
  // 1st order vars



  RJigsawVariables[ "RJVars_PP_Mass"           ] = m_MPP;
  RJigsawVariables[ "RJVars_PP_InvGamma"       ] = m_PP_VisShape;
  RJigsawVariables[ "RJVars_PP_dPhiBetaR"      ] = PP->GetDeltaPhiBoostVisible();
  RJigsawVariables[ "RJVars_PP_dPhiVis"        ] = PP->GetDeltaPhiVisible();
  RJigsawVariables[ "RJVars_PP_CosTheta"       ] = m_cosPP;
  RJigsawVariables[ "RJVars_PP_dPhiDecayAngle" ] = m_dphiVP ; // I think ...
  RJigsawVariables[ "RJVars_PP_VisShape"       ] = m_PP_VisShape;
  RJigsawVariables[ "RJVars_PP_MDeltaR"        ] = m_MDR;

  RJigsawVariables[ "RJVars_P1_Mass"           ] = Pa->GetMass();
  RJigsawVariables[ "RJVars_P1_CosTheta"       ] = m_cosP; //same as Pa->GetCosDecayAngle(*Ia)
  RJigsawVariables[ "RJVars_P2_Mass"           ] = Pb->GetMass();
  RJigsawVariables[ "RJVars_P2_CosTheta"       ] = Pb->GetCosDecayAngle(*Ib); //I think ...-100;
  RJigsawVariables[ "RJVars_I1_Depth"          ] = Pa->GetFrameDepth   (*Ia);
  RJigsawVariables[ "RJVars_I2_Depth"          ] = Pb->GetFrameDepth   (*Ib);

  RJigsawVariables["RJVars_dphiPV1a"  ] = m_dphiPV1a;
  RJigsawVariables["RJVars_cosV1a"    ] = m_cosV1a;
  RJigsawVariables["RJVars_dphiCV2a"  ] = m_dphiCV2a;
  RJigsawVariables["RJVars_cosV2a"    ] = m_cosV2a;
  RJigsawVariables["RJVars_dphiPV1b" ]  = m_dphiPV1b;
  RJigsawVariables["RJVars_cosV1b"   ]  = m_cosV1b;
  RJigsawVariables["RJVars_dphiCV2b" ]  = m_dphiCV2b;
  RJigsawVariables["RJVars_cosV2b"]	= m_cosV2b;

  RJigsawVariables[ "RJVars_V1_N" ]        = -100; // VIS_R->GetNElementsInFrame(*VS[i]);
  RJigsawVariables[ "RJVars_V2_N" ]        = -100;
  RJigsawVariables[ "RJVars_MP"          ] = m_MP;
  RJigsawVariables[ "RJVars_DeltaBetaGG" ] = DeltaBetaGG;
  RJigsawVariables[ "RJVars_dphiVG"      ] = PP->GetDeltaPhiDecayVisible();

  RJigsawVariables[ "RJVars_QCD_dPhiR"    ] = -100;
  RJigsawVariables[ "RJVars_QCD_Rpt"      ] = m_RPT;
  RJigsawVariables[ "RJVars_QCD_Rsib"    ]  = m_Rsib;
  RJigsawVariables[ "RJVars_QCD_Delta1"   ] = m_deltaQCD;

  RJigsawVariables["RJVars_H2PP"]      = m_H2PP ;
  RJigsawVariables["RJVars_H3PP"]      = m_H3PP;
  RJigsawVariables["RJVars_H4PP"]      = m_H4PP;
  RJigsawVariables["RJVars_H6PP"]      = m_H6PP;
  RJigsawVariables["RJVars_H2Pa"]      = m_H2Pa;
  RJigsawVariables["RJVars_H2Pb"]      = m_H2Pb;
  RJigsawVariables["RJVars_H3Pa"]      = m_H3Pa;
  RJigsawVariables["RJVars_H3Pb"]      = m_H3Pb;
  RJigsawVariables["RJVars_H4Pa"]      = m_H4Pa;
  RJigsawVariables["RJVars_H4Pb"]      = m_H4Pb;
  RJigsawVariables["RJVars_H5Pa"]      = m_H5Pa;
  RJigsawVariables["RJVars_H5Pb"]      = m_H5Pb;
  RJigsawVariables["RJVars_H2Ca"]      = m_H2Ca;
  RJigsawVariables["RJVars_H2Cb"]      = m_H2Cb;
  RJigsawVariables["RJVars_H3Ca"]      = m_H3Ca;
  RJigsawVariables["RJVars_H3Cb"]      = m_H3Cb;
  RJigsawVariables["RJVars_HT4PP"]     = -100 ; //m_HT4PP;
  RJigsawVariables["RJVars_HT6PP"]     = -100 ; //m_HT6PP;
  RJigsawVariables["RJVars_minH3P"]    = -100 ; //m_minH3P;
  RJigsawVariables["RJVars_sangle"]    = -100;
  RJigsawVariables["RJVars_dangle"]    = -100;
  RJigsawVariables["RJVars_ddphiPC"]   = -100;
  RJigsawVariables["RJVars_sdphiPC"]   = -100;
  RJigsawVariables["RJVars_dH2o3P"]    = -100;
  RJigsawVariables["RJVars_RPZ_HT4PP"] = -100;
  RJigsawVariables["RJVars_RPZ_HT6PP"] = -100;

  // end
  ////////////////////////////////////////////////////////////////////////////////



  return;

}




void PhysObjProxyUtils::RazorVariables(const std::vector<JetProxy>& jets,
					 Double_t metx,
					 Double_t mety,
					 double &gaminvRp1 ,
					 double &shatR ,
					 double &mdeltaR ,
					 double &cosptR ,
					 double &Minv2 ,
					 double &Einv ,
					 double & gamma_R,
					 double &dphi_BETA_R ,
					 double &dphi_leg1_leg2 ,
					 double &costhetaR ,
					 double &dphi_BETA_Rp1_BETA_R,
					 double &gamma_Rp1,
					 double &Eleg1,
					 double &Eleg2,
					 double &costhetaRp1)
{
  if ( jets.size() < 2 ) {
    gaminvRp1 = -999.;
    shatR = -999.;
    mdeltaR = -999.;
    cosptR = -999.;
    Minv2 = -999.;
    Einv = -999.;
    gamma_R= -999.;
    dphi_BETA_R = -999.;
    dphi_leg1_leg2 = -999.;
    costhetaR = -999.;
    dphi_BETA_Rp1_BETA_R= -999.;
    gamma_Rp1= -999.;
    Eleg1= -999.;
    Eleg2= -999.;
    costhetaRp1 = -999.;
    return;
  }


  //=============================================================
  // Step 1: make megajet
  //=============================================================
  //This code is adapted from a code from CMS
  //https://twiki.cern.ch/twiki/bin/view/CMSPublic/RazorLikelihoodHowTo

  //To minimize the change in the code, the vector of JetProxy is converted to a vector a TLorentzVector
  std::vector<TLorentzVector> myjets;
  for(size_t ijet=0; ijet<jets.size(); ijet++)
    {
      TLorentzVector jet;
      jet.SetPtEtaPhiM(jets[ijet].Pt(),
		       jets[ijet].Eta(),
		       jets[ijet].Phi(),
		       jets[ijet].M());
      myjets.push_back(jet);
    }



  //build megajets
  //  vector<TLorentzVector> mynewjets;
  TLorentzVector J1, J2;
  //  bool foundGood = false;
  size_t N_comb = 1;
  //for(size_t i = 0; i < myjets.size(); i++)
  for(size_t i = 0; i < myjets.size() && i<15; i++)//code very slow if there are many jets
    {
      N_comb *= 2;
    }
  double M_min = 9999999999999999999999999.0;
  int j_count;
  for(size_t i = 1; i < N_comb-1; i++)
  //  for(size_t i = 1; i+1 < N_comb; i++)
    {
      TLorentzVector  j_temp1, j_temp2;
      int itemp = i;
      j_count = N_comb/2;
      int count = 0;

      while(j_count > 0)
	{

	  TLorentzVector TLorentzJets_count = myjets[count];
	  if(itemp/j_count == 1)
	    {
	      j_temp1 += TLorentzJets_count;
	    }
	  else
	    {
	      j_temp2 += TLorentzJets_count;
	    }

	      itemp -= j_count*(itemp/j_count);
	      j_count /= 2;
	      count++;
	}

      double M_temp = j_temp1.M2()+j_temp2.M2();
      // cout << j_temp1.M2()<< " " << j_temp2.M2()<< " "<< M_temp << " "  << endl;

      // smallest mass
      if(M_temp < M_min)
	{
	  M_min = M_temp;
	  J1 = j_temp1;
	  J2 = j_temp2;
	}
    }

  if(J2.Pt() > J1.Pt())
    {
      TLorentzVector temp = J1;
      J1 = J2;
      J2 = temp;
    }
  //  mynewjets.push_back(J1);
  //  mynewjets.push_back(J2);


  //=============================================================
  // Step 2: compute superrazor variables
  //=============================================================
  //based on code provided privately by L. Lee

  TVector3 MET(metx, mety, 0.0);

  J1.SetVectM(J1.Vect(),0.0);
  J2.SetVectM(J2.Vect(),0.0);

  TVector3 vBETA_z = (1./(J1.E()+J2.E()))*(J1+J2).Vect();
  vBETA_z.SetX(0.0);
  vBETA_z.SetY(0.0);


  //transformation from lab frame to approximate rest frame along beam-axis
  J1.Boost(-vBETA_z);
  J2.Boost(-vBETA_z);

  TVector3 pT_CM = (J1+J2).Vect() + MET;
  pT_CM.SetZ(0.0); //should be redundant...

  Minv2 = (J1+J2).M2();
  Einv = sqrt(MET.Mag2()+Minv2);

  //////////////////////
  // definition of shatR
  //////////////////////
  TLorentzVector J1J2 = J1+J2;

  shatR = sqrt( ((J1J2).E()+Einv)*((J1J2).E()+Einv) - pT_CM.Mag2() );

  TVector3 vBETA_R = (1./sqrt(pT_CM.Mag2() + shatR*shatR))*pT_CM;
  gamma_R = 1./sqrt(1.-vBETA_R.Mag2());


  //transformation from lab frame to R frame
  J1.Boost(-vBETA_R);
  J2.Boost(-vBETA_R);


  dphi_BETA_R = ((J1J2).Vect()).DeltaPhi(vBETA_R);
  dphi_leg1_leg2 = J1.Vect().DeltaPhi(J2.Vect());
  costhetaR =  fabs((J1J2).Vect().Dot(vBETA_R)/((J1J2).Vect().Mag()*vBETA_R.Mag()));


  /////////////
  //
  // R-frame
  //
  /////////////

  TVector3 vBETA_Rp1 = (1./(J1.E()+J2.E()))*(J1.Vect() - J2.Vect());

  ////////////////////////
  // definition of gaminvRp1
  ////////////////////////
  gaminvRp1 = sqrt(1.-vBETA_Rp1.Mag2() );
  gamma_Rp1 = 1/sqrt(1.-vBETA_Rp1.Mag2() );
  dphi_BETA_Rp1_BETA_R = vBETA_Rp1.DeltaPhi(vBETA_R);

  //transformation from R frame to R+1 frames
  J1.Boost(-vBETA_Rp1);
  J2.Boost(vBETA_Rp1);
  //////////////
  //
  // R+1-frames
  //
  //////////////

  ///////////////////////
  // definition of mdeltaR
  ////////////////////////
  mdeltaR = J1.E()+J2.E();

  ///////////////////////
  // definition of cosptR
  ////////////////////////
  cosptR = pT_CM.Mag()/sqrt(pT_CM.Mag2()+mdeltaR * mdeltaR);


  Eleg1 = J1.E();
  Eleg2 = J2.E();
  costhetaRp1 = fabs(J1.Vect().Dot(vBETA_Rp1)/(J1.Vect().Mag()*vBETA_Rp1.Mag()));





 }

/*
double PhysObjProxyUtils::MT2(const std::vector<JetProxy>& jets,const TVector2& MissingET) const
{
  double mT2=0;
  if (jets.size()<2) return -1;

  std::vector<TLorentzVector> jets_tmp;
  for (size_t i = 0 ; i <jets.size(); i++ ) {
    jets_tmp.push_back(TLorentzVector(*dynamic_cast<const TLorentzVector*>(&(jets[i]))));
  }

  OxbridgeKinetics::OKTwoParents ok2p;
  TVector2 ptmiss(MissingET.X(),MissingET.Y());
  double Minvis = 0.;
  ok2p.addVis(jets_tmp[0],1);
  ok2p.addVis(jets_tmp[1],2);
  ok2p.setPtMiss(ptmiss);
  ok2p.setMinvis(Minvis);

  mT2 = ok2p.calcM2T();

  return mT2;
}
*/

void PhysObjProxyUtils::GetAlphaISRVar(const std::vector<JetProxy>& jets, double met, std::vector<double>& alpha_vec) const
{
  alpha_vec.clear();
  for (size_t ijet=0; ijet<jets.size() ; ijet++) {
    if ( jets[ijet].Pt() < 50000.) continue;
    double alpha = std::min(jets[ijet].Pt(),met)/std::max(jets[ijet].Pt(),met);
    alpha_vec.push_back(alpha);
  }
  return;
}

void PhysObjProxyUtils::GetMinPtDistinctionISR(const std::vector<JetProxy>& jets, std::vector<double>& minPtDistinction_vec) const
{
  minPtDistinction_vec.clear();
  for (size_t isrcandidate=0; isrcandidate<jets.size() ; isrcandidate++) {
    if ( jets[isrcandidate].Pt() < 50000.) continue;
    double minDist=99.;
    for (size_t jjet=0; jjet<jets.size() ; jjet++) {
      if (jjet==isrcandidate) continue;
      double max2jet_pt=std::max(jets[jjet].Pt(),jets[isrcandidate].Pt());
      double min2jet_pt=std::min(jets[jjet].Pt(),jets[isrcandidate].Pt());
      double fraction=max2jet_pt/min2jet_pt;
      minDist=std::min(minDist,fraction);
    }
    minPtDistinction_vec.push_back(minDist);
  }
  return;

}

void PhysObjProxyUtils::GetMinDeltaFraction(const std::vector<JetProxy>& jets, std::vector<double>& minDeltaFrac_vec) const
{
  minDeltaFrac_vec.clear();
  for (size_t isrcandidate=0; isrcandidate<jets.size() ; isrcandidate++) {
    if (jets[isrcandidate].Pt() < 50000.) continue;
    double delta_isrcand = jets[isrcandidate].M()/jets[isrcandidate].Pt();
    double minDeltaFrac=99.;
    for (size_t jjet=0; jjet<jets.size() ; jjet++) {
      if (jjet==isrcandidate) continue;
      double delta_jjet = jets[jjet].M()/jets[jjet].Pt();
      double mindelta = std::min(delta_isrcand,delta_jjet);
      double maxdelta = std::max(delta_isrcand,delta_jjet);
      double fraction = maxdelta/mindelta;
      minDeltaFrac = std::min(minDeltaFrac,fraction);
    }
    minDeltaFrac_vec.push_back(minDeltaFrac);
  }
  return;
}


void PhysObjProxyUtils::GetMinRapidityGap(const std::vector<JetProxy>& jets, std::vector<double>& minRapGap_vec) const
{
  minRapGap_vec.clear();
  for (size_t isrcandidate=0; isrcandidate<jets.size() ; isrcandidate++) {
    if (jets[isrcandidate].Pt() < 50000.) continue;
    double minRapGap=99.;
    for (size_t jjet=0; jjet<jets.size() ; jjet++) {
      if (jjet==isrcandidate) continue;
      double diffrap = std::abs(jets[isrcandidate].Eta()-jets[jjet].Eta());
      minRapGap = std::min(minRapGap,diffrap);
    }
    minRapGap_vec.push_back(minRapGap);
  }
  return;
}

void PhysObjProxyUtils::GetMaxRapidityOtherJets(const std::vector<JetProxy>& jets, std::vector<double>& maxRapOtherJets_vec) const
{
  maxRapOtherJets_vec.clear();
  for (size_t isrcandidate=0; isrcandidate<jets.size() ; isrcandidate++) {
    if (jets[isrcandidate].Pt() < 50000.) continue;
    double maxRap=-99.;
    for (size_t jjet=0; jjet<jets.size() ; jjet++) {
      if (jjet==isrcandidate) continue;
      maxRap = std::max(maxRap,jets[jjet].Eta());
    }
    maxRapOtherJets_vec.push_back(maxRap);
  }
  return;
}

void PhysObjProxyUtils::GetdPhiJetMet(const std::vector<JetProxy>& jets, double met_phi, std::vector<double>& dPhiJetMet_vec) const
{
  dPhiJetMet_vec.clear();
  for (size_t isrcandidate=0; isrcandidate<jets.size() ; isrcandidate++) {
    if (jets[isrcandidate].Pt() < 50000.) continue;
    double dphi = std::abs(jets[isrcandidate].Phi()-met_phi);
    dPhiJetMet_vec.push_back(dphi);
  }
  return;
}

void PhysObjProxyUtils::GetISRJet(const std::vector<JetProxy>& jets,
				  std::vector<size_t>& isr_jet_indices,
				  double met,
				  double phi_met,
				  std::string signal,
				  bool usealpha) const
{
  using std::cout;
  using std::endl;

  std::vector<size_t> tmp_isr_jet_indices;
  std::vector<bool> maxpt_crit_vector;
  std::vector<bool> maxdelta_crit_vector;
  std::vector<bool> rapid_crit_vector;
  if (!(signal=="gluino" || signal=="squark")) {
    cout << " [GetISRJet]: a signal type other than 'squark' or 'gluino' is requested. Please use one of those two. " << endl;
    cout << " Exiting without tagging an ISR jet " << endl;
  }


  for (size_t isrcandidate=0; isrcandidate<jets.size() ; isrcandidate++) {
    if (jets[isrcandidate].Pt() < 50000.) continue;  // Only consider jets with Pt > 50 GeV
    double delta_isrjet = jets[isrcandidate].M()/jets[isrcandidate].Pt();
    bool pretag_crit=true;
    bool maxpt_crit=true;
    bool rapid_crit=true;
    bool maxdelta_crit=true;
    for (size_t jjet=0; jjet<jets.size() ; jjet++) {
      if (jjet==isrcandidate) continue;

      double max2jet_pt = std::max(jets[jjet].Pt(),jets[isrcandidate].Pt());
      double min2jet_pt = std::min(jets[jjet].Pt(),jets[isrcandidate].Pt());
      double diffrap = std::abs(jets[isrcandidate].Eta()-jets[jjet].Eta());
      double delta_jjet = jets[jjet].M()/jets[jjet].Pt();
      double max2jet_delta = std::max(delta_isrjet,delta_jjet);
      double min2jet_delta = std::min(delta_isrjet,delta_jjet);

      // pretag criteria:

      if (max2jet_pt/min2jet_pt<=2.0) maxpt_crit = false;
      if (diffrap<=1.0) rapid_crit = false;
      if (max2jet_delta/min2jet_delta<=1.5) maxdelta_crit = false;
    }
    if (signal=="gluino") {
      if (maxpt_crit || rapid_crit || maxdelta_crit) pretag_crit = true;
    }
    if (signal=="squark") {
      if ( rapid_crit) pretag_crit = true;
    }
    if (!pretag_crit) continue;

    // Only isr candidate jets which have survived pre-tag conditions
    bool jjetrapidity_crit = true;
    bool isr_rapidity_diff_crit = true;
    for (size_t jjet = 0; jjet<jets.size() ; jjet++) {
      if (jjet==isrcandidate) continue;
      if (std::abs(jets[jjet].Eta())>=2.0) jjetrapidity_crit = false;
      if (std::abs(jets[isrcandidate].Eta()-jets[jjet].Eta())<=0.5) isr_rapidity_diff_crit = false;
    }
    bool isrjetrapidity_crit = false;
    if (std::abs(jets[isrcandidate].Eta())>1.0) isrjetrapidity_crit = true;
    bool deltaphi_crit = false;
    double dphi = std::abs(jets[isrcandidate].Phi()-phi_met);
    if (dphi>=TMath::Pi()) dphi = 2*TMath::Pi()-dphi;
    if (dphi>2.0) deltaphi_crit = true;
    bool alpha_crit = false;
    double alpha = std::min(jets[isrcandidate].Pt(),met)/std::max(jets[isrcandidate].Pt(),met);
    if (alpha>0.4) alpha_crit = true;

    // tagging:
    if (signal=="gluino") {
      if (jjetrapidity_crit && isrjetrapidity_crit && isr_rapidity_diff_crit && deltaphi_crit && ( ( alpha_crit && usealpha ) || (!usealpha))) {
	tmp_isr_jet_indices.push_back(isrcandidate);

	maxpt_crit_vector.push_back(maxpt_crit);
	maxdelta_crit_vector.push_back(maxdelta_crit);
	rapid_crit_vector.push_back(rapid_crit);
      }
    }
    if (signal=="squark") {
      if (jjetrapidity_crit && isrjetrapidity_crit && deltaphi_crit && ( ( alpha_crit && usealpha ) || (!usealpha))) {
	tmp_isr_jet_indices.push_back(isrcandidate);
	rapid_crit_vector.push_back(rapid_crit);
      }
    }
  }


  if (tmp_isr_jet_indices.size() > 1) {
    if (signal=="gluino") {
      std::vector<int> which_criteria;
      std::vector<size_t> whichjet_crit1;
      std::vector<size_t> whichjet_crit2;
      std::vector<size_t> whichjet_crit3;
      for (size_t ijet=0; ijet <tmp_isr_jet_indices.size() ; ijet++) {
	if (maxpt_crit_vector[ijet]==true) whichjet_crit1.push_back(tmp_isr_jet_indices[ijet]);
	if (maxpt_crit_vector[ijet]==false && rapid_crit_vector[ijet]==true)  whichjet_crit1.push_back(tmp_isr_jet_indices[ijet]);
	if (maxpt_crit_vector[ijet]==false && rapid_crit_vector[ijet]==false && maxdelta_crit_vector[ijet]==true)  whichjet_crit1.push_back(tmp_isr_jet_indices[ijet]);
      }
      if (whichjet_crit1.size()==1) isr_jet_indices.push_back(whichjet_crit1[0]);
      else if (whichjet_crit1.size()==0 && whichjet_crit2.size()==1 ) isr_jet_indices.push_back(whichjet_crit2[0]);
      else if (whichjet_crit1.size()>1 ) {
	double maxpt = 0;
	size_t thisjet = 0;
	for (size_t ijet=0; ijet < whichjet_crit1.size() ; ijet++){
	  if (jets[whichjet_crit1[ijet]].Pt() > maxpt){
	    maxpt =  jets[whichjet_crit1[ijet]].Pt();
	    thisjet = ijet;
	  }
	}
	isr_jet_indices.push_back(whichjet_crit1[thisjet]);

	return;
      }

      else if (whichjet_crit1.size()==0 && whichjet_crit2.size()>1 )   return;
      else if (whichjet_crit1.size()==0 && whichjet_crit2.size()==0 && whichjet_crit3.size()>1 ) {
	double maxdelta = 0;
	size_t thisjet = 0;
	for (size_t ijet=0; ijet < whichjet_crit3.size() ; ijet++){
	  if (jets[whichjet_crit3[ijet]].M()/jets[whichjet_crit3[ijet]].Pt()> maxdelta) {
	    maxdelta = jets[whichjet_crit3[ijet]].M()/jets[whichjet_crit3[ijet]].Pt();
	    thisjet = ijet;
	  }
	}
	isr_jet_indices.push_back(whichjet_crit3[thisjet]);

	return;
      }
    }

    else if (signal=="squark") {
      // if more than 1 jet passes these selections: none of them is tagged as an ISR jet!

      return;
    }
  }
  else if (tmp_isr_jet_indices.size() == 1) {
    isr_jet_indices = tmp_isr_jet_indices;
  }
  return;
}


bool PhysObjProxyUtils::CosmicMuon(const std::vector<MuonProxy>& muons) const
{
  for ( size_t i = 0; i < muons.size(); ++i ) {
    if ( muons[i].isCosmic() ) return true;
  }
  return false;
}


bool PhysObjProxyUtils::isbadMETmuon(const std::vector<MuonProxy>& muons,
				     float MET, const TVector2& MissingET) const
{
  bool isbadmetmuon=false;
  TVector2 MissingETMuon(0,0);
  for ( size_t iMu = 0; iMu < muons.size(); ++iMu ) {
    MissingETMuon -= muons[iMu].Vect().XYvector();
  }
  double METMuon = MissingETMuon.Mod();

  double MET_muon_ratio = METMuon/MET*std::cos(MissingETMuon.Phi()-MissingET.Phi()) ;
  if(MET_muon_ratio>=0.5) isbadmetmuon=true;

  return isbadmetmuon;
}


bool PhysObjProxyUtils::badTileVeto(const std::vector<JetProxy>& jets, const TVector2& MissingET) const
{
  bool isDeadTile=false;
  for ( std::vector<JetProxy>::const_iterator itjet = jets.begin();
	itjet != jets.end(); itjet++ ) {
    if ( ! (*itjet).jet() ) continue;
    double jet_pt = (*itjet).Pt();
    if(jet_pt<40000.) continue;
    double jet_phi = (*itjet).Phi();
    double  jet_BCH_CORR_JET = 0.;
    (*itjet).jet()->getAttribute(xAOD::JetAttribute::BchCorrJet,jet_BCH_CORR_JET);
    double MET_phi = MissingET.Phi();
    if(std::acos(std::cos( jet_phi-MET_phi )  )<0.3 && jet_BCH_CORR_JET>0.05) isDeadTile=true;
  }
  return isDeadTile ;
}

bool PhysObjProxyUtils::chfVeto(const std::vector<JetProxy>& jets) const
{
  bool shouldbecleaned=false;
  for ( std::size_t i = 0; i < std::min((std::size_t)2,jets.size()); ++i ) {
    if ( !jets[i].jet() ) continue; // skip lepton disguised as jet

    std::vector<float> sumPtTrk;
    jets[i].jet()->getAttribute(xAOD::JetAttribute::SumPtTrkPt500,sumPtTrk); // FIXME or SumPtTrkPt1000 ??
    float chf = sumPtTrk[0]/jets[i].Pt();

    float emf = 0.;
    jets[i].jet()->getAttribute(xAOD::JetAttribute::EMFrac,emf);

    if ( jets[i].Pt()>100000. && chf<0.02 && std::abs(jets[i].Eta())<2. ) shouldbecleaned=true;

    if ( jets[i].Pt()>100000. && chf<0.05 && std::abs(jets[i].Eta())<2. && emf>0.9) shouldbecleaned=true;
  }
  return shouldbecleaned;
}

bool PhysObjProxyUtils::chfTileVeto(const std::vector<JetProxy>& jets) const
{
  bool shouldbecleaned=false;
  for ( std::size_t i = 0; i < std::min((std::size_t)2,jets.size()); ++i ) {
    bool isIn=false;
    if ( jets[i].Eta()<-0.6 && jets[i].Eta()>-1.0 &&
	 jets[i].Phi()<-0.6 && jets[i].Phi()>-1.0 ) isIn=true;

    std::vector<float> sumPtTrk;
    jets[i].jet()->getAttribute(xAOD::JetAttribute::SumPtTrkPt500,sumPtTrk); // FIXME or SumPtTrkPt1000 ??
    float chf = sumPtTrk[0]/jets[i].Pt();

    float emf = 0.;
    jets[i].jet()->getAttribute(xAOD::JetAttribute::EMFrac,emf);

    if ( isIn && jets[i].Pt()>100000. &&  chf<0.20 && emf<0.3 ) shouldbecleaned=true;
  }
  return shouldbecleaned;
}


void PhysObjProxyUtils::FillNTExtraVars(NTExtraVars& extrantv,
					double MET_Track,
					double MET_Track_phi,
					double mT2,
					double mT2_noISR,
					double Ap)
{
  extrantv.Reset();
  extrantv.mettrack = MET_Track * 0.001;
  extrantv.mettrack_phi = MET_Track_phi;
  extrantv.mT2=mT2 * 0.001;
  extrantv.mT2_noISR=mT2_noISR * 0.001;
  extrantv.Ap=Ap;
}

// void PhysObjProxyUtils::FillNTExtraVarsTriggerBits(NTExtraVars& extrantv,
// 						   long triggers){
//   extrantv.triggers = triggers;
// }

void PhysObjProxyUtils::FillNTRJigsawVars(NTRJigsawVars& rjigsawntv,
              std::map<TString,float> & RJigsawVariables
          )
{
  // rjigsawntv.Reset();

  //std::cout << "In filling function----------------" << std::endl;
  //std::cout << RJigsawVariables["RJVars_P_0_CosTheta"] << " -----------------------" << std::endl;
 rjigsawntv.RJVars_PP_Mass           = RJigsawVariables[ "RJVars_PP_Mass"           ]; //] = m_MPP];
 rjigsawntv.RJVars_PP_InvGamma       = RJigsawVariables[ "RJVars_PP_InvGamma"       ]; //] = m_PP_VisShape];
 rjigsawntv.RJVars_PP_dPhiBetaR      = RJigsawVariables[ "RJVars_PP_dPhiBetaR"      ]; //] = PP->GetDeltaPhiBoostVisible()];
 rjigsawntv.RJVars_PP_dPhiVis        = RJigsawVariables[ "RJVars_PP_dPhiVis"        ]; //] = PP->GetDeltaPhiVisible()];
 rjigsawntv.RJVars_PP_CosTheta       = RJigsawVariables[ "RJVars_PP_CosTheta"       ]; //] = m_cosPP];
 rjigsawntv.RJVars_PP_dPhiDecayAngle = RJigsawVariables[ "RJVars_PP_dPhiDecayAngle" ]; //] = m_dphiVP ]; // I think ...
 rjigsawntv.RJVars_PP_VisShape       = RJigsawVariables[ "RJVars_PP_VisShape"       ]; //] = m_PP_VisShape];
 rjigsawntv.RJVars_PP_MDeltaR        = RJigsawVariables[ "RJVars_PP_MDeltaR"        ]; //] = m_MDR];
 rjigsawntv.RJVars_P1_Mass           = RJigsawVariables[ "RJVars_P1_Mass"           ];
 rjigsawntv.RJVars_P1_CosTheta       = RJigsawVariables[ "RJVars_P1_CosTheta"       ]; //] = Pa->GetMass();
 rjigsawntv.RJVars_P2_Mass           = RJigsawVariables[ "RJVars_P2_Mass"           ]; //] = m_cosP; //same as Pa->GetCosDecayAngle(*Ia)
 rjigsawntv.RJVars_P2_CosTheta       = RJigsawVariables[ "RJVars_P2_CosTheta"       ]; //] = Pa->GetMass();
 rjigsawntv.RJVars_I1_Depth          = RJigsawVariables[ "RJVars_I1_Depth"          ]; //] = Pb->GetCosDecayAngle(*Ib); //I think ...-100;
 rjigsawntv.RJVars_I2_Depth          = RJigsawVariables[ "RJVars_I2_Depth"          ]; //] = Pa->GetFrameDepth   (*Ia);
 rjigsawntv.RJVars_V1_N              = RJigsawVariables[ "RJVars_V1_N" ]; //] = Pb->GetFrameDepth   (*Ib);
 rjigsawntv.RJVars_V2_N              = RJigsawVariables[ "RJVars_V2_N" ];
 rjigsawntv.RJVars_MP                = RJigsawVariables["RJVars_MP"];

 rjigsawntv.RJVars_dphiPV1a          = RJigsawVariables["RJVars_dphiPV1a"  ]; //] = m_dphiPV1a;
 rjigsawntv.RJVars_cosV1a            = RJigsawVariables["RJVars_cosV1a"    ]; //] = m_cosV1a;
 rjigsawntv.RJVars_dphiCV2a          = RJigsawVariables["RJVars_dphiCV2a"  ]; //] = m_dphiCV2a;
 rjigsawntv.RJVars_cosV2a            = RJigsawVariables["RJVars_cosV2a"    ]; //] = m_cosV2a;
 rjigsawntv.RJVars_dphiPV1b          = RJigsawVariables["RJVars_dphiPV1b" ] ; // = m_dphiPV1b;
 rjigsawntv.RJVars_cosV1b            = RJigsawVariables["RJVars_cosV1b"   ] ; // = m_cosV1b;
 rjigsawntv.RJVars_dphiCV2b          = RJigsawVariables["RJVars_dphiCV2b" ]  ; //= m_dphiCV2b;
 rjigsawntv.RJVars_cosV2b            = RJigsawVariables["RJVars_cosV2b"]; //	= m_cosV2b;

 rjigsawntv.RJVars_DeltaBetaGG = RJigsawVariables["RJVars_DeltaBetaGG" ]; //] = DeltaBetaGG;
 rjigsawntv.RJVars_dphiVG     = RJigsawVariables["RJVars_dphiVG"      ]; //] = PP->GetDeltaPhiDecayVisible();
 rjigsawntv.RJVars_QCD_dPhiR  = RJigsawVariables["RJVars_QCD_dPhiR"    ]; //] = -100;
 rjigsawntv.RJVars_QCD_Rpt    = RJigsawVariables["RJVars_QCD_Rpt"      ]; //] = m_RPT;
 rjigsawntv.RJVars_QCD_Rsib   = RJigsawVariables["RJVars_QCD_Rsib"    ];
 rjigsawntv.RJVars_QCD_Delta1 = RJigsawVariables["RJVars_QCD_Delta1"   ]; //] = m_deltaQCD;

 rjigsawntv.RJVars_H2PP      = RJigsawVariables["RJVars_H2PP"];
 rjigsawntv.RJVars_H3PP      = RJigsawVariables["RJVars_H3PP"];
 rjigsawntv.RJVars_H4PP      = RJigsawVariables["RJVars_H4PP"];
 rjigsawntv.RJVars_H6PP      = RJigsawVariables["RJVars_H6PP"];
 rjigsawntv.RJVars_H2Pa      = RJigsawVariables["RJVars_H2Pa"];
 rjigsawntv.RJVars_H2Pb      = RJigsawVariables["RJVars_H2Pb"];
 rjigsawntv.RJVars_H3Pa      = RJigsawVariables["RJVars_H3Pa"];
 rjigsawntv.RJVars_H3Pb      = RJigsawVariables["RJVars_H3Pb"];
 rjigsawntv.RJVars_H4Pa      = RJigsawVariables["RJVars_H4Pa"];
 rjigsawntv.RJVars_H4Pb      = RJigsawVariables["RJVars_H4Pb"];
 rjigsawntv.RJVars_H5Pa      = RJigsawVariables["RJVars_H5Pa"];
 rjigsawntv.RJVars_H5Pb      = RJigsawVariables["RJVars_H5Pb"];
 rjigsawntv.RJVars_H2Ca      = RJigsawVariables["RJVars_H2Ca"];
 rjigsawntv.RJVars_H2Cb      = RJigsawVariables["RJVars_H2Cb"];
 rjigsawntv.RJVars_H3Ca      = RJigsawVariables["RJVars_H3Ca"];
 rjigsawntv.RJVars_H3Cb      = RJigsawVariables["RJVars_H3Cb"];
 rjigsawntv.RJVars_HT4PP     = RJigsawVariables["RJVars_HT4PP"];
 rjigsawntv.RJVars_HT6PP     = RJigsawVariables["RJVars_HT6PP"];
 rjigsawntv.RJVars_minH3P    = RJigsawVariables["RJVars_minH3P"];
 rjigsawntv.RJVars_sangle    = RJigsawVariables["RJVars_sangle"];
 rjigsawntv.RJVars_dangle    = RJigsawVariables["RJVars_dangle"];
 rjigsawntv.RJVars_ddphiPC   = RJigsawVariables["RJVars_ddphiPC"];
 rjigsawntv.RJVars_sdphiPC   = RJigsawVariables["RJVars_sdphiPC"];
 rjigsawntv.RJVars_dH2o3P    = RJigsawVariables["RJVars_dH2o3P"];
 rjigsawntv.RJVars_RPZ_HT4PP = RJigsawVariables["RJVars_RPZ_HT4PP"];
 rjigsawntv.RJVars_RPZ_HT6PP = RJigsawVariables["RJVars_RPZ_HT6PP"];

}




void PhysObjProxyUtils::FillNTVars(NTVars& ntv,
				   unsigned int RunNumber,
				   unsigned int EventNumber,
				   unsigned int LumiBlockNumber,
				   unsigned int veto,
				   float weight,
				   std::vector<float>& normWeight,
				   std::vector<float>& pileupWeight,
				   float genWeight,
				   float ttbarWeightHT,
				   float ttbarWeightPt2,
				   float ttbarAvgPt,
				   float WZweight,
				   std::vector<float>& bTagWeight,
				   std::vector<float>& cTagWeight,
				   int nBJet,
				   int nCJet,
				   double MissingEt,
				   double METPhi,
				   double* Meff,
				   double meffincl,
				   double minDphi,
				   double RemainingminDPhi,
				   const std::vector<JetProxy>& good_jets,
				   int hardproc,
				   unsigned int cleaning,
				   float timing,
				   const std::vector<float>& jetSmearSystW,
				   const std::vector<float>* flaggedtau,
				   float tauMt,
				   float SherpaBugMET,
				   bool isTruth,
				   std::vector<TauProxy> baseline_taus,
				   std::vector<TauProxy> signal_taus)
{
  ntv.Reset();

  ntv.RunNumber = RunNumber;
  ntv.EventNumber = EventNumber;
  ntv.LumiBlockNumber = LumiBlockNumber;
  ntv.veto = veto;
  ntv.weight = weight;
  ntv.pileupWeight     = (pileupWeight.size() >= 1) ? pileupWeight.at(0) : 1.;
  ntv.pileupWeightUp   = (pileupWeight.size() >= 2) ? pileupWeight.at(1) : 1.;
  ntv.pileupWeightDown = (pileupWeight.size() >= 3) ? pileupWeight.at(2) : 1.;
  ntv.genWeight = genWeight;
  ntv.ttbarWeightHT = ttbarWeightHT;
  ntv.ttbarWeightPt2 = ttbarWeightPt2;
  ntv.ttbarAvgPt = ttbarAvgPt * 0.001;
  ntv.WZweight = WZweight;
  ntv.bTagWeight     = (bTagWeight.size() >= 1) ? bTagWeight.at(0) : 1.;
  ntv.bTagWeightBUp   = (bTagWeight.size() >= 5) ? bTagWeight.at(4) : 1.;
  ntv.bTagWeightBDown = (bTagWeight.size() >= 2) ? bTagWeight.at(1) : 1.;
  ntv.bTagWeightCUp   = (bTagWeight.size() >= 6) ? bTagWeight.at(5) : 1.;
  ntv.bTagWeightCDown = (bTagWeight.size() >= 3) ? bTagWeight.at(2) : 1.;
  ntv.bTagWeightLUp   = (bTagWeight.size() >= 7) ? bTagWeight.at(6) : 1.;
  ntv.bTagWeightLDown = (bTagWeight.size() >= 4) ? bTagWeight.at(3) : 1.;
  ntv.cTagWeight     = (cTagWeight.size() >= 1) ? cTagWeight.at(0) : 1.;
  ntv.cTagWeightBUp   = (cTagWeight.size() >= 5) ? cTagWeight.at(4) : 1.;
  ntv.cTagWeightBDown = (cTagWeight.size() >= 2) ? cTagWeight.at(1) : 1.;
  ntv.cTagWeightCUp   = (cTagWeight.size() >= 6) ? cTagWeight.at(5) : 1.;
  ntv.cTagWeightCDown = (cTagWeight.size() >= 3) ? cTagWeight.at(2) : 1.;
  ntv.cTagWeightLUp   = (cTagWeight.size() >= 7) ? cTagWeight.at(6) : 1.;
  ntv.cTagWeightLDown = (cTagWeight.size() >= 4) ? cTagWeight.at(3) : 1.;
  ntv.nBJet = nBJet;
  ntv.nCJet = nCJet;
  ntv.MET = MissingEt * 0.001;
  ntv.METPhi = METPhi;
  ntv.deltaPhi = minDphi;
  ntv.deltaPhiRemaining=RemainingminDPhi;

  ntv.MeffIncl= meffincl * 0.001;
  ntv.normWeight=(normWeight.size() >= 1) ? normWeight.at(0) : 1.;
  ntv.normWeightUp=(normWeight.size() >= 2) ? normWeight.at(1) : 1.;
  ntv.normWeightDown=(normWeight.size() >= 3) ? normWeight.at(2) : 1.;

  ntv.hardproc=hardproc;
  ntv.cleaning=cleaning;
  ntv.timing=timing;

  ntv.SherpaBugMET = SherpaBugMET * 0.001;

  // pdgid of incoming partons
  // FIXME
  /*
  if (!m_isData) {
    if (m_data->mcevt_pdf_id1 && m_data->mcevt_pdf_id1->size())
      ntv.pdfId1 = m_data->mcevt_pdf_id1->at(0);
    if (m_data->mcevt_pdf_id2 && m_data->mcevt_pdf_id2->size())
      ntv.pdfId2 = m_data->mcevt_pdf_id2->at(0);
  }
  */


  ntv.Njet = 0;
  TLorentzVector jet1TLV;
  TLorentzVector jet2TLV;
  for ( size_t jet0=0; jet0<good_jets.size(); jet0++)
  {
    const JetProxy& thisjet = good_jets[jet0];
    float pt = thisjet.Pt();
    float eta = (pt > 0.) ? thisjet.Eta() : 0.;
    float phi = (pt > 0.) ? thisjet.Phi() : 0.;
    float m = (pt > 0.) ? thisjet.M() : 0.;

    float emf = 0.;
    float chf = 1.;
    double jetbtag = -999.;  // MV1 not available on xAOD
    if ( thisjet.jet() && !isTruth ) {
      thisjet.jet()->getAttribute(xAOD::JetAttribute::EMFrac,emf);
      std::vector<float> sumPtTrk;
      thisjet.jet()->getAttribute(xAOD::JetAttribute::SumPtTrkPt500,sumPtTrk); // FIXME or SumPtTrkPt1000 ??
      chf = (pt > 0.) ? sumPtTrk[0]/pt : 0;
      thisjet.jet()->btagging()->MVx_discriminant("MV2c20", jetbtag);
    }
    int jetflav = 0;
    float jettagU = -999.;
    float jettagB = -999.;
    float jettagC = -999.;
    float jetFracSamplingMax = 0.;
    int jetFracSamplingMaxIndex = 0.;
    if ( thisjet.jet() ) {
      thisjet.jet()->getAttribute(xAOD::JetAttribute::FracSamplingMax,jetFracSamplingMax);
      thisjet.jet()->getAttribute(xAOD::JetAttribute::FracSamplingMaxIndex,jetFracSamplingMaxIndex);
    }

    if (pt > 40000.) {
      ntv.jetPt.push_back(pt * 0.001);
      ntv.jetEta.push_back(eta);
      ntv.jetPhi.push_back(phi);
      ntv.jetM.push_back(m * 0.001);
      ntv.jetBTag.push_back(jetbtag);
      ntv.jetFlav.push_back(jetflav);
      ntv.jetTagU.push_back(jettagU);
      ntv.jetTagB.push_back(jettagB);
      ntv.jetTagC.push_back(jettagC);
      ntv.jetFracSamplingMax.push_back(jetFracSamplingMax);
      ntv.jetFracSamplingMaxIndex.push_back(jetFracSamplingMaxIndex);
      ntv.Njet++;
    }
    switch (jet0)
    {
    case 0:
      jet1TLV.SetPtEtaPhiM(pt,eta,phi,m);
      ntv.emfjet0 = emf;
      ntv.chfjet0 = chf;
      break;
    case 1:
      jet2TLV.SetPtEtaPhiM(pt,eta,phi,m);
      ntv.emfjet1 = emf;
      ntv.chfjet1 = chf;
      break;
    }
  }
  ntv.jetSmearSystW = jetSmearSystW;

  ntv.tauN = signal_taus.size();
  ntv.tauLooseN = baseline_taus.size();
  //std::cout<<"tauN="<<ntv.tauN<<" tauLooseN="<<ntv.tauLooseN<<std::endl;
  for ( size_t tau0=0; tau0<baseline_taus.size(); tau0++)
  {
    const TauProxy& thistau = baseline_taus[tau0];
    ntv.tauPt.push_back(thistau.Pt() * 0.001);
    ntv.tauEta.push_back(thistau.Eta());
    ntv.tauPhi.push_back(thistau.Phi());
    float sf, sfStatUp, sfStatDown, sfSystUp, sfSystDown;
    thistau.getSF(sf,sfStatUp,sfStatDown,sfSystUp,sfSystDown);
    ntv.tauLooseSF.push_back(sf);
    ntv.tauLooseSFStatUp.push_back(sfStatUp);
    ntv.tauLooseSFStatDown.push_back(sfStatDown);
    ntv.tauLooseSFSystUp.push_back(sfSystUp);
    ntv.tauLooseSFSystDown.push_back(sfSystDown);
    //std::cout<<" tauLooseSF = "<<sf<<"+"<<sfStatUp<<"-"<<sfStatDown<<"+"<<sfSystUp<<"-"<<sfSystDown<<std::endl;
  }
}

void PhysObjProxyUtils::FillTriggerBits(NTVars& ntv,
					long trigger)
{
  ntv.triggerBits.push_back(trigger);
}



void PhysObjProxyUtils::FillNTReclusteringVars(NTReclusteringVars& RTntv,
					       const std::vector<JetProxy>& good_jets,
					       std::vector<float> vReclJetMass, std::vector<float> vReclJetPt,
                                               std::vector<float> vReclJetEta, std::vector<float> vReclJetPhi,
                                               std::vector<float> vD2,std::vector<bool> visWmedium,
                                               std::vector<bool> visWtight, std::vector<bool> visZmedium,
                                               std::vector<bool> visZtight)
{
  RTntv.Reset();

  //Reclustering:
  PhysObjProxyUtils::ReclJets myRT;
  const float fCut=0.1;  //trimming cut
  myRT=Recluster(good_jets, 40000., fCut, 1.0);
  RTntv.RTjets10SubJetIndeces = myRT.recl_jets_subInds;

  std::vector<float> pts =  myRT.recl_jets_Pt;
  for ( size_t i = 0; i < pts.size(); i++ ) {pts[i] *= 0.001;}
  std::vector<float> masses =  myRT.recl_jets_M;
  for ( size_t i = 0; i < masses.size(); i++ ) {masses[i] *= 0.001;}

  RTntv.RTjetM = masses;
  RTntv.RTjetPt = pts;
  RTntv.RTjetEta = myRT.recl_jets_Eta;
  RTntv.RTjetPhi = myRT.recl_jets_Phi;

  int NWcandidates= 0;
  for (unsigned int iRT=0; iRT< RTntv.RTjetM.size(); ++iRT){
    if(RTntv.RTjetM[iRT]>60000. && RTntv.RTjetM[iRT]<100000.) NWcandidates++;
  }
  RTntv.NWcandidates= NWcandidates;

  // NEW RECLUSTERING

  RTntv.nJetsRecl  = vReclJetMass.size();
  RTntv.ReclJetMass = vReclJetMass;
  RTntv.ReclJetPt = vReclJetPt;
  RTntv.ReclJetPhi = vReclJetPhi;
  RTntv.ReclJetEta = vReclJetEta;
  RTntv.D2 = vD2;
  RTntv.isWmedium = visWmedium;
  RTntv.isWtight = visWtight;
  RTntv.isZmedium = visZmedium;
  RTntv.isZtight = visZtight;
}

PhysObjProxyUtils::ReclJets PhysObjProxyUtils::Recluster(const std::vector<JetProxy>& small_jets, double PTcut, double fcut, double jetRad){
  vector<int> NumOfSubJets;
  vector<fastjet::PseudoJet> particles = ObjsToPJ(small_jets);
  fastjet::JetDefinition fJetDef(fastjet::antikt_algorithm, jetRad,fastjet::E_scheme, fastjet::Best);
  fastjet::ClusterSequence hardClustSeq(particles, fJetDef);
  vector<fastjet::PseudoJet> StandardJets = fastjet::sorted_by_pt(hardClustSeq.inclusive_jets(PTcut));
  vector<TLorentzVector> RCjets;
  for (unsigned int i=0; i<StandardJets.size(); ++i){
    TLorentzVector sub = TLorentzVector();
    sub.SetPtEtaPhiE(StandardJets[i].pt(), StandardJets[i].eta(), StandardJets[i].phi(), StandardJets[i].e());
    vector<fastjet::PseudoJet> constituents = StandardJets[i].constituents();
    for(unsigned int iCons = 0; iCons < constituents.size(); ++iCons){
      //Do something with the small radius jets.
    }
    RCjets.push_back(sub);
  }

  vector<TLorentzVector> RTjets;
  RTjets.clear();
  vector<vector<int> > RTjets_small_jets_inds;
  vector<float> RTjets_M;

  //Now for my trimming on the re-clustered jets
  for (unsigned int i=0; i<StandardJets.size(); ++i){
    int NumSubJets=0;
    TLorentzVector trimmedjet = TLorentzVector();
    vector<fastjet::PseudoJet> constituents = StandardJets[i].constituents();
    //float emfrac_recalculated=0;
    vector<int> this_jet_subinds;
    for(unsigned int iCons = 0; iCons < constituents.size(); ++iCons){
      TLorentzVector subjet = TLorentzVector();
      subjet.SetPtEtaPhiE(constituents[iCons].pt(), constituents[iCons].eta(), constituents[iCons].phi(), constituents[iCons].e());
      //float emfrac= FindEmFrac(small_jets, subjet.Pt());
      if (subjet.Pt() > fcut*RCjets[i].Pt()){
        //emfrac_recalculated+= subjet.E()*emfrac;
	trimmedjet+=subjet;
        NumSubJets+=1;
        for (unsigned int iSub = 0; iSub < small_jets.size(); ++iSub){
            if (fabs(small_jets[iSub].Pt()-subjet.Pt())>0.00001) continue;
            if (fabs(small_jets[iSub].Phi()-subjet.Phi())>0.01) continue;
            if (fabs(small_jets[iSub].Eta()-subjet.Eta())>0.01) continue;
            this_jet_subinds.push_back(iSub);
        }
      }
    }
    //emfrac_recalculated= emfrac_recalculated/trimmedjet.E();
    //if (emfrac_recalculated<0.99) {
      NumOfSubJets.push_back(NumSubJets);
      RTjets.push_back(trimmedjet);
      RTjets_small_jets_inds.push_back(this_jet_subinds);
    //}
  }

  //finally sort RTjets according to pt
  vector<int> sorted_indexes= GetSortedJetIndexes(RTjets);
  PhysObjProxyUtils::ReclJets myRecl;
  for (unsigned int i=0; i<sorted_indexes.size(); ++i){
    myRecl.recl_jets_tlv.push_back(RTjets[sorted_indexes[i]]);
    myRecl.recl_jets_subInds.push_back(RTjets_small_jets_inds[sorted_indexes[i]]);
    myRecl.recl_jets_Pt.push_back(RTjets[sorted_indexes[i]].Pt());
    myRecl.recl_jets_Eta.push_back(RTjets[sorted_indexes[i]].Eta());
    myRecl.recl_jets_Phi.push_back(RTjets[sorted_indexes[i]].Phi());
    myRecl.recl_jets_M.push_back(RTjets[sorted_indexes[i]].M());
    myRecl.sub_jets.push_back(NumOfSubJets[sorted_indexes[i]]);
    if(i==0){
      //fill leading jet constituents
      vector<fastjet::PseudoJet> constituents = StandardJets[sorted_indexes[i]].constituents();
      for(unsigned int iCons = 0; iCons < constituents.size(); ++iCons){
        TLorentzVector subjet = TLorentzVector();
        subjet.SetPtEtaPhiE(constituents[iCons].pt(), constituents[iCons].eta(), constituents[iCons].phi(), constituents[iCons].e());
	if (subjet.Pt() > fcut*StandardJets[sorted_indexes[i]].pt()){
	  myRecl.recljet1_smalljets_tlv.push_back(subjet);
	}
      }
      //extract two latest protojets
      fastjet::PseudoJet parent1, parent2;
      if (!hardClustSeq.has_parents(StandardJets[sorted_indexes[i]],parent1,parent2)) {
          myRecl.recljet1_massdrop= 0;
	 continue;
      }
      vector<fastjet::PseudoJet> constituents1 = parent1.constituents();
      TLorentzVector parentTrim1(0,0,0,0),parentTrim2(0,0,0,0);
      for(unsigned int iCons = 0; iCons < constituents1.size(); ++iCons){
        TLorentzVector subjet = TLorentzVector();
        subjet.SetPtEtaPhiE(constituents1[iCons].pt(), constituents1[iCons].eta(),
	constituents1[iCons].phi(), constituents1[iCons].e());
	if (subjet.Pt() > fcut*StandardJets[sorted_indexes[i]].pt()){
	  parentTrim1+=subjet;
	}
      }
      vector<fastjet::PseudoJet> constituents2 = parent2.constituents();
      for(unsigned int iCons = 0; iCons < constituents2.size(); ++iCons){
        TLorentzVector subjet = TLorentzVector();
        subjet.SetPtEtaPhiE(constituents2[iCons].pt(), constituents2[iCons].eta(),
	constituents2[iCons].phi(), constituents2[iCons].e());
	if (subjet.Pt() > fcut*StandardJets[sorted_indexes[i]].pt()){
	  parentTrim2+=subjet;
	}
      }
      myRecl.recljet1_massdrop= TMath::Max(parentTrim1.M(),parentTrim2.M())/RTjets[sorted_indexes[i]].M();
    }
  }

  return myRecl;
}

std::vector<int> PhysObjProxyUtils::GetSortedJetIndexes(const std::vector<TLorentzVector> jets)
{
   float vec_pt[500];

   for (UInt_t iJet=0; iJet< jets.size(); ++iJet){
       vec_pt[iJet]= jets[iJet].Pt();
   }

   // now obtain list of sorted indexes
   int sorted_index[500];
   TMath::Sort((int)jets.size(),vec_pt,sorted_index);

   vector<int> sorted;
   if(!sorted.empty()) sorted.clear();
   for (UInt_t iJet=0; iJet< jets.size(); ++iJet){
       sorted.push_back(sorted_index[iJet]);
   }
   return sorted;

}


std::vector<fastjet::PseudoJet> PhysObjProxyUtils::ObjsToPJ(const std::vector<JetProxy>& good_jets){
  vector<fastjet::PseudoJet> out;
  for(unsigned int i = 0; i < good_jets.size(); i++){
      const JetProxy& thisjet = good_jets[i];
      fastjet::PseudoJet newjet(thisjet.Px(), thisjet.Py(), thisjet.Pz(), thisjet.E());
      out.push_back(newjet);
  }
  return out;
}
