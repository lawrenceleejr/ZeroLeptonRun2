#ifndef PhysObjProxyUtils_h
#define PhysObjProxyUtils_h

class JetProxy;
class MuonProxy;
class TauProxy;
class NTVars;
class NTReclusteringVars;
class NTExtraVars;
class NTRJigsawVars;
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

#include <vector>
#include <map>
#include "fastjet/PseudoJet.hh"

#include "RestFrames/RestFrames.hh"

//-----------------------------------------------------------------------
// A collection of utility function for variables based on Proxy objects
//----------------------------------------------------------------------
class PhysObjProxyUtils
{
 public:
  PhysObjProxyUtils(bool IsData);
  ~PhysObjProxyUtils() ;

  void EnergyWeightedTime(const std::vector<JetProxy>& jets, std::vector<float>& time) const;

  bool badTileVeto(const std::vector<JetProxy>& jets, const TVector2& MissingET) const;

  bool chfVeto(const std::vector<JetProxy>& jets) const;
  bool chfTileVeto(const std::vector<JetProxy>& jets) const;

  double SmallestdPhi(const std::vector<JetProxy>& jets, double met_phi) const;

  double SmallestRemainingdPhi(const std::vector<JetProxy>& jets, double met_phi) const;

  void GetAlphaISRVar(const std::vector<JetProxy>& jets, double met,
		      std::vector<double>& alpha_vec) const;
  void GetMinPtDistinctionISR(const std::vector<JetProxy>& jets,
			      std::vector<double>& minPtDistinction_vec) const;
  void GetMinDeltaFraction(const std::vector<JetProxy>& jets,
			   std::vector<double>& minDeltaFrac_vec) const;
  void GetMinRapidityGap(const std::vector<JetProxy>& jets,
			 std::vector<double>& minRapGap_vec) const;
  void GetMaxRapidityOtherJets(const std::vector<JetProxy>& jets,
			       std::vector<double>& maxRapOtherJets_vec) const;
  void GetdPhiJetMet(const std::vector<JetProxy>& jets, double met_phi,
		     std::vector<double>& dPhiJetMet_vec) const;
  void GetISRJet(const std::vector<JetProxy>& jets,
		 std::vector<size_t>& isr_jet_indices,
		 double met,
		 double phi_met,
		 std::string signal,
		 bool usealpha) const;

  double Meff(const std::vector<JetProxy>& jets, size_t njets, double met, double jetPtCut, double extraJetPtCut);

  void ComputeSphericity(const std::vector<JetProxy>& jets, double & S, double & ST, double & Ap);

  void RazorVariables(const std::vector<JetProxy>& jets,
			Double_t metx,
			Double_t  mety,
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
			double &costhetaRp1);

  void RJigsawInit();

  void CalculateRJigsawVariables(const std::vector<JetProxy>& jets,
			Double_t metx,
			Double_t mety,
			std::map<TString,float>& RJigsawVariables,
			Double_t jetPtCut=0.);

  //double MT2(const std::vector<JetProxy>& jets,const TVector2& MissingET) const;

  bool CosmicMuon(const std::vector<MuonProxy>& muons) const;


  bool isbadMETmuon(const std::vector<MuonProxy>& muons,
		    float MET, const TVector2& MissingET) const;


  void FillNTVars(NTVars& ntv,
		  unsigned int RunNumber, unsigned int EventNumber,
		  unsigned int LumiBlockNumber,
		  unsigned int veto, float weight,
		  std::vector<float>& normWeight,
		  std::vector<float>& pileupWeight, float genWeight,
		  float ttbarWeightHT, float ttbarWeightPt2,
		  float ttbarAvgPt, float WZweight,
		  std::vector<float>& bTagWeight,
		  std::vector<float>& cTagWeight,int nBJet, int nCJet,
		  double MissingEt, double METPhi, double* Meff,
		  double meffincl, double minDphi, double RemainingminDPhi,
		  const std::vector<JetProxy>& good_jets,
		  int hardproc, unsigned int cleaning, float timing,
		  const std::vector<float>& jetSmearSystW,
		  const std::vector<float>* flaggedtau = NULL,
		  float tauMt = 0.f, float SherpaBugMET = 0.f,
		  bool istruth = false,
		  std::vector<TauProxy> baseline_taus = std::vector<TauProxy>(),
		  std::vector<TauProxy> signal_taus   = std::vector<TauProxy>());

  void FillNTReclusteringVars(NTReclusteringVars& RTntv,
		  const std::vector<JetProxy>& good_jets);

  void FillNTExtraVars(NTExtraVars& extrantv,
		       double MET_Track,
		       double MET_Track_phi,
		       double MT2,
		       double MT2_noISR,
		       double Ap);


  void FillNTRJigsawVars(NTRJigsawVars& rjigsawntv, std::map<TString,float> & RJigsawVariables );

  class ReclJets {
      public:
        std::vector<TLorentzVector> recl_jets_tlv;
        std::vector<int> sub_jets;
	std::vector<TLorentzVector> recljet1_smalljets_tlv;
        std::vector<std::vector<int>> recl_jets_subInds;
        std::vector<float> recl_jets_Pt;
        std::vector<float> recl_jets_Eta;
        std::vector<float> recl_jets_Phi;
        std::vector<float> recl_jets_M;

	double recljet1_massdrop;

	ReclJets(){  //constructor
	  recl_jets_tlv.clear();
	  sub_jets.clear();
          recl_jets_subInds.clear();
	  recljet1_smalljets_tlv.clear();
	  recljet1_massdrop= 0.;
          recl_jets_Pt.clear();
          recl_jets_Eta.clear();
          recl_jets_Phi.clear();
          recl_jets_M.clear();

	}
   };
  PhysObjProxyUtils::ReclJets Recluster(const std::vector<JetProxy>& good_jets, double PTcut=15, double fcut=0.05, double jetRad=1.0);

  std::vector<fastjet::PseudoJet> ObjsToPJ(const std::vector<JetProxy>& good_jets);
  std::vector<int> GetSortedJetIndexes(const std::vector<TLorentzVector> jets);


 private:
  bool m_IsData;

  TRandom3 m_random;

  RestFrames::LabRecoFrame  LAB; //!
  RestFrames::DecayRecoFrame  PP; //!
  RestFrames::DecayRecoFrame  Pa; //!
  RestFrames::DecayRecoFrame  Pb; //!
  RestFrames::DecayRecoFrame  Ca; //!
  RestFrames::DecayRecoFrame  Cb; //!
  RestFrames::SelfAssemblingRecoFrame  SAV1a; //!
  RestFrames::SelfAssemblingRecoFrame  SAV1b; //!
  RestFrames::SelfAssemblingRecoFrame  SAV2a; //!
  RestFrames::SelfAssemblingRecoFrame  SAV2b; //!
  RestFrames::VisibleRecoFrame  V1a; //!
  RestFrames::VisibleRecoFrame  V1b; //!
  RestFrames::VisibleRecoFrame  V2a; //!
  RestFrames::VisibleRecoFrame  V2b; //!
  RestFrames::InvisibleRecoFrame  Ia; //!
  RestFrames::InvisibleRecoFrame  Ib; //!

  RestFrames::InvisibleGroup              INV              ; //!
  RestFrames::CombinatoricGroup           VIS              ; //!
  RestFrames::SetMassInvJigsaw            InvMassJigsaw    ; //!
  RestFrames::SetRapidityInvJigsaw        InvRapidityJigsaw; //!
  RestFrames::ContraBoostInvJigsaw        InvCBJigsaw      ; //!
  RestFrames::MinMassesCombJigsaw         CombPPJigsaw     ; //!
  RestFrames::MinMassesCombJigsaw         CombPaJigsaw     ; //!
  RestFrames::MinMassesCombJigsaw         CombPbJigsaw     ; //!
  RestFrames::LabRecoFrame                LAB_bkg          ; //!
  RestFrames::SelfAssemblingRecoFrame     S_bkg            ; //!
  RestFrames::VisibleRecoFrame            V_bkg            ; //!
  RestFrames::InvisibleRecoFrame          I_bkg            ; //!
  RestFrames::InvisibleGroup              INV_bkg          ; //!
  RestFrames::CombinatoricGroup           VIS_bkg          ; //!
  RestFrames::SetMassInvJigsaw            InvMass_bkg      ; //!
  RestFrames::SetRapidityInvJigsaw        InvRapidity_bkg  ; //!
  RestFrames::InvisibleGroup              INV              ; //!
  RestFrames::CombinatoricGroup           VIS              ; //!
  RestFrames::SetMassInvJigsaw            InvMassJigsaw    ; //!
  RestFrames::SetRapidityInvJigsaw        InvRapidityJigsaw; //!
  RestFrames::ContraBoostInvJigsaw        InvCBJigsaw      ; //!
  RestFrames::MinMassesCombJigsaw         CombPPJigsaw     ; //!
  RestFrames::MinMassesCombJigsaw         CombPaJigsaw     ; //!
  RestFrames::MinMassesCombJigsaw         CombPbJigsaw     ; //!
  RestFrames::LabRecoFrame                LAB_bkg          ; //!
  RestFrames::SelfAssemblingRecoFrame     S_bkg            ; //!
  RestFrames::VisibleRecoFrame            V_bkg            ; //!
  RestFrames::InvisibleRecoFrame          I_bkg            ; //!
  RestFrames::InvisibleGroup              INV_bkg          ; //!
  RestFrames::CombinatoricGroup           VIS_bkg          ; //!
  RestFrames::SetMassInvJigsaw            InvMass_bkg      ; //!
  RestFrames::SetRapidityInvJigsaw        InvRapidity_bkg  ; //!
};

#endif
