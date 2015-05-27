
#include "ZeroLeptonRun2/ZeroLeptonCRZ.h"
#include "ZeroLeptonRun2/PhysObjProxies.h"
#include "ZeroLeptonRun2/PtOrder.h"


#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TActiveStore.h"
#include "xAODRootAccess/TStore.h"
#include "xAODTracking/Vertex.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODEgamma/Electron.h"
#include "xAODMuon/Muon.h"
#include "PATInterfaces/SystematicSet.h"
#include "cafe/Processor.h"
#include "cafe/Controller.h"
#include "cafe/Config.h"
#include "TDirectory.h"
#include "TVector2.h"
#include "TLorentzVector.h"

#include <iostream>
#include <stdexcept>
#include <algorithm>

ZeroLeptonCRZ::ZeroLeptonCRZ(const char *name)
  : cafe::Processor(name), 
    m_tree(0), 
    m_stringRegion("CRZ_SRAll"), 
    m_doSmallNtuple(true),
    m_IsData(false),
    m_IsTruth(false),
    m_IsSignal(false),
    m_UseSystematics(false),
    m_period(INVALID),
    m_isMuonChannel(false),
    m_isElectronChannel(false),
    m_suffix(""),
    m_physobjsFiller(0),
    m_physobjsFillerTruth(0),
    m_cutVal(),
    m_proxyUtils(m_IsData),
    m_ZLUtils(m_IsData, NotADerivation),
    m_counter(0),
    m_counterRepository("",false,0),
    m_treeRepository(),
    m_derivationTag(INVALID_Derivation)
{
  cafe::Config config(name);
  m_IsData = config.get("IsData",false);
  m_IsSignal = config.get("IsSignal",false);
  m_IsTruth = config.get("IsTruth",false);
  if ( m_IsData ) {
    m_isMuonChannel = config.get("IsMuonChannel",false);
    m_isElectronChannel = config.get("IsElectronChannel",false);;
  }
  else {
    m_isMuonChannel = true;
    m_isElectronChannel = true;
  }
  m_UseSystematics = config.get("UseSystematics",false);

  m_period = periodFromString(config.get("Period","p13tev"));
  if ( m_period == p7tev ) throw(std::domain_error("ZeroLeptonCRZ does not support the 7tev run period"));
  if ( m_period == INVALID ) throw(std::domain_error("ZeroLeptonCRZ: invalid run period specified"));

  m_derivationTag = derivationTagFromString(config.get("DerivationTag",""));
  if ( m_derivationTag == INVALID_Derivation ) throw(std::domain_error("ZeroLeptonSR: invalid derivation tag specified"));

  std::string cutfile = config.get("cutfile","None");
  if ( cutfile == "None" ) throw(std::domain_error("ZeroLeptonCRZ: invalid cut file specified"));
  m_cutVal.ReadCutValues(cutfile);

  m_suffix = config.get("suffix","");
  m_physobjsFiller = new PhysObjProxyFiller(20000.f,10000.f,10000.f,m_suffix);
  m_physobjsFillerTruth = new PhysObjProxyFillerTruth(20000.f,20000.f,10000.f,m_suffix);
  m_proxyUtils = PhysObjProxyUtils(m_IsData);

  m_ZLUtils = ZeroLeptonUtils(m_IsData, m_derivationTag);
}

ZeroLeptonCRZ::~ZeroLeptonCRZ()
{
  if ( !m_UseSystematics && m_counter ) delete m_counter;
  if ( m_physobjsFiller ) delete m_physobjsFiller;
  if ( m_physobjsFillerTruth ) delete m_physobjsFillerTruth;
}

TTree* ZeroLeptonCRZ::bookTree(const std::string& treename)
{
  const char* name(treename.c_str());
  TTree* tree = new TTree(name,"ZeroLepton final optimisation");
  tree->SetDirectory(getDirectory());
  bookNTVars(tree,m_ntv,false);
  bookNTExtraVars(tree,m_extrantv);
  bookNTReclusteringVars(tree,m_RTntv);
  bookNTCRZVars(tree,m_crzntv);
  return tree;
}

TTree* ZeroLeptonCRZ::getTree(const std::string& treename)
{
  std::map<std::string,TTree*>::const_iterator pos = m_treeRepository.find(treename);
  if ( pos == m_treeRepository.end() ) {
    pos = m_treeRepository.insert(std::make_pair(treename,bookTree(treename))).first;
  }
  return pos->second;
}

void ZeroLeptonCRZ::begin()
{
  std::string sSR = m_stringRegion;
  if(m_doSmallNtuple) { 
    sSR+="NT";
  }

  if ( m_UseSystematics ) {
    m_counterRepository = CounterRepository("ZeroLeptonCounter"+m_stringRegion,m_IsSignal,getDirectory());
  }
  else {
    m_counter = new Counter("ZeroLeptonCounter"+m_stringRegion,40,m_IsSignal);
    if(m_doSmallNtuple) m_tree = bookTree(sSR);
  }

}


bool ZeroLeptonCRZ::processEvent(xAOD::TEvent& event)
{
  // access the transient store
  xAOD::TStore* store = xAOD::TActiveStore::store();
  if ( m_UseSystematics ) {
    CP::SystematicSet* currentSyst = 0;
    if ( ! store->retrieve(currentSyst, "CurrentSystematicSet").isSuccess() ) throw std::runtime_error("Could not retrieve CurrentSystematicSet");
    m_counter = m_counterRepository.counter(currentSyst->name());
    if (currentSyst->name() == "" ) {
      m_tree = getTree(m_stringRegion+"NT");
    }
    else {
      m_tree = getTree(m_stringRegion+"NT_"+currentSyst->name());
    }
  }

  // eventInfo
  const xAOD::EventInfo* eventInfo = 0;
  if ( ! event.retrieve( eventInfo, "EventInfo").isSuccess() ) throw std::runtime_error("Could not retrieve EventInfo");
  uint32_t RunNumber = eventInfo->runNumber();
  unsigned long long EventNumber = eventInfo->eventNumber();
  uint32_t mc_channel_number = 0;
  if ( ! m_IsData ) mc_channel_number = eventInfo->mcChannelNumber();

  // global event weight
  float weight = 1.f;

  // get generator weight
  float genWeight = 1.f;
  if ( !m_IsData ) { 
    genWeight = eventInfo->mcEventWeight(0);
    //out() << " gen weight " << genWeight << std::endl;
    weight *= genWeight;
  }

  // get pileup weights
  std::vector<float>* pileupWeights = 0;
  if ( !m_IsData && !m_IsTruth  ) {
    if ( !store->retrieve< std::vector<float> >(pileupWeights,"pileupWeights").isSuccess() ) throw std::runtime_error("could not retrieve pileupWeights");
    //out() << " pileup weight " << (*pileupWeights)[0] << std::endl;
    //weight *= (*pileupWeights)[0];
  }
  else {
    static std::vector<float> dummy(3,0.);
    pileupWeights = &dummy;
  }

  // hardproc (see SUSYTools)
  // FIXME : to be implemented
  int trueTopo = 0;
  if ( m_IsSignal ) {
    unsigned int* finalstate = 0;
    if ( !store->retrieve<unsigned int>(finalstate,"HardProcess").isSuccess() ) throw std::runtime_error("could not retrieve HardProcess");
    trueTopo = *finalstate;
  }

  // counters
  int incr=0;
  m_counter->increment(1.,incr++,"NbOfEvents",trueTopo);
  m_counter->increment(weight,incr++,"runNumber",trueTopo);
  // FIXME do something with bin 38 & 39 ?!


  // Normalisation weight, e.g. MC cross-section vs luminosity
  // FIXME : to be implemented ?
  std::vector<float> normWeight(3,0.);

  unsigned int veto = 0;
  // MC event veto (e.g. to remove sample phase space overlap)
  if ( ! m_IsData && (m_period == p8tev || m_period == p13tev) && !m_IsTruth  ) {
    unsigned int* pveto = 0;
    if ( !store->retrieve<unsigned int>(pveto,"mcVetoCode").isSuccess() ) throw std::runtime_error("could not retrieve mcVetoCode");
    veto = *pveto;
    bool* mcaccept = 0;
    if ( !store->retrieve<bool>(mcaccept,"mcAccept").isSuccess() ) throw std::runtime_error("could not retrieve mcaccept");
    if ( ! *mcaccept ) return true;
  }
  m_counter->increment(weight,incr++,"Truth filter",trueTopo);

  // Good run list
  if ( m_IsData ) {
    bool* passGRL = 0;
    if ( !store->retrieve<bool>(passGRL,"passGRL").isSuccess() ) throw std::runtime_error("could not retrieve passGRL");
    if ( ! *passGRL ) return true;
  }
  m_counter->increment(weight,incr++,"GRL",trueTopo);

  // HFor veto
  // FIXME : to be implemented ... not in xAOD yet
  m_counter->increment(weight,incr++,"hfor veto",trueTopo);

  // Trigger selection
  // FIXME : no trigger information in xAOD yet 
  m_counter->increment(weight,incr++,"Trigger",trueTopo);

  // These jets have overlap removed
  std::vector<JetProxy> good_jets, bad_jets, b_jets, c_jets;
  if(! m_IsTruth){
    m_physobjsFiller->FillJetProxies(good_jets,bad_jets,b_jets);
  }
  if(m_IsTruth){
    m_physobjsFillerTruth->FillJetProxies(good_jets,b_jets);
  }
  std::vector<float> btag_weight(7,1.); // not implemented in SUSYTools
  std::vector<float> ctag_weight(7,1.); // not implemented in SUSYTools

  // isolated_xxx have overlap removed
  std::vector<ElectronProxy> baseline_electrons, isolated_baseline_electrons, isolated_signal_electrons;
  if(! m_IsTruth){
    m_physobjsFiller->FillElectronProxies(baseline_electrons, isolated_baseline_electrons, isolated_signal_electrons);
  }
  std::vector<ElectronTruthProxy> baseline_electrons_truth, isolated_baseline_electrons_truth, isolated_signal_electrons_truth;
  if(m_IsTruth){
    m_physobjsFillerTruth->FillElectronProxies(baseline_electrons_truth, isolated_baseline_electrons_truth, isolated_signal_electrons_truth);
  }
  // keep only signal electrons with Pt>25GeV
  for ( std::vector<ElectronProxy>::iterator it = isolated_signal_electrons.begin();
	it != isolated_signal_electrons.end(); ) {
    if ( it->Pt() < 25000. ) it = isolated_signal_electrons.erase(it);
    else it++;
  }
  // FIXME : trigger matching
  std::vector<ElectronProxy> trigmatched_electrons = isolated_signal_electrons;

  // isolated_xxx have overlap removed
  std::vector<MuonProxy> baseline_muons, isolated_baseline_muons, isolated_signal_muons;
  if(! m_IsTruth){
    m_physobjsFiller->FillMuonProxies(baseline_muons, isolated_baseline_muons, isolated_signal_muons);
  }
  std::vector<MuonTruthProxy> baseline_muons_truth, isolated_baseline_muons_truth, isolated_signal_muons_truth;
  if(m_IsTruth){
    m_physobjsFillerTruth->FillMuonProxies(baseline_muons_truth, isolated_baseline_muons_truth, isolated_signal_muons_truth);
  }
  // keep only signal muons with Pt>25GeV
  for ( std::vector<MuonProxy>::iterator it = isolated_signal_muons.begin();
	it != isolated_signal_muons.end(); ) {
    if ( it->Pt() < 25000. ) it = isolated_signal_muons.erase(it);
    else it++;
  }
  // FIXME : trigger matching
  std::vector<MuonProxy> trigmatched_muons = isolated_signal_muons;

  // missing ET
  TVector2* missingET = 0;
  if(! m_IsTruth){
    if ( ! store->retrieve<TVector2>(missingET,"SUSYMET"+m_suffix).isSuccess() ) throw std::runtime_error("could not retrieve SUSYMET"+m_suffix);
  }
  if(m_IsTruth){
    if ( ! store->retrieve<TVector2>(missingET,"TruthMET"+m_suffix).isSuccess() ) throw std::runtime_error("could not retrieve TruthMET"+m_suffix);
  }
  double MissingEt = missingET->Mod();


  // bad jet veto
  if ( !bad_jets.empty() ) return true;
  m_counter->increment(weight,incr++,"JetCleaning",trueTopo);

  // LAr, Tile, reco problems in data
  if ( m_IsData ) {
    bool* badDetectorQuality = 0 ; 
    if ( !store->retrieve<bool>(badDetectorQuality,"badDetectorQuality").isSuccess() ) throw std::runtime_error("could not retrieve badDetectorQuality");
    if ( *badDetectorQuality ) return true;
  }
  m_counter->increment(weight,incr++,"Detector cleaning",trueTopo);

  // MET track
  double MET_Track = -999.;
  double MET_Track_phi = -999.;
  if(! m_IsTruth){
    m_ZLUtils.trackMET(event, MET_Track, MET_Track_phi);
  }
  
  // primary vertex cut
  const xAOD::Vertex* primVertex = 0;
  if(! m_IsTruth){
    primVertex = ZeroLeptonUtils::GetPrimVtx(event);
    if ( !primVertex ||  !( primVertex->nTrackParticles() > 4) ) return true;
  }
  m_counter->increment(weight,incr++,"Vertex Cut",trueTopo);

  //FIXME isBadMuon not yet implement in SUSYObjDef_xAOD

  // Cosmic muon cut
  if ( m_proxyUtils.CosmicMuon(isolated_baseline_muons) ) return true;
  m_counter->increment(weight,incr++,"CosmicMuons",trueTopo);

  // FIXME fake lepton bkg estimate
  // do we need to re-implement that ?

  std::vector<TLorentzVector> leptonTLVs;
  std::vector<int> leptonCharges;
  if ( m_isMuonChannel && 
       isolated_baseline_muons.size()==2 && 
       isolated_baseline_electrons.size()==0 ) {
    leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_baseline_muons[0]))));
    leptonCharges.push_back((int)(isolated_baseline_muons[0].muon()->charge()));
    leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_baseline_muons[1]))));
    leptonCharges.push_back((int)(isolated_baseline_muons[1].muon()->charge()));
  } 
  else if (  m_isElectronChannel && 
             isolated_baseline_electrons.size()==2 && 
             isolated_baseline_muons.size()==0) {
    leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_baseline_electrons[0]))));
    leptonCharges.push_back((int)(isolated_baseline_electrons[0].electron()->charge()));
    leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_baseline_electrons[1]))));
    leptonCharges.push_back((int)(isolated_baseline_electrons[1].electron()->charge()));
  }
  if ( leptonTLVs.empty() ) return true;
  TLorentzVector dileptonTLV = leptonTLVs[0]+leptonTLVs[1];
  if ( leptonCharges[0]*leptonCharges[1] > 0 ) return true;
  m_counter->increment(weight,incr++,"2 OS Baseline Leptons",trueTopo);

  // leading lepton is signal lepton and trigger matched
  /*
  if ( m_isMuonChannel ) {
    if ( isolated_signal_muons.empty() ) return true;
    if ( trigmatched_muons.empty() ) return true;
    if ( isolated_signal_muons[0].muon() != trigmatched_muons[0].muon() ) return true;
  }
  if ( m_isElectronChannel ) {
    if ( isolated_signal_electrons.empty() ) return true;
    if ( trigmatched_electrons.empty() ) return true;
    if ( isolated_signal_electrons[0].electron() != trigmatched_electrons[0].electron() ) return true;
  }
  // both leading leptons are signal leptons
  if ( m_isMuonChannel && isolated_signal_muons.size()!=2 ) return true;
  if ( m_isElectronChannel && isolated_signal_electrons.size()!=2 ) return true;
  m_counter->increment(weight,incr++,"2 OS Signal Leptons",trueTopo);
  */

  double InvMassLepPair = dileptonTLV.M();
  if (InvMassLepPair < 66000 || InvMassLepPair > 116000) return true;
  m_counter->increment(weight,incr++,"M_ll Cut",trueTopo);


  // FIXME: Apply Lepton scale factors
  
  // Add lepton to jets (SR) or MET (VR)
  TVector2 missingETPrime =  *missingET;
  missingETPrime = missingETPrime +
    TVector2(leptonTLVs[0].Px(),leptonTLVs[0].Py())  + 
    TVector2(leptonTLVs[1].Px(),leptonTLVs[1].Py());
  double MissingEtPrime = missingETPrime.Mod();
  

  // bad muons for MET cut: based on non isolated muons
  // FIXME do something special with isbadMETmuon when there are signal muons
  //if ( m_proxyUtils.isbadMETmuon(baseline_muons, MissingEt, *missingET) ) return true;
  //m_counter->increment(weight,incr++,"IsBadMETMuon",trueTopo);

  // bad Tile cut
  if ( m_proxyUtils.badTileVeto(good_jets,*missingET)) return true;
  m_counter->increment(weight,incr++,"Bad Tile Veto",trueTopo);


  // Negative-cell cleaning cut
  bool HasNegCell = 0 ; 
  if(! m_IsTruth){
    HasNegCell = m_ZLUtils.NegCellCleaning(event,*missingET);
  }
  //out() << " NegCell " << HasNegCell << std::endl;
  if ( HasNegCell ) return true;
  m_counter->increment(weight,incr++,"Negative-cell cleaning",trueTopo);



  // Topology cut
  // Selection
  bool inSRmono=false;
  bool inSR1=false;
  bool inSR2=false;
  bool inSR3=false;
  bool inSR4=false;
  bool inSR5=false;
  bool inSR6=false;

  if (good_jets.size()<1) return true;  
  if ((good_jets.size()>=1)) inSRmono=true; 
  if ((good_jets.size()>=2)) {inSR1=true; inSR2=true;}
  if ((good_jets.size()>=3)) inSR3=true;
  if ((good_jets.size()>=4)) inSR4=true;
  if ((good_jets.size()>=5)) inSR5=true;
  if ((good_jets.size()>=6)) inSR6=true;
  if (!(inSRmono||inSR1||inSR2||inSR3||inSR4||inSR5||inSR6)) return true;
  m_counter->increment(weight,incr++,"SignalRegion vs. Number of jets",trueTopo);
  
  // jet timing cut
  std::vector<float> time;
  m_proxyUtils.EnergyWeightedTime(good_jets,time);
  if ((good_jets.size()>=2) && (fabs(time[0])>5)) { inSR1=false; }
  if ((good_jets.size()>=2) && (fabs(time[0])>5)) { inSR2=false; }
  if ((good_jets.size()>=3) && (fabs(time[1])>5)) { inSR3=false; }
  if ((good_jets.size()>=4) && (fabs(time[2])>5)) { inSR4=false; }
  if ((good_jets.size()>=5) && (fabs(time[3])>5)) { inSR5=false; }
  if ((good_jets.size()>=6) && (fabs(time[4])>5)) { inSR6=false; }
  if (!(inSRmono||inSR1||inSR2||inSR3||inSR4||inSR5||inSR6)) return true;
  m_counter->increment(weight,incr++,"Energy weighted time ",trueTopo);

  // MissingET cut
  if (!(MissingEtPrime > m_cutVal.m_cutEtMiss)) return true;
  m_counter->increment(weight,incr++,"MET cut",trueTopo);

  // Leading jet Pt cut
  if (!(good_jets[0].Pt() > m_cutVal.m_cutJetPt0)) return true; 
  m_counter->increment(weight,incr++,"1 jet Pt > 130 GeV Selection",trueTopo);

  // Other jet pT cuts depending on selection
  inSRmono = (MissingEt > m_cutVal.m_cutEtMiss1Jet);
  inSR1 = false;
  inSR2 = false;
  inSR3 = false;
  inSR4 = false;
  inSR5 = false;
  inSR6 = false;
  
  if ((good_jets.size() > 1) && 
      (good_jets[1].Pt() > m_cutVal.m_cutJetPt1)) {
    inSR1 = true; inSR2 = true;
  }
  if (inSR2 && (good_jets.size() > 2) &&
      (good_jets[2].Pt() > m_cutVal.m_cutJetPt2)) inSR3 = true;
  if (inSR3 && (good_jets.size() > 3) &&
      (good_jets[3].Pt() > m_cutVal.m_cutJetPt3)) inSR4 = true;
  if (inSR4 && (good_jets.size() > 4) &&
      (good_jets[4].Pt() > m_cutVal.m_cutJetPt4)) inSR5 = true;
  if (inSR5 && (good_jets.size() > 5) &&
      (good_jets[5].Pt() > m_cutVal.m_cutJetPt5)) inSR6 = true;
  if (!(inSRmono||inSR1||inSR2||inSR3||inSR4||inSR5||inSR6)) return true;
  m_counter->increment(weight,incr++,"jet Pt Selection",trueTopo);


  // Calculate variables for ntuple -----------------------------------------
  double phi_met = TMath::ATan2(missingETPrime.Y(),missingETPrime.X());
  double minDphi = m_proxyUtils.SmallestdPhi(good_jets,phi_met);
  double RemainingminDPhi = m_proxyUtils.SmallestRemainingdPhi(good_jets,phi_met);
  //out() << " minDphi " << minDphi << " " << RemainingminDPhi << std::endl;

  //out() << " MissingEt " << MissingEt << std::endl;
  //out() << " jet Pt ";
  //for ( size_t i = 0; i<good_jets.size(); ++i ) out() << " " << good_jets[i].Pt();
  //out() << std::endl;
  double Meff[6];
  for ( size_t i = 0; i<6; ++i ) {
    Meff[i] = m_proxyUtils.Meff(good_jets,
				std::max<size_t>(1,i+1),
				MissingEtPrime,
				m_cutVal.m_cutJetPt1,
				m_cutVal.m_cutJetPt4);
    //out() << " Meff[" << i <<"] " << Meff[i] << std::endl;
  }
  double meffincl = m_proxyUtils.Meff(good_jets,
				      good_jets.size(),
				      MissingEtPrime,
				      m_cutVal.m_cutJetPt4,
				      m_cutVal.m_cutJetPt4);
  //out() << " MeffInc " << meffincl << std::endl;

  // ttbar reweighting not available yet in SUSYOBJDef_xAOD
  float ttbarWeightHT = 1.;
  float ttbarWeightPt2 = 1.;
  float ttbarAvgPt = 0.;
  float HT = meffincl -  MissingEt;

  // Sherpa MassiveCB W/Z reweighting : not implemented yet in SUSYOBJDef_xAOD
  float WZweight = 1.;



  double mT2=-9; 
  double mT2_noISR=-9; 
  //if (good_jets.size()>=2) mT2 = m_proxyUtils.MT2(good_jets,missingETPrime);

  //Super Razor variables
  double gaminvRp1 =-999;
  double shatR =-999;
  double mdeltaR =-999;
  double cosptR =-999;
  double Minv2 =-999;
  double Einv =-999;
  double  gamma_R=-999;
  double dphi_BETA_R =-999; 
  double dphi_leg1_leg2 =-999; 
  double costhetaR =-999;
  double dphi_BETA_Rp1_BETA_R=-999;
  double gamma_Rp1=-999;
  double Eleg1=-999;
  double Eleg2=-999; 
  double costhetaRp1=-999;

  m_proxyUtils.RazorVariables(good_jets, 
			      missingETPrime.X(),
			      missingETPrime.Y(),
			      gaminvRp1 ,
			      shatR ,
			      mdeltaR ,
			      cosptR ,
			      Minv2 ,
			      Einv ,
			      gamma_R,
			      dphi_BETA_R , 
			      dphi_leg1_leg2 , 
			      costhetaR ,
			      dphi_BETA_Rp1_BETA_R,
			      gamma_Rp1,
			      Eleg1,
			      Eleg2, 
			      costhetaRp1);

  double Sp,ST,Ap=-1;
  m_proxyUtils.ComputeSphericity(good_jets, Sp,ST,Ap);


  if(m_doSmallNtuple) { 
    unsigned int runnum = RunNumber;
    if ( ! m_IsData ) runnum = mc_channel_number;

    std::vector<float> jetSmearSystW;

    // other cleaning tests
    unsigned int cleaning = 0;
    if (fabs(time[0]) > 5) cleaning += 2;  //t2j use for all SRs

    if(! m_IsTruth){
      // FIXME why not in CRWT ?
      //bool chfTileVeto =  m_proxyUtils.chfTileVeto(good_jets);
      //if ( chfTileVeto ) cleaning += 4;
      
      bool chfVeto = m_proxyUtils.chfVeto(good_jets);
      if ( chfVeto ) cleaning += 8;
    }

    m_proxyUtils.FillNTVars(m_ntv, runnum, EventNumber, veto, weight, normWeight, *pileupWeights, genWeight,ttbarWeightHT,ttbarWeightPt2,ttbarAvgPt,WZweight, btag_weight, ctag_weight, b_jets.size(), c_jets.size(), MissingEtPrime, phi_met, Meff, meffincl, minDphi, RemainingminDPhi, good_jets, trueTopo, cleaning, time[0],jetSmearSystW,0, 0., 0.,m_IsTruth);

    m_proxyUtils.FillNTExtraVars(m_extrantv, MET_Track, MET_Track_phi, mT2, mT2_noISR, gaminvRp1, shatR, mdeltaR, cosptR, gamma_R,dphi_BETA_R, dphi_leg1_leg2, costhetaR, dphi_BETA_Rp1_BETA_R, gamma_Rp1, costhetaRp1, Ap);
      

    if( !m_IsTruth ){
      m_proxyUtils.FillNTReclusteringVars(m_RTntv,good_jets);
    }
    
    FillCRZVars(m_crzntv, leptonTLVs, *missingET, leptonCharges);

    m_tree->Fill();
  }
  return true;
}

void ZeroLeptonCRZ::finish()
{
  if ( m_UseSystematics ) {
    out() << m_counterRepository << std::endl;
  } 
  else {
    out() << *m_counter << std::endl;
  }
}

void ZeroLeptonCRZ::FillCRZVars(NTCRZVars& crzvars, std::vector<TLorentzVector>& leptons, const TVector2& metv, std::vector<int> lepsigns)
{
  crzvars.Reset();
  crzvars.lep1sign = lepsigns.at(0);
  crzvars.lep1Pt  = (leptons.at(0)).Pt();
  crzvars.lep1Eta = (leptons.at(0)).Eta();
  crzvars.lep1Phi = (leptons.at(0)).Phi();

  crzvars.lep2sign = lepsigns.at(1);
  crzvars.lep2Pt  = (leptons.at(1)).Pt();
  crzvars.lep2Eta = (leptons.at(1)).Eta();
  crzvars.lep2Phi = (leptons.at(1)).Phi();


  //double met = std::sqrt(metv.Px()*metv.Px()+metv.Py()*metv.Py());
  crzvars.mll = (leptons.at(0)+leptons.at(1)).M();

  double zpx = leptons.at(0).Px()+leptons.at(1).Py();
  double zpy = leptons.at(0).Py()+leptons.at(1).Py();
  crzvars.Zpt = std::sqrt(zpx*zpx+zpy*zpy);
}



ClassImp(ZeroLeptonCRZ);

