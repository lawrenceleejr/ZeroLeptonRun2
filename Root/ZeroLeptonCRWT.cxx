
#include "ZeroLeptonRun2/ZeroLeptonCRWT.h"
#include "ZeroLeptonRun2/PhysObjProxies.h"
#include "ZeroLeptonRun2/PtOrder.h"
#include "ZeroLeptonRun2/CleaningHelper.h"

#include "ZeroLeptonRun2/BosonTagging.h"

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
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthParticle.h"

#include <iostream>
#include <stdexcept>
#include <algorithm>

ZeroLeptonCRWT::ZeroLeptonCRWT(const char *name)
  : cafe::Processor(name),
    m_tree(0),
    m_stringRegion("CRWT_SRAll"),
    m_doSmallNtuple(true),
    m_fillTRJigsawVars(true),
    m_fillReclusteringVars(true),
    m_doRecl(false),
    m_IsData(false),
    m_IsTruth(false),
    m_IsSignal(false),
    m_DoSystematics(false),
    m_period(INVALID),
    m_isVR(false),
    // m_isMuonChannel(false),
    // m_isElectronChannel(false),
    m_LowPtLepton(false),
    m_doFake(false),
    m_suffix(""),
    m_suffixRecl(""),
    m_physobjsFiller(0),
    m_physobjsFillerTruth(0),
    m_cutVal(),
    m_proxyUtils(m_IsData),
    m_ZLUtils(m_IsData),
    m_counter(0),
    m_counterRepository("",false,0),
    m_treeRepository()
{
  cafe::Config config(name);
  m_fillTRJigsawVars = config.get("fillTRJigsawVars",true);
  m_IsData = config.get("IsData",false);
  m_IsSignal = config.get("IsSignal",false);
  m_suffixRecl = config.get("suffixRecl","");
  m_doRecl = config.get("doRecl",false);
  m_IsTruth = config.get("IsTruth",false);
  m_DoSystematics = config.get("DoSystematics",false);
  m_doFake = config.get("doFake", false);

  m_period = periodFromString(config.get("Period","p13tev"));
  if ( m_period == p7tev ) throw(std::domain_error("ZeroLeptonCRWT does not support the 7tev run period"));
  if ( m_period == INVALID ) throw(std::domain_error("ZeroLeptonCRWT: invalid run period specified"));

  // if ( m_IsData && m_period == p8tev ) {
  //   m_isMuonChannel = config.get("IsMuonChannel",false);
  //   m_isElectronChannel = config.get("IsElectronChannel",false);;
  // }
  // else {
  //   m_isMuonChannel = true;
  //   m_isElectronChannel = true;
  // }

  m_LowPtLepton = config.get("LowPtLepton",false);
  if ( m_LowPtLepton ) m_stringRegion = "CRWTLPT_SRAll";

  if( m_doFake ){
    std::string RootCoreBin = std::getenv("ROOTCOREBIN");

    m_elecRealEff = new LeptonEfficiency;
    m_elecRealEff->initEfficiencyFile( RootCoreBin+"/data/ZeroLeptonRun2/electronRealEfficiency_Z_v05.root");
    m_elecRealEff->setHistName("el_pt_el_eta_HLT_e24_lhmedium_iloose_L1EM18VH_1e_signal_mll80100_signal");
    m_elecRealEff->setPtMax(300000.);
    m_elecRealEff->setPtMin(10000.);
    m_elecRealEff->setEtaMax(2.47);

    m_muonRealEff = new LeptonEfficiency;
    m_muonRealEff->initEfficiencyFile( RootCoreBin+"/data/ZeroLeptonRun2/muonRealEfficiency_Z_v05.root");
    m_muonRealEff->setHistName("mu_pt_mu_eta_HLT_mu24_iloose_L1MU15_1mu_signal_mll80100_signal");
    m_muonRealEff->setPtMax(300000.);
    m_muonRealEff->setPtMin(10000.);
    m_muonRealEff->setEtaMax(2.7);

    m_elecFakeEff = new LeptonEfficiency;
    m_elecFakeEff->initEfficiencyFile( RootCoreBin+"/data/ZeroLeptonRun2/electronFakeRate_data_v05.root");
    m_elecFakeEff->setHistName("el_pt_el_eta_HLT_e24_lhvloose_L1EM18VH_mt40_MET30_signal" );
    m_elecFakeEff->setPtMax(400000.);
    m_elecFakeEff->setPtMin(25000.);
    m_elecFakeEff->setEtaMax(2.47);

    m_muonFakeEff = new LeptonEfficiency;
    m_muonFakeEff->initEfficiencyFile( RootCoreBin+"/data/ZeroLeptonRun2/muonFakeRate_data_v05.root");
    m_muonFakeEff->setHistName("mu_pt_mu_eta_HLT_mu24_mt40_MET30_signal");
    m_muonFakeEff->setPtMax(60000.);
    m_muonFakeEff->setPtMin(25000.);
    m_muonFakeEff->setEtaMax(2.4);
  }

  m_isVR = config.get("IsVR",false);
  if ( m_isVR ) {
    m_stringRegion = "VRWT_SRAll";
    if ( m_LowPtLepton ) m_stringRegion = "VRWTLPT_SRAll";
  }

  std::string cutfile = config.get("cutfile","None");
  if ( cutfile == "None" ) throw(std::domain_error("ZeroLeptonCRWT: invalid cut file specified"));
  m_cutVal.ReadCutValues(cutfile);

  m_suffix = config.get("suffix","");

  m_suffixSyst = "";

  m_physobjsFiller = new PhysObjProxyFiller(20000.f,10000.f,10000.f,25000.f,m_suffix,m_doRecl,m_suffixRecl,m_suffixSyst);
  m_physobjsFillerTruth = new PhysObjProxyFillerTruth(20000.f,10000.f,10000.f,25000.f,m_suffix);
  m_proxyUtils = PhysObjProxyUtils(m_IsData);

  m_ZLUtils = ZeroLeptonUtils(m_IsData);
}

ZeroLeptonCRWT::~ZeroLeptonCRWT()
{
  if ( !m_DoSystematics && m_counter ) delete m_counter;
  if ( m_physobjsFiller ) delete m_physobjsFiller;
  if ( m_physobjsFillerTruth ) delete m_physobjsFillerTruth;
}

TTree* ZeroLeptonCRWT::bookTree(const std::string& treename)
{
  const char* name(treename.c_str());
  TTree* tree = new TTree(name,"ZeroLepton final optimisation");
  tree->SetDirectory(getDirectory());
  bookNTVars(tree,m_ntv,false);
  if ( m_fillReclusteringVars ) bookNTReclusteringVars(tree,m_RTntv);
  bookNTCRWTVars(tree,m_crwtntv);
  bookNTExtraVars(tree,m_extrantv);
  if ( m_fillTRJigsawVars) bookNTRJigsawVars(tree,m_rjigsawntv);
  return tree;
}

TTree* ZeroLeptonCRWT::getTree(const std::string& treename)
{
  std::map<std::string,TTree*>::const_iterator pos = m_treeRepository.find(treename);
  if ( pos == m_treeRepository.end() ) {
    pos = m_treeRepository.insert(std::make_pair(treename,bookTree(treename))).first;
  }
  return pos->second;
}

void ZeroLeptonCRWT::begin()
{
  std::string sSR = m_stringRegion;
  if(m_doSmallNtuple) {
    sSR+="NT";
  }

  if ( m_DoSystematics ) {
    m_counterRepository = CounterRepository("ZeroLeptonCounter"+m_stringRegion,m_IsSignal,getDirectory());
  }
  else {
    m_counter = new Counter("ZeroLeptonCounter"+m_stringRegion,40,m_IsSignal);
    if(m_doSmallNtuple) m_tree = bookTree(sSR);
  }

  if (  m_fillTRJigsawVars ) {    m_proxyUtils.RJigsawInit(); }

}


bool ZeroLeptonCRWT::processEvent(xAOD::TEvent& event)
{
  // access the transient store
  xAOD::TStore* store = xAOD::TActiveStore::store();
  std::string systag = "";
  if ( m_DoSystematics ) {
    CP::SystematicSet* currentSyst = 0;
    if ( ! store->retrieve(currentSyst, "CurrentSystematicSet").isSuccess() ) throw std::runtime_error("Could not retrieve CurrentSystematicSet");
    m_counter = m_counterRepository.counter(currentSyst->name());
    std::string sysname = currentSyst->name();
    if (sysname != "" ) systag = "_"+sysname+"_";
    if ( sysname == "" ) {
      m_tree = getTree(m_stringRegion+"NT");
      m_physobjsFiller->setSuffix(m_suffix);
      m_physobjsFiller->setSuffixSyst(m_suffixSyst);
    }
    else {
      m_tree = getTree(m_stringRegion+"NT_"+sysname);
      m_physobjsFiller->setSuffix(m_suffix+systag);
      m_physobjsFiller->setSuffixSyst(m_suffixSyst+systag);
    }
  }

  // eventInfo
  const xAOD::EventInfo* eventInfo = 0;
  if ( ! event.retrieve( eventInfo, "EventInfo").isSuccess() ) throw std::runtime_error("Could not retrieve EventInfo");
  uint32_t RunNumber = eventInfo->runNumber();
  unsigned long long EventNumber = eventInfo->eventNumber();
  uint32_t LumiBlockNumber = eventInfo->lumiBlock();
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
  unsigned long long PRWHash = 0;
  if ( !m_IsData && !m_IsTruth ) {
    if ( !store->retrieve< std::vector<float> >(pileupWeights,"pileupWeights").isSuccess() ) throw std::runtime_error("could not retrieve pileupWeights");
    //out() << " pileup weight " << (*pileupWeights)[0] << std::endl;
    //weight *= (*pileupWeights)[0];
    unsigned long long *storedPRWHash = 0;
    if ( !store->retrieve< unsigned long long >(storedPRWHash, "PRWHash").isSuccess() ) throw std::runtime_error("could not retrieve PRWHash");
    PRWHash = *storedPRWHash;
  }
  else {
    static std::vector<float> dummy(3,1.);
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
  if ( ! m_IsData && (m_period == p8tev || m_period == p13tev2015 || m_period == p13tev2016) && !m_IsTruth ) {
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

  // primary vertex cut
  const xAOD::Vertex* primVertex = 0;
  if(! m_IsTruth){
    primVertex = ZeroLeptonUtils::GetPrimVtx(event);
    if ( !primVertex ) return true;
  }
  m_counter->increment(weight,incr++,"Vertex Cut",trueTopo);

  // HFor veto
  // FIXME : to be implemented ... not in xAOD yet
  m_counter->increment(weight,incr++,"hfor veto",trueTopo);

  // Trigger selection

  bool passEltrigger=false;
  bool passMutrigger=false;
  double prescaleEl = 1.;
  double prescaleMu = 1.;
  if(! m_IsTruth){
    if ( m_LowPtLepton ) {
      if ( m_IsData && m_period == p13tev2015 ) {
        if( !(int)eventInfo->auxdata<char>("HLT_xe70")==1) return true;
      } else if ( m_IsData && m_period == p13tev2016 ){
        if( !(int)eventInfo->auxdata<char>("HLT_xe80_tc_lcw_L1XE50")==1) return true;
      } else {
        if( !(int)eventInfo->auxdata<char>("HLT_xe70")==1) return true;
        if( !(int)eventInfo->auxdata<char>("HLT_xe80_tc_lcw_L1XE50")==1) return true;
      }
    }
    else {

      // Higher threshold electron triggers
      if( (int)eventInfo->auxdata<char>("HLT_e60_lhmedium")==1  ||
          (int)eventInfo->auxdata<char>("HLT_e120_lhloose")==1  ||
          (int)eventInfo->auxdata<char>("HLT_e60_lhmedium_nod0")==1 ||
          (int)eventInfo->auxdata<char>("HLT_e140_lhloose_nod0")==1  )  passEltrigger = true;

      if ( m_IsData && m_period == p13tev2015 ) {
        if ( (int)eventInfo->auxdata<char>("HLT_e24_lhmedium_L1EM20VH")==1 ) passEltrigger = true;
      } else if ( m_IsData && m_period == p13tev2016 )  {
        if ( (int)eventInfo->auxdata<char>("HLT_e24_lhtight_nod0_ivarloose")==1 ) passEltrigger = true;
      } else {
        if ( (int)eventInfo->auxdata<char>("HLT_e24_lhmedium_L1EM18VH")==1 ) passEltrigger = true;
        if ( (int)eventInfo->auxdata<char>("HLT_e24_lhtight_nod0_ivarloose")==1 ) passEltrigger = true;
      }

      // Muon Triggers
      if((int)eventInfo->auxdata<char>("HLT_mu50")==1               ||
         (int)eventInfo->auxdata<char>("HLT_mu40")==1               )
        passMutrigger = true;
      if ( m_IsData && m_period == p13tev2015 ) {
        if ( (int)eventInfo->auxdata<char>("HLT_mu20_iloose_L1MU15")==1 ) passMutrigger = true;
      } else if ( m_IsData && m_period == p13tev2016 )  {
        if ( (int)eventInfo->auxdata<char>("HLT_mu24_ivarloose")==1 ) passMutrigger = true;
        if ( (int)eventInfo->auxdata<char>("HLT_mu24_ivarmedium")==1 ) passMutrigger = true;
      } else {
        if ( (int)eventInfo->auxdata<char>("HLT_mu20_iloose_L1MU15")==1 ) passMutrigger = true;
        if ( (int)eventInfo->auxdata<char>("HLT_mu24_ivarloose_L1MU15")==1 ) passMutrigger = true;
        if ( (int)eventInfo->auxdata<char>("HLT_mu24_ivarmedium")==1 ) passMutrigger = true;
      }

      if( !(passEltrigger || passMutrigger) ){
        if( m_doFake ){
          if( (int)eventInfo->auxdata<char>("HLT_e24_lhvloose_L1EM20VH")==1 ) prescaleEl = 35.;
          else prescaleEl = 0.;

          if( (int)eventInfo->auxdata<char>("HLT_mu20_L1MU15")==1 ) prescaleMu = 10.;
          else prescaleMu = 0.;

          if( prescaleEl<1. && prescaleMu<1. ) return true;
          else weight = 0.;
        }else{
          return true;
        }
      }
    }
  }
  m_counter->increment(weight,incr++,"Trigger",trueTopo);

  // These jets have overlap removed
  std::vector<JetProxy> good_jets, bad_jets, b_jets, c_jets, good_jets_recl;
  std::vector<JetProxy> good_fat_jets, bad_fat_jets;
  std::vector<float> vD2, vD2_fat;
  std::vector<bool> visWmedium_fat ;
  if(! m_IsTruth){
    m_physobjsFiller->FillJetProxies(good_jets,bad_jets,b_jets);
    m_physobjsFiller->FillFatJetProxies(good_fat_jets,bad_fat_jets,vD2_fat,visWmedium_fat);

    if(m_doRecl) m_physobjsFiller->FillJetReclProxies(good_jets_recl,vD2);
  }
  if(m_IsTruth){
    m_physobjsFillerTruth->FillJetProxies(good_jets,b_jets);
  }

  double leptonPtCut = 25000.;
  if ( m_LowPtLepton ) leptonPtCut = 10000.;

  // isolated_xxx have overlap removed
  std::vector<ElectronProxy> baseline_electrons, isolated_baseline_electrons, isolated_signal_electrons;
  if(! m_IsTruth){
    m_physobjsFiller->FillElectronProxies(baseline_electrons, isolated_baseline_electrons, isolated_signal_electrons);
  }
  std::vector<ElectronTruthProxy> baseline_electrons_truth, isolated_baseline_electrons_truth, isolated_signal_electrons_truth;
  if(m_IsTruth){
    m_physobjsFillerTruth->FillElectronProxies(baseline_electrons_truth, isolated_baseline_electrons_truth, isolated_signal_electrons_truth);
  }
  // keep only signal electrons with Pt>25GeV for lepton triggers
  for ( std::vector<ElectronProxy>::iterator it = isolated_signal_electrons.begin();
	it != isolated_signal_electrons.end(); ) {
    if ( it->Pt() < leptonPtCut ) it = isolated_signal_electrons.erase(it);
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
    if ( it->Pt() < leptonPtCut ) it = isolated_signal_muons.erase(it);
    else it++;
  }
  // FIXME : trigger matching
  std::vector<MuonProxy> trigmatched_muons = isolated_signal_muons;

  std::vector<TauProxy> baseline_taus, signal_taus;
  if(! m_IsTruth){
    m_physobjsFiller->FillTauProxies(baseline_taus, signal_taus);
  }

  // Boson Tagging

  std::vector<float> vReclJetMass ;
  std::vector<float> vReclJetPt;
  std::vector<float> vReclJetEta;
  std::vector<float> vReclJetPhi;
  std::vector<bool> visWtight ;
  std::vector<bool> visWmedium ;
  std::vector<bool> visZtight ;
  std::vector<bool> visZmedium ;

  if(m_doRecl){
    for ( size_t j0=0; j0<good_jets_recl.size(); j0++){
      vReclJetMass.push_back(good_jets_recl[j0].M());
      vReclJetPt.push_back(good_jets_recl[j0].Pt());
      vReclJetEta.push_back(good_jets_recl[j0].Eta());
      vReclJetPhi.push_back(good_jets_recl[j0].Phi());

      float jetpt = good_jets_recl[j0].Pt();
      float jetm  = good_jets_recl[j0].M();
      float jetD2 = vD2.at(j0);

      BosonTagging BT;

      bool isWmedium = BT.ReturnTag(1,jetpt,jetm,jetD2);
      bool isWtight  = BT.ReturnTag(2,jetpt,jetm,jetD2);
      bool isZmedium = BT.ReturnTag(3,jetpt,jetm,jetD2);
      bool isZtight  = BT.ReturnTag(4,jetpt,jetm,jetD2);

      visWmedium.push_back(isWmedium);
      visWtight.push_back(isWtight);
      visZmedium.push_back(isZmedium);
      visZtight.push_back(isZtight);
    }
  }

  // missing ET
  TVector2* missingET = 0;
  TVector2* missingET_TST = 0;
  if(! m_IsTruth){
    if ( ! store->retrieve<TVector2>(missingET,"SUSYMET"+m_suffix+systag).isSuccess() ) throw std::runtime_error("could not retrieve SUSYMET"+m_suffix+systag);
    if ( ! store->retrieve<TVector2>(missingET_TST,"SUSYMETTST"+m_suffix+systag).isSuccess() ) throw std::runtime_error("could not retrieve SUSYMET"+m_suffix+systag);
  }
  if(m_IsTruth){
    if ( ! store->retrieve<TVector2>(missingET,"TruthMET"+m_suffix).isSuccess() ) throw std::runtime_error("could not retrieve TruthMET"+m_suffix);
  }
  double MissingEt = missingET->Mod();

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

  // select between electron and muon channel
  bool isonlymuonchannel = false;
  bool signalmuonistrigmatched = false;
  if ( !isolated_signal_muons.empty() ) { //FIXME replace by trigger condition
    isonlymuonchannel = true;
    if ( !isolated_signal_muons.empty() && !trigmatched_muons.empty() &&
	 isolated_signal_muons[0].muon() == trigmatched_muons[0].muon() ) {
      signalmuonistrigmatched = true;
    }
  }
  bool signalelecistrigmatched =false;
  if ( !isonlymuonchannel && !isolated_signal_electrons.empty() &&
       !trigmatched_electrons.empty() &&
       isolated_signal_electrons[0].electron() == trigmatched_electrons[0].electron() ) {
    signalelecistrigmatched = true;
  }
  bool oneLepton = false;
  int lep1Signal = 0;
  bool oneBaseLepton = false;
  TLorentzVector leptonTLV;
  int leptonCharge = 0;


  //if(m_IsTruth){ // GERALDINE
  //  signalelecistrigmatched = true;
  //  signalmuonistrigmatched = true;
  //}
  if ( true &&
       (( !m_IsTruth && isolated_baseline_muons.size()==1)      || (m_IsTruth && isolated_baseline_muons_truth.size()==1)) &&
       (( !m_IsTruth && isolated_signal_muons.size()==1)        || (m_IsTruth && isolated_signal_muons_truth.size()==1)) &&
       (( !m_IsTruth && trigmatched_muons.size()==1)            || m_IsTruth) &&
       (( !m_IsTruth && signalmuonistrigmatched)                || m_IsTruth) &&
       (( !m_IsTruth && isolated_baseline_electrons.size()==0 ) || (m_IsTruth && isolated_baseline_electrons_truth.size()==0)) ) {
    oneLepton = true;
    lep1Signal = 1;
    // to add for truth ?
    if(!m_IsTruth){
      leptonTLV = *(dynamic_cast<TLorentzVector*>(&(isolated_signal_muons[0])));
      leptonCharge = (int)(isolated_signal_muons[0].muon()->charge())*13;
    }
    if(m_IsTruth){
      leptonTLV = *(dynamic_cast<TLorentzVector*>(&(isolated_signal_muons_truth[0])));
      leptonCharge = (int)(isolated_signal_muons_truth[0].muontruth()->charge())*13;
    }
  }
  else if (  true &&
             (( !m_IsTruth && isolated_baseline_electrons.size()==1) || (m_IsTruth && isolated_baseline_electrons_truth.size()==1)) &&
	     (( !m_IsTruth && isolated_signal_electrons.size()==1)   || (m_IsTruth && isolated_signal_electrons_truth.size()==1))    &&
             (( !m_IsTruth && trigmatched_electrons.size()==1)       || m_IsTruth ) &&
	     (( !m_IsTruth && signalelecistrigmatched)               || m_IsTruth) &&
             (( !m_IsTruth && isolated_baseline_muons.size()==0)     || (m_IsTruth && isolated_baseline_muons_truth.size()==0))) {
    oneLepton = true;
    lep1Signal = 1;
    if(!m_IsTruth){
      leptonTLV = *(dynamic_cast<TLorentzVector*>(&(isolated_signal_electrons[0])));
      leptonCharge = (int)(isolated_signal_electrons[0].electron()->charge())*11;
    }
    if(m_IsTruth){
      leptonTLV = *(dynamic_cast<TLorentzVector*>(&(isolated_signal_electrons_truth[0])));
      leptonCharge = (int)(isolated_signal_electrons_truth[0].eltruth()->charge())*11;
    }
  }
  if ( !oneLepton ){
    if( m_doFake ){
      if ( true &&
	   ( isolated_baseline_muons.size()==1 &&  isolated_baseline_electrons.size()==0) ){
	oneBaseLepton = true;
	leptonTLV = *(dynamic_cast<TLorentzVector*>(&(isolated_baseline_muons[0])));
	leptonCharge = (int)(isolated_baseline_muons[0].muon()->charge())*13;
      }
      else if (  true &&
		 ( isolated_baseline_electrons.size()==1  && isolated_baseline_muons.size()==0)){
	oneBaseLepton = true;
	leptonTLV = *(dynamic_cast<TLorentzVector*>(&(isolated_baseline_electrons[0])));
	leptonCharge = (int)(isolated_baseline_electrons[0].electron()->charge())*11;
      }
      if( !oneBaseLepton ){
	return true;
      }
    }else{
      return true;
    }
  }
  m_counter->increment(weight,incr++,"1 Lepton",trueTopo);

  // lepton isolation

  float topoetcone20 = 0 ;
  float ptvarcone30 = 0 ;
  float ptvarcone20 = 0 ;

  if(!m_IsTruth){
    if(true && isolated_signal_electrons.size()==1){
      topoetcone20 = (isolated_signal_electrons[0].electron())->isolation(xAOD::Iso::topoetcone20);
      ptvarcone30  = (isolated_signal_electrons[0].electron())->isolation(xAOD::Iso::ptvarcone30);
      ptvarcone20  = (isolated_signal_electrons[0].electron())->isolation(xAOD::Iso::ptvarcone20);
    }

    if(true && isolated_signal_muons.size()==1){
      topoetcone20 = (isolated_signal_muons[0].muon())->isolation(xAOD::Iso::topoetcone20);
      ptvarcone30  = (isolated_signal_muons[0].muon())->isolation(xAOD::Iso::ptvarcone30);
      ptvarcone20  = (isolated_signal_muons[0].muon())->isolation(xAOD::Iso::ptvarcone20);
    }
  }

  // calculate fake lepton weight
  //
  std::vector<float> fake_weight(3,0.);
  if( m_doFake ){
    Double_t fakeEff = 0.;
    Double_t fakeEffUp = 0.;
    Double_t fakeEffDown = 0.;
    Double_t realEff = 0.;
    Double_t signalWeight = 0.;
    Double_t baseWeight = 1.;
    if(true && isolated_baseline_electrons.size()==1){
      fakeEff = m_elecFakeEff->getEfficiency( leptonTLV.Pt(), leptonTLV.Eta(), 0);
      fakeEffUp = m_elecFakeEff->getEfficiency( leptonTLV.Pt(), leptonTLV.Eta(), 1);
      fakeEffDown = m_elecFakeEff->getEfficiency( leptonTLV.Pt(), leptonTLV.Eta(), 2);
      realEff = m_elecRealEff->getEfficiency( leptonTLV.Pt(), leptonTLV.Eta(), 0);
      baseWeight *= prescaleEl;
      if( isolated_signal_electrons.size()==1 ) signalWeight = 1.;
    }
    else if(true && isolated_baseline_muons.size()==1){
      fakeEff = m_elecFakeEff->getEfficiency( leptonTLV.Pt(), leptonTLV.Eta(), 0);
      fakeEffUp = m_elecFakeEff->getEfficiency( leptonTLV.Pt(), leptonTLV.Eta(), 1);
      fakeEffDown = m_elecFakeEff->getEfficiency( leptonTLV.Pt(), leptonTLV.Eta(), 2);
      realEff = m_elecRealEff->getEfficiency( leptonTLV.Pt(), leptonTLV.Eta(), 0);
      baseWeight *= prescaleMu;
      if( isolated_signal_muons.size()==1 ) signalWeight = 1.;
    }
    fake_weight.at(0) = Nfake(realEff, fakeEff, signalWeight, baseWeight-signalWeight );
    fake_weight.at(1) = Nfake(realEff, fakeEffUp, signalWeight, baseWeight-signalWeight );
    fake_weight.at(2) = Nfake(realEff, fakeEffDown, signalWeight, baseWeight-signalWeight );
  }

  // Apply Lepton scale factors
  float muSF = eventInfo->auxdecor<float>("muSF");
  if ( muSF != 0.f ) weight *= muSF;
  float elSF = eventInfo->auxdecor<float>("elSF");
  if ( elSF != 0.f ) weight *= elSF;

  // Add lepton to jets (SR) or MET (VR)
  TVector2 missingETPrime =  *missingET;
  if ( m_isVR ) {
    missingETPrime = missingETPrime + TVector2(leptonTLV.Px(),leptonTLV.Py());
  }
  else {
    good_jets.push_back(JetProxy(leptonTLV,true,true,true,true,false));
    std::sort(good_jets.begin(),good_jets.end(),PtOrder<JetProxy>);
  }
  double MissingEtPrime = missingETPrime.Mod();

  // MT cut
  double mt = std::sqrt( 2.*leptonTLV.Pt()*MissingEt *
			 (1.-(leptonTLV.Px()*missingET->Px() + leptonTLV.Py()*missingET->Py())/(leptonTLV.Pt()*MissingEt)) );
  if(!(mt >30000 && mt<100000)) return true;
  m_counter->increment(weight,incr++,"MT cut",trueTopo);

  if (good_jets.size()<1) return true;
  m_counter->increment(weight,incr++,"At least one jet",trueTopo);

  // jet timing cut
  std::vector<float> time;
  m_proxyUtils.EnergyWeightedTime(good_jets,time);

  // MissingET cut
  if (!(MissingEtPrime > m_cutVal.m_cutEtMiss)) return true;
  m_counter->increment(weight,incr++,"MET cut",trueTopo);

  // Leading jet Pt cut
  if (!(good_jets[0].Pt() > m_cutVal.m_cutJetPt0)) return true;
  m_counter->increment(weight,incr++,"1 jet Pt > 130 GeV Selection",trueTopo);

  // leave counter to keep same cutflow
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
  //float HT = meffincl -  MissingEt;

  // Sherpa MassiveCB W/Z reweighting : not implemented yet in SUSYOBJDef_xAOD
  float WZweight = 1.;
  if( !m_IsData ){
    WZweight = eventInfo->auxdecor<float>("WZweight");
  }

  // Fill ISR variables. These vectors have for each jet > 50GeV the ISR variables in them.
  std::vector<size_t> isr_jet_indices;
  std::vector<std::vector<double> > ISRvars;
  std::vector<double> ISRvar_alpha;
  m_proxyUtils.GetAlphaISRVar(good_jets,MissingEtPrime,ISRvar_alpha);
  std::vector<double> ISRvar_minPtDistinction;
  m_proxyUtils.GetMinPtDistinctionISR(good_jets,ISRvar_minPtDistinction);
  std::vector<double> ISRvar_minDeltaFraction;
  m_proxyUtils.GetMinDeltaFraction(good_jets,ISRvar_minDeltaFraction);
  std::vector<double> ISRvar_minRapidityGap;
  m_proxyUtils.GetMinRapidityGap(good_jets,ISRvar_minRapidityGap);
  std::vector<double> ISRvar_maxRapidityOtherJets;
  m_proxyUtils.GetMaxRapidityOtherJets(good_jets,ISRvar_maxRapidityOtherJets);
  std::vector<double> ISRvar_dPhiJetMET;
  m_proxyUtils.GetdPhiJetMet(good_jets,phi_met,ISRvar_dPhiJetMET);
  m_proxyUtils.GetISRJet(good_jets,isr_jet_indices,MissingEtPrime,phi_met,"squark",false);
  ISRvars.push_back(ISRvar_alpha);
  ISRvars.push_back(ISRvar_minPtDistinction);
  ISRvars.push_back(ISRvar_minDeltaFraction);
  ISRvars.push_back(ISRvar_minRapidityGap);
  ISRvars.push_back(ISRvar_maxRapidityOtherJets);
  ISRvars.push_back(ISRvar_dPhiJetMET);
  std::vector<JetProxy> nonISR_jets = good_jets;
  if (isr_jet_indices.size()==1) {
    nonISR_jets.erase(nonISR_jets.begin()+isr_jet_indices[0],nonISR_jets.end());
  }

  /* not used now
  double mT2=-9;
  if (good_jets.size()>=2) mT2 = m_proxyUtils.MT2(good_jets,missingETPrime);
  double mT2_noISR=-9;
  if (nonISR_jets.size()>=2) mT2_noISR = m_proxyUtils.MT2(nonISR_jets,missingETPrime);
  out() << " mT2 " << mT2 << " " << mT2_noISR << std::endl;
  */

  std::map<TString,float> RJigsawVariables;
  if (  m_fillTRJigsawVars ) {
    m_proxyUtils.CalculateRJigsawVariables(good_jets,
					   missingETPrime.X(),
					   missingETPrime.Y(),
					   RJigsawVariables,
             			m_cutVal.m_cutRJigsawJetPt);
  }


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

  CleaningHelper cleaningHelper;

  if(m_doSmallNtuple) {
    unsigned int runnum = RunNumber;
    if ( ! m_IsData  && ! m_IsTruth) runnum = mc_channel_number;

    std::vector<float> jetSmearSystW;

    // other cleaning tests
    if(!m_IsTruth){

      // bad jet veto
      if ( !bad_jets.empty() ) cleaningHelper.cleaning.at("badJetVeto") = true;

      // bad muon veto
      for ( size_t i = 0; i < isolated_baseline_muons.size(); i++) {
	if ( isolated_baseline_muons[i].passOVerlapRemoval() &&
	     isolated_baseline_muons[i].isBad() ) {
	  cleaningHelper.cleaning.at("badMuonVeto") = true;
	  break;
	}
      }
      // Cosmic muon cut
      if ( m_proxyUtils.CosmicMuon(isolated_baseline_muons) )  cleaningHelper.cleaning.at("cosmicMuonVeto") = true;

      // bad muons for MET cut: based on non isolated muons
      if (  m_proxyUtils.isbadMETmuon(baseline_muons, MissingEt, *missingET) ) cleaningHelper.cleaning.at("badMetMuonVeto") = true;

      // bad Tile cut
      if ( !m_IsTruth && m_proxyUtils.badTileVeto(good_jets,*missingET)) cleaningHelper.cleaning.at("badTileVeto") = true;

      // average timing of 2 leading jets
      if (fabs(time[0]) > 5) cleaningHelper.cleaning.at("leadingJetTimingVeto") = true;
      // FIXME why not in CRWT ?
      //bool chfTileVeto =  m_proxyUtils.chfTileVeto(good_jets);
      //if ( chfTileVeto ) cleaning += 4;

      if ( m_proxyUtils.chfVeto(good_jets)) cleaningHelper.cleaning.at("chfVeto") = true;;

      bool * failMetCleaning = nullptr;
      if ( !store->retrieve<bool>(failMetCleaning,"failMetCleaning").isSuccess() ) throw std::runtime_error("could not retrieve failMetCleaning");
      if ( *failMetCleaning) cleaningHelper.cleaning.at("metTSTCleaningVeto") = true;
    }

    float dPhiBadTile = m_proxyUtils.dPhiBadTile(good_jets,*missingET);

    bool isNCBEvent = false;
    if ( m_IsData ) {
      bool* NCBEventFlag = 0;
      if ( !store->retrieve<bool>(NCBEventFlag,"NCBEventFlag").isSuccess() ) throw std::runtime_error("could not retrieve NCBEventFlag");
      isNCBEvent = *NCBEventFlag;
    }

    unsigned long const cleaning = cleaningHelper.finalCleaning();

    m_proxyUtils.FillNTVars(m_ntv, runnum, EventNumber, LumiBlockNumber, veto, weight, normWeight, *pileupWeights, PRWHash, genWeight,ttbarWeightHT,ttbarWeightPt2,ttbarAvgPt,WZweight, b_jets.size(), c_jets.size(), MissingEtPrime, phi_met, missingET_TST->Mod(), missingET_TST->Phi(), Meff, meffincl, minDphi, RemainingminDPhi, good_jets, good_fat_jets, vD2_fat, visWmedium_fat, trueTopo, cleaning, time[0],jetSmearSystW,0, 0.,dPhiBadTile,isNCBEvent,m_IsTruth,baseline_taus,signal_taus);

    if ( systag == ""  && !m_IsTruth) {
      std::vector<float>* p_systweights = 0;
      if ( ! store->retrieve(p_systweights,"event_weights"+m_suffix).isSuccess() ) throw std::runtime_error("Could not retrieve event_weights"+m_suffix);
      m_ntv.systWeights = *p_systweights;

      std::vector<float>* p_btagSystweights = 0;
      if ( ! store->retrieve(p_btagSystweights,"btag_weights"+m_suffix).isSuccess() ) throw std::runtime_error("Could not retrieve btag_weights"+m_suffix);
      m_ntv.btagSystWeights = *p_btagSystweights;
    }

    if( !m_IsTruth && m_fillReclusteringVars){
      m_proxyUtils.FillNTReclusteringVars(m_RTntv, good_jets,vReclJetMass,vReclJetPt,vReclJetEta,vReclJetPhi,vD2,visWmedium, visWtight, visZmedium, visZtight);
    }

    FillCRWTVars(m_crwtntv,leptonTLV,*missingET,leptonCharge,ptvarcone20,ptvarcone30,topoetcone20,lep1Signal,fake_weight);

    m_proxyUtils.FillNTExtraVars(m_extrantv, MET_Track, MET_Track_phi, Ap);

    if (  m_fillTRJigsawVars ) m_proxyUtils.FillNTRJigsawVars(m_rjigsawntv, RJigsawVariables );

    m_tree->Fill();
  }
  return true;
}

void ZeroLeptonCRWT::finish()
{
  if ( m_DoSystematics ) {
    out() << m_counterRepository << std::endl;
  }
  else {
    out() << *m_counter << std::endl;
  }
}


void ZeroLeptonCRWT::FillCRWTVars(NTCRWTVars& crwtvars, const TLorentzVector& lepton, const TVector2& metv, int lepsign,
				  float ptvarcone20, float ptvarcone30, float topoetcone20, int lep1Signal, const std::vector<float> fakeWeight)
{
  crwtvars.Reset();
  crwtvars.lep1sign = lepsign;
  crwtvars.lep1Pt  = lepton.Pt() * 0.001;
  crwtvars.lep1Eta = lepton.Eta();
  crwtvars.lep1Phi = lepton.Phi();
  crwtvars.lep1ptvarcone20 = ptvarcone20 * 0.001;
  crwtvars.lep1ptvarcone30 = ptvarcone30 * 0.001 ;
  crwtvars.lep1topoetcone20 = topoetcone20 * 0.001 ;
  crwtvars.lep1Signal = lep1Signal;

  double met = std::sqrt(metv.Px()*metv.Px()+metv.Py()*metv.Py());
  crwtvars.mt = std::sqrt(2.*lepton.Pt() * met * (1. - (lepton.Px()*metv.Px() + lepton.Py()*metv.Py())/(lepton.Pt()*met))) * 0.001;

  TLorentzVector metTLV;
  metTLV.SetPxPyPzE(metv.Px(),metv.Py(),0,met);
  crwtvars.dphilMET = lepton.DeltaPhi(metTLV);

  crwtvars.Weta = (lepton+metTLV).Eta() ;

  double wpx = lepton.Px()+metv.Px();
  double wpy = lepton.Py()+metv.Py();
  crwtvars.Wpt = std::sqrt(wpx*wpx+wpy*wpy) * 0.001;

  if(fakeWeight.size()>0) crwtvars.fakeWeight = fakeWeight.at(0);
  if(fakeWeight.size()>1) crwtvars.fakeWeightUp = fakeWeight.at(1);
  if(fakeWeight.size()>2) crwtvars.fakeWeightDown = fakeWeight.at(2);
}

ClassImp(ZeroLeptonCRWT);

