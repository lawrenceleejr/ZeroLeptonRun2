
#include "ZeroLeptonRun2/ZeroLeptonCR3L.h"
#include "ZeroLeptonRun2/PhysObjProxies.h"
#include "ZeroLeptonRun2/PtOrder.h"

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

ZeroLeptonCR3L::ZeroLeptonCR3L(const char *name)
  : cafe::Processor(name),
    m_tree(0),
    m_stringRegion("CR3L_SRAll"),
    m_doSmallNtuple(true),
    m_fillTRJigsawVars(true),
    m_fillReclusteringVars(true),
    m_IsData(false),
    m_IsSignal(false),
    m_IsTruth(false),
    m_isVR(false),
    m_doRecl(false),
    m_DoSystematics(false),
    m_period(INVALID),
    m_isMuonChannel(false),
    m_isElectronChannel(false),
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
  m_doRecl = config.get("doRecl",false);
  m_suffixRecl = config.get("suffixRecl","");
  m_IsTruth = config.get("IsTruth",false);
  m_DoSystematics = config.get("DoSystematics",false);
  m_doFake = config.get("doFake", false);

  m_period = periodFromString(config.get("Period","p13tev"));
  if ( m_period == p7tev ) throw(std::domain_error("ZeroLeptonCR3L does not support the 7tev run period"));
  if ( m_period == INVALID ) throw(std::domain_error("ZeroLeptonCR3L: invalid run period specified"));

  if ( m_IsData && m_period == p8tev ) {
    m_isMuonChannel = config.get("IsMuonChannel",false);
    m_isElectronChannel = config.get("IsElectronChannel",false);;
  }
  else {
    m_isMuonChannel = true;
    m_isElectronChannel = true;
  }

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
  if ( m_isVR ) m_stringRegion = "VR3L_SRAll";

  std::string cutfile = config.get("cutfile","None");
  if ( cutfile == "None" ) throw(std::domain_error("ZeroLeptonCR3L: invalid cut file specified"));
  m_cutVal.ReadCutValues(cutfile);

  m_suffix = config.get("suffix","");
  m_suffixRecl = config.get("suffixRecl","");
  m_suffixSyst = "test";
  m_physobjsFiller = new PhysObjProxyFiller(20000.f,10000.f,10000.f,25000.f,m_suffix,m_doRecl,m_suffixRecl,m_suffixSyst);
  m_physobjsFillerTruth = new PhysObjProxyFillerTruth(20000.f,10000.f,10000.f,25000.f,m_suffix);
  m_proxyUtils = PhysObjProxyUtils(m_IsData);

  m_ZLUtils = ZeroLeptonUtils(m_IsData);
}

ZeroLeptonCR3L::~ZeroLeptonCR3L()
{
  if ( !m_DoSystematics && m_counter ) delete m_counter;
  if ( m_physobjsFiller ) delete m_physobjsFiller;
  if ( m_physobjsFillerTruth ) delete m_physobjsFillerTruth;
}

TTree* ZeroLeptonCR3L::bookTree(const std::string& treename)
{
  const char* name(treename.c_str());
  TTree* tree = new TTree(name,"ZeroLepton final optimisation");
  tree->SetDirectory(getDirectory());
  bookNTVars(tree,m_ntv,false);
  bookNTExtraVars(tree,m_extrantv);
  if ( m_fillTRJigsawVars) bookNTRJigsawVars(tree,m_rjigsawntv);
  if ( m_fillReclusteringVars) bookNTReclusteringVars(tree,m_RTntv);
  bookNTCR3LVars(tree,m_cr3lntv);
  return tree;
}

TTree* ZeroLeptonCR3L::getTree(const std::string& treename)
{
  std::map<std::string,TTree*>::const_iterator pos = m_treeRepository.find(treename);
  if ( pos == m_treeRepository.end() ) {
    pos = m_treeRepository.insert(std::make_pair(treename,bookTree(treename))).first;
  }
  return pos->second;
}

void ZeroLeptonCR3L::begin()
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


bool ZeroLeptonCR3L::processEvent(xAOD::TEvent& event)
{
  // access the transient store
  xAOD::TStore* store = xAOD::TActiveStore::store();
  std::string systag = "";
  if ( m_DoSystematics ) {
    CP::SystematicSet* currentSyst = 0;
    if ( ! store->retrieve(currentSyst, "CurrentSystematicSet").isSuccess() ) throw std::runtime_error("Could not retrieve CurrentSystematicSet");
    std::string sysname = currentSyst->name();
    if (sysname != "" ) systag = "_"+sysname+"_";
    m_counter = m_counterRepository.counter(sysname);
    if (sysname == "" ) {
      m_tree = getTree(m_stringRegion+"NT");
    }
    else {
      m_tree = getTree(m_stringRegion+"NT_"+sysname);
      m_physobjsFiller->setSuffix(m_suffix+systag);
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
  if ( !m_IsData && !m_IsTruth  ) {
    if ( !store->retrieve< std::vector<float> >(pileupWeights,"pileupWeights").isSuccess() ) throw std::runtime_error("could not retrieve pileupWeights");
    //out() << " pileup weight " << (*pileupWeights)[0] << std::endl;
    //weight *= (*pileupWeights)[0];
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

  if(! m_IsTruth){
    bool passEltrigger=false;
    bool passMutrigger=false;
    if( (int)eventInfo->auxdata<char>("HLT_e60_lhmedium")==1  ||
	(int)eventInfo->auxdata<char>("HLT_e120_lhloose")==1     )
      passEltrigger = true;
    if ( m_IsData ) {
      if ( (int)eventInfo->auxdata<char>("HLT_e24_lhmedium_L1EM20VH")==1 ) passEltrigger = true;
    }
    else {
      if ( (int)eventInfo->auxdata<char>("HLT_e24_lhmedium_L1EM18VH")==1 ) passEltrigger = true;
    }
    if((int)eventInfo->auxdata<char>("HLT_mu20_iloose_L1MU15")==1 ||
       (int)eventInfo->auxdata<char>("HLT_mu50")==1)
      passMutrigger = true;
    if( !(passEltrigger || passMutrigger) ) return true;
  }

  m_counter->increment(weight,incr++,"Trigger",trueTopo);

  // These jets have overlap removed
  std::vector<JetProxy> good_jets, bad_jets, b_jets, c_jets;
  std::vector<JetProxy> good_fat_jets, bad_fat_jets;
  std::vector<float> vD2, vD2_fat;
  std::vector<bool> visWmedium_fat ;
  if(! m_IsTruth){
    m_physobjsFiller->FillJetProxies(good_jets,bad_jets,b_jets);
    m_physobjsFiller->FillFatJetProxies(good_fat_jets,bad_fat_jets,vD2_fat,visWmedium_fat);
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
  int nel=0;
  int neltruth=0;
  for ( std::vector<ElectronProxy>::iterator it = isolated_signal_electrons.begin();
  	it != isolated_signal_electrons.end(); ) {
    if ( it->Pt() < 25000. && nel<1 ) it = isolated_signal_electrons.erase(it);
    else it++;
    nel++;
  }
  for ( std::vector<ElectronTruthProxy>::iterator itt = isolated_signal_electrons_truth.begin();
        itt != isolated_signal_electrons_truth.end(); ) {
    if ( itt->Pt() < 25000. && neltruth<1 ) itt = isolated_signal_electrons_truth.erase(itt);
    else itt++;
    neltruth++;
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
  int nmu=0;
  int nmutruth=0;
  for ( std::vector<MuonProxy>::iterator it = isolated_signal_muons.begin();
	it != isolated_signal_muons.end(); ) {
    if ( it->Pt() < 25000. && nmu<1 ) it = isolated_signal_muons.erase(it);
    else it++;
    nmu++;
  }
  for ( std::vector<MuonTruthProxy>::iterator itt = isolated_signal_muons_truth.begin();
        itt != isolated_signal_muons_truth.end(); ) {
    if ( itt->Pt() < 25000. && nmutruth<1 ) itt = isolated_signal_muons_truth.erase(itt);
    else itt++;
    nmutruth++;
  }
  // FIXME : trigger matching
  std::vector<MuonProxy> trigmatched_muons = isolated_signal_muons;

  std::vector<TauProxy> baseline_taus, signal_taus;
  if(! m_IsTruth){
    m_physobjsFiller->FillTauProxies(baseline_taus, signal_taus);
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
  //double MissingEt = missingET->Mod();

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

  // FIXME fake lepton bkg estimate
  // do we need to re-implement that ?

  std::vector<TLorentzVector> leptonTLVs;
  std::vector<int> leptonCharges;
  std::vector<int> leptonSignal(3,1);

  int nW=1000;

  if ( m_isMuonChannel &&
       ( (!m_isVR && !m_IsTruth && (isolated_signal_muons.size()==3 || (isolated_signal_muons.size()==2 && isolated_signal_electrons.size()==1)) )
	 || ( (m_isVR || m_doFake) && !m_IsTruth && isolated_signal_muons.size()==2 && isolated_signal_electrons.size()==0 && (isolated_baseline_muons.size()==3 || isolated_baseline_electrons.size()==1) )
	 || (!m_isVR && m_IsTruth && (isolated_signal_muons_truth.size()==3 || (isolated_signal_muons_truth.size()==2 && isolated_signal_electrons_truth.size()==1)) )
	 || ( m_isVR && m_IsTruth &&  isolated_signal_muons_truth.size()==2 && isolated_signal_electrons_truth.size()==0 && (isolated_baseline_muons_truth.size()==3 || isolated_baseline_electrons_truth.size()==1) ))){

    if(!m_IsTruth){
      leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_muons[0]))));
      leptonCharges.push_back((int)(isolated_signal_muons[0].muon()->charge())*13);
      leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_muons[1]))));
      leptonCharges.push_back((int)(isolated_signal_muons[1].muon()->charge())*13);

      if(isolated_signal_muons.size()==3){
	leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_muons[2]))));
	leptonCharges.push_back((int)(isolated_signal_muons[2].muon()->charge())*13);
      } // 3 signal muons
      else if(isolated_signal_muons.size()==2){
	if(isolated_signal_electrons.size()==1){
	  leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_electrons[0]))));
	  leptonCharges.push_back((int)(isolated_signal_electrons[0].electron()->charge())*11);
	} // not validation region
	else if( m_isVR || m_doFake ){
	  if(isolated_baseline_muons.size()==3){
	    for(int i=0;i<3;i++){
	      if(isolated_baseline_muons[i].Pt()==isolated_signal_muons[0].Pt() || isolated_baseline_muons[i].Pt()==isolated_signal_muons[1].Pt())
		continue;
	      leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_baseline_muons[i]))));
	      leptonCharges.push_back((int)(isolated_baseline_muons[i].muon()->charge())*13);
	      nW = i ;
	      leptonSignal.at(2) = 0;
	    }
	  } // 3 baseline m
	  else if(isolated_baseline_electrons.size()==1){
	    leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_baseline_electrons[0]))));
	    leptonCharges.push_back((int)(isolated_baseline_electrons[0].electron()->charge())*11);
	    leptonSignal.at(2) = 0;
	  } // 1 baseline e
	} // validation region
      } // 2 signal muons
    } // not truth
    if(m_IsTruth){
      leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_muons_truth[0]))));
      leptonCharges.push_back((int)(isolated_signal_muons_truth[0].muontruth()->charge())*13);
      leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_muons_truth[1]))));
      leptonCharges.push_back((int)(isolated_signal_muons_truth[1].muontruth()->charge())*13);

      if(isolated_signal_muons_truth.size()==3){
	leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_muons_truth[2]))));
	leptonCharges.push_back((int)(isolated_signal_muons_truth[2].muontruth()->charge())*13);
      } // 3 signal muons
      if(isolated_signal_muons_truth.size()==2){
	if(!m_isVR){
	  leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_electrons_truth[0]))));
	  leptonCharges.push_back((int)(isolated_signal_electrons_truth[0].eltruth()->charge())*11);
	} // not validation region
	if( m_isVR){
	  if(isolated_baseline_muons_truth.size()==3){
            for(int i=0;i<3;i++){
              if(isolated_baseline_muons_truth[i].Pt()==isolated_signal_muons_truth[0].Pt() || isolated_baseline_muons_truth[i].Pt()==isolated_signal_muons_truth[1].Pt())
                continue;
              leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_baseline_muons_truth[i]))));
              leptonCharges.push_back((int)(isolated_baseline_muons_truth[i].muontruth()->charge())*13);
	      nW=i;
            }
          } // 3 baseline mu
          if(isolated_baseline_electrons_truth.size()==1){
            leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_baseline_electrons_truth[0]))));
            leptonCharges.push_back((int)(isolated_baseline_electrons_truth[0].eltruth()->charge())*11);
          } // 1 baseline e

	} // validation region
      } // 2 signal muons
    } // m_IsTruth
  }

  // ELECTRONS

  if ( m_isElectronChannel &&
       (  (!m_isVR && !m_IsTruth && (isolated_signal_electrons.size()==3 || (isolated_signal_electrons.size()==2 && isolated_signal_muons.size()==1)) )
	  || ( (m_isVR || m_doFake) && !m_IsTruth && isolated_signal_electrons.size()==2 && isolated_signal_muons.size()==0 && (isolated_baseline_electrons.size()==3 || isolated_baseline_muons.size()==1) )
	  || (!m_isVR && m_IsTruth && (isolated_signal_electrons_truth.size()==3 || (isolated_signal_electrons_truth.size()==2 && isolated_signal_muons_truth.size()==1)) )
	  || ( m_isVR && m_IsTruth &&  isolated_signal_electrons_truth.size()==2 && isolated_signal_muons_truth.size()==0 && (isolated_baseline_electrons_truth.size()==3 || isolated_baseline_muons_truth.size()==1) ))){

    if(!m_IsTruth){
      leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_electrons[0]))));
      leptonCharges.push_back((int)(isolated_signal_electrons[0].electron()->charge())*11);
      leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_electrons[1]))));
      leptonCharges.push_back((int)(isolated_signal_electrons[1].electron()->charge())*11);

      if(isolated_signal_electrons.size()==3){
        leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_electrons[2]))));
        leptonCharges.push_back((int)(isolated_signal_electrons[2].electron()->charge())*11);
      } // 3 signal e
      else if(isolated_signal_electrons.size()==2){
        if(isolated_signal_muons.size()==1){
          leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_muons[0]))));
          leptonCharges.push_back((int)(isolated_signal_muons[0].muon()->charge())*13);
        } // not validation region
        else if( m_isVR || m_doFake){
          if(isolated_baseline_electrons.size()==3){
            for(int i=0;i<3;i++){
              if(isolated_baseline_electrons[i].Pt()==isolated_signal_electrons[0].Pt() || isolated_baseline_electrons[i].Pt()==isolated_signal_electrons[1].Pt())
                continue;
              leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_baseline_electrons[i]))));
              leptonCharges.push_back((int)(isolated_baseline_electrons[i].electron()->charge())*11);
	      nW=i;
	      leptonSignal.at(2) = 0;
            }
          } // 3 baseline e
          else if(isolated_baseline_muons.size()==1){
            leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_baseline_muons[0]))));
            leptonCharges.push_back((int)(isolated_baseline_muons[0].muon()->charge())*13);
	    leptonSignal.at(2) = 0;
          } // 1 baseline m
	} // validation
      } // 2 signal e
    } // not truth
    if(m_IsTruth){
      leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_electrons_truth[0]))));
      leptonCharges.push_back((int)(isolated_signal_electrons_truth[0].eltruth()->charge())*11);
      leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_electrons_truth[1]))));
      leptonCharges.push_back((int)(isolated_signal_electrons_truth[1].eltruth()->charge())*11);

      if(isolated_signal_electrons_truth.size()==3){
        leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_electrons_truth[2]))));
        leptonCharges.push_back((int)(isolated_signal_electrons_truth[2].eltruth()->charge())*11);
      } // 3 signal electrons
      if(isolated_signal_electrons_truth.size()==2){
        if(!m_isVR){
          leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_muons_truth[0]))));
          leptonCharges.push_back((int)(isolated_signal_muons_truth[0].muontruth()->charge())*13);
        }// not validation region
        if( m_isVR){
          if(isolated_baseline_electrons_truth.size()==3){
            for(int i=0;i<3;i++){
              if(isolated_baseline_electrons_truth[i].Pt()==isolated_signal_electrons_truth[0].Pt() || isolated_baseline_electrons_truth[i].Pt()==isolated_signal_electrons_truth[1].Pt())
                continue;
              leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_baseline_electrons_truth[i]))));
              leptonCharges.push_back((int)(isolated_baseline_electrons_truth[i].eltruth()->charge())*11);
	      nW=i;
            }
          } // 3 baseline electrons
          if(isolated_baseline_muons_truth.size()==1){
            leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_baseline_muons_truth[0]))));
            leptonCharges.push_back((int)(isolated_baseline_muons_truth[0].muontruth()->charge())*13);
          } // 1 baseline muon
	} // validation region
      } // 2 signal e
    } // isTruth
  } // Electrons


  if ( leptonTLVs.empty() ) return true;

  // calculate fake lepton weight
  //
  std::vector<float> fake_weight(3,0.);
  if( m_doFake ){
    std::vector<Double_t> fakeEff(3, 0.);
    std::vector<Double_t> fakeEffUp(3, 0.);
    std::vector<Double_t> fakeEffDown(3, 0.);
    std::vector<Double_t> realEff(3, 0.);
    Int_t nt(0);
    Double_t nttt(0), nttl(0), ntlt(0), nltt(0);

    for(UInt_t i=0; i<leptonTLVs.size() && i<leptonCharges.size() && i<leptonSignal.size(); i++){
      if( abs(leptonCharges.at(i))==11 ){
	fakeEff.at(i) = m_elecFakeEff->getEfficiency( leptonTLVs.at(i).Pt(), leptonTLVs.at(i).Eta(), 0);
	fakeEffUp.at(i) = m_elecFakeEff->getEfficiency( leptonTLVs.at(i).Pt(), leptonTLVs.at(i).Eta(), 1);
	fakeEffDown.at(i) = m_elecFakeEff->getEfficiency( leptonTLVs.at(i).Pt(), leptonTLVs.at(i).Eta(), 2);
	realEff.at(i) = m_elecRealEff->getEfficiency( leptonTLVs.at(i).Pt(), leptonTLVs.at(i).Eta(), 0);
      }else if( abs(leptonCharges.at(i))==13 ){
	fakeEff.at(i) = m_muonFakeEff->getEfficiency( leptonTLVs.at(i).Pt(), leptonTLVs.at(i).Eta(), 0);
	fakeEffUp.at(i) = m_muonFakeEff->getEfficiency( leptonTLVs.at(i).Pt(), leptonTLVs.at(i).Eta(), 1);
	fakeEffDown.at(i) = m_muonFakeEff->getEfficiency( leptonTLVs.at(i).Pt(), leptonTLVs.at(i).Eta(), 2);
	realEff.at(i) = m_muonRealEff->getEfficiency( leptonTLVs.at(i).Pt(), leptonTLVs.at(i).Eta(), 0);
      }
      if( leptonSignal.at(i)==1 ) nt++;
    }
    if( nt==3 ) nttt = 1.;
    else if( nt==2 ){
      if( leptonSignal.at(0)==0 ) nltt = 1.;
      else if( leptonSignal.at(1)==0 ) ntlt = 1.;
      else if( leptonSignal.at(2)==0 ) nttl = 1.;
    }
    fake_weight.at(0) = Nrrf_rfr_frr(realEff.at(0), realEff.at(1), realEff.at(2), fakeEff.at(0), fakeEff.at(1), fakeEff.at(2), nttt, nttl, ntlt, nltt);
    fake_weight.at(1) = Nrrf_rfr_frr(realEff.at(0), realEff.at(1), realEff.at(2), fakeEffUp.at(0), fakeEffUp.at(1), fakeEffUp.at(2), nttt, nttl, ntlt, nltt);
    fake_weight.at(2) = Nrrf_rfr_frr(realEff.at(0), realEff.at(1), realEff.at(2), fakeEffDown.at(0), fakeEffDown.at(1), fakeEffDown.at(2), nttt, nttl, ntlt, nltt);
  }

  int lepfromW = 0 ;
  //float lepptfromW = 0 ;

  TLorentzVector dileptonTLV ;
  TLorentzVector oneleptonTLV;
  // CASE WITH THREE IDENTICAL LEPTONS - not needed for VR, as we have only two OS signal leptons
  if( !m_isVR && ((!m_IsTruth && ((isolated_signal_electrons.size()==3 || isolated_signal_muons.size()==3) || (m_doFake && (isolated_baseline_electrons.size()==3 || isolated_baseline_muons.size()==3))))
		  || (m_IsTruth && (isolated_signal_electrons_truth.size()==3 || isolated_signal_muons_truth.size()==3)))){

    // Remove e+e+e+ and e-e-e- events, not consistent with Z(ee) boson!
    if(fabs(leptonCharges[0]+leptonCharges[1]+leptonCharges[2])==3) return true;

    double mass1=0;
    double mass2=0;
    double mass3=0;
    double massZ = 91187.6;

    if(leptonCharges[0]*leptonCharges[1]<0){
      mass1 = (leptonTLVs[0]+leptonTLVs[1]).M();
    }
    if(leptonCharges[0]*leptonCharges[2]<0){
      mass2 = (leptonTLVs[0]+leptonTLVs[2]).M();
    }
    if(leptonCharges[1]*leptonCharges[2]<0){
      mass3 = (leptonTLVs[1]+leptonTLVs[2]).M();
    }

    double diff1 = abs(massZ-mass1);
    double diff2 = abs(massZ-mass2);
    double diff3 = abs(massZ-mass3);

    if(diff1<diff2 && diff1<diff3){
      dileptonTLV = leptonTLVs[0]+leptonTLVs[1];
      oneleptonTLV = leptonTLVs[2];
      lepfromW = 3 ;
      if ( leptonCharges[0]*leptonCharges[1] > 0 ) return true;
    }
    if(diff2<diff1 && diff2<diff3){
      dileptonTLV = leptonTLVs[0]+leptonTLVs[2];
      oneleptonTLV = leptonTLVs[1];
      lepfromW = 2 ;
      if ( leptonCharges[0]*leptonCharges[2] > 0 ) return true;
    }
    if(diff3<diff1 && diff3<diff2){
      dileptonTLV = leptonTLVs[1]+leptonTLVs[2];
      oneleptonTLV = leptonTLVs[0];
      lepfromW = 1 ;
      if ( leptonCharges[1]*leptonCharges[2] > 0 ) return true;
    }
  }
  else{
    dileptonTLV = leptonTLVs[0]+leptonTLVs[1];
    oneleptonTLV = leptonTLVs[2];
    lepfromW = 3 ;
    if ( leptonCharges[0]*leptonCharges[1] > 0 ) return true;
  }

  m_counter->increment(weight,incr++,"2 OSSF Signal Leptons",trueTopo);

  if(leptonCharges.size()!=3) return true;
  m_counter->increment(weight,incr++,"3 Signal Leptons",trueTopo);

  // Lepton isolation #################################

  std::vector<float> vtopoetcone20 ;
  std::vector<float> vptvarcone30 ;
  std::vector<float> vptvarcone20 ;

  float e1=0;
  float p1=0;
  float p2=0;

  if(!m_IsTruth){
    if(m_isElectronChannel && isolated_signal_electrons.size()>=2){
      e1 = (isolated_signal_electrons[0].electron())->isolation(xAOD::Iso::topoetcone20);
      vtopoetcone20.push_back(e1);
      e1 = (isolated_signal_electrons[1].electron())->isolation(xAOD::Iso::topoetcone20);
      vtopoetcone20.push_back(e1);

      p1 = (isolated_signal_electrons[0].electron())->isolation(xAOD::Iso::ptvarcone20);
      vptvarcone20.push_back(p1);
      p1 = (isolated_signal_electrons[1].electron())->isolation(xAOD::Iso::ptvarcone20);
      vptvarcone20.push_back(p1);

      p2 = (isolated_signal_electrons[0].electron())->isolation(xAOD::Iso::ptvarcone30);
      vptvarcone30.push_back(p2);
      p2 = (isolated_signal_electrons[1].electron())->isolation(xAOD::Iso::ptvarcone30);
      vptvarcone30.push_back(p2);

      if(isolated_signal_electrons.size()==3){
	e1 = (isolated_signal_electrons[2].electron())->isolation(xAOD::Iso::topoetcone20);
	vtopoetcone20.push_back(e1);
	p1 = (isolated_signal_electrons[2].electron())->isolation(xAOD::Iso::ptvarcone20);
	vptvarcone20.push_back(p1);
	p2 = (isolated_signal_electrons[2].electron())->isolation(xAOD::Iso::ptvarcone30);
	vptvarcone30.push_back(p2);
      }
      else if(isolated_signal_electrons.size()==2 && isolated_signal_muons.size()==1){
	e1 = (isolated_signal_muons[0].muon())->isolation(xAOD::Iso::topoetcone20);
	vtopoetcone20.push_back(e1);
	p1 = (isolated_signal_muons[0].muon())->isolation(xAOD::Iso::ptvarcone20);
	vptvarcone20.push_back(p1);
	p2 = (isolated_signal_muons[0].muon())->isolation(xAOD::Iso::ptvarcone30);
	vptvarcone30.push_back(p2);
      }
      else if(isolated_baseline_electrons.size()==3 && (m_isVR || m_doFake)){
        e1 = (isolated_baseline_electrons[nW].electron())->isolation(xAOD::Iso::topoetcone20);
        vtopoetcone20.push_back(e1);
        p1 = (isolated_baseline_electrons[nW].electron())->isolation(xAOD::Iso::ptvarcone20);
        vptvarcone20.push_back(p1);
        p2 = (isolated_baseline_electrons[nW].electron())->isolation(xAOD::Iso::ptvarcone30);
        vptvarcone30.push_back(p2);
      }
      else if(isolated_baseline_muons.size()==1 && (m_isVR || m_doFake)){
	e1 = (isolated_baseline_muons[0].muon())->isolation(xAOD::Iso::topoetcone20);
	vtopoetcone20.push_back(e1);
        p1 = (isolated_baseline_muons[0].muon())->isolation(xAOD::Iso::ptvarcone20);
        vptvarcone20.push_back(p1);
        p2 = (isolated_baseline_muons[0].muon())->isolation(xAOD::Iso::ptvarcone30);
        vptvarcone30.push_back(p2);
      }
    } // electron

    else if(m_isMuonChannel && isolated_signal_muons.size()>=2){
      e1 = (isolated_signal_muons[0].muon())->isolation(xAOD::Iso::topoetcone20);
      vtopoetcone20.push_back(e1);
      e1 = (isolated_signal_muons[1].muon())->isolation(xAOD::Iso::topoetcone20);
      vtopoetcone20.push_back(e1);

      p1 = (isolated_signal_muons[0].muon())->isolation(xAOD::Iso::ptvarcone20);
      vptvarcone20.push_back(p1);
      p1 = (isolated_signal_muons[1].muon())->isolation(xAOD::Iso::ptvarcone20);
      vptvarcone20.push_back(p1);

      p2 = (isolated_signal_muons[0].muon())->isolation(xAOD::Iso::ptvarcone30);
      vptvarcone30.push_back(p2);
      p2 = (isolated_signal_muons[1].muon())->isolation(xAOD::Iso::ptvarcone30);
      vptvarcone30.push_back(p2);

      if(isolated_signal_muons.size()==3){
	e1 = (isolated_signal_muons[2].muon())->isolation(xAOD::Iso::topoetcone20);
	vtopoetcone20.push_back(e1);
	p1 = (isolated_signal_muons[2].muon())->isolation(xAOD::Iso::ptvarcone20);
	vptvarcone20.push_back(p1);
	p2 = (isolated_signal_muons[2].muon())->isolation(xAOD::Iso::ptvarcone30);
	vptvarcone30.push_back(p2);
      }
      else if(isolated_signal_muons.size()==2 && isolated_signal_electrons.size()==1){
	e1 = (isolated_signal_electrons[0].electron())->isolation(xAOD::Iso::topoetcone20);
	vtopoetcone20.push_back(e1);
	p1 = (isolated_signal_electrons[0].electron())->isolation(xAOD::Iso::ptvarcone20);
	vptvarcone20.push_back(p1);
	p2 = (isolated_signal_electrons[0].electron())->isolation(xAOD::Iso::ptvarcone30);
	vptvarcone30.push_back(p2);
      }
      else if(isolated_baseline_muons.size()==3 && (m_isVR || m_doFake)){
	e1 = (isolated_baseline_muons[nW].muon())->isolation(xAOD::Iso::topoetcone20);
	vtopoetcone20.push_back(e1);
	p1 = (isolated_baseline_muons[nW].muon())->isolation(xAOD::Iso::ptvarcone20);
	vptvarcone20.push_back(p1);
	p2 = (isolated_baseline_muons[nW].muon())->isolation(xAOD::Iso::ptvarcone30);
	vptvarcone30.push_back(p2);
      }
      else if(isolated_baseline_electrons.size()==1 && (m_isVR || m_doFake)){
	e1 = (isolated_baseline_electrons[0].electron())->isolation(xAOD::Iso::topoetcone20);
	vtopoetcone20.push_back(e1);
	p1 = (isolated_baseline_electrons[0].electron())->isolation(xAOD::Iso::ptvarcone20);
	vptvarcone20.push_back(p1);
	p2 = (isolated_baseline_electrons[0].electron())->isolation(xAOD::Iso::ptvarcone30);
	vptvarcone30.push_back(p2);
      }
    } // muon
  } // not truth

  // Apply Lepton scale factors
  float muSF = eventInfo->auxdecor<float>("muSF");
  if ( muSF != 0.f ) weight *= muSF;
  float elSF = eventInfo->auxdecor<float>("elSF");
  if ( elSF != 0.f ) weight *= elSF;

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
  double Zpt = dileptonTLV.Pt();
  if (InvMassLepPair < 66000 || InvMassLepPair > 116000) return true;
  m_counter->increment(weight,incr++,"M_ll Cut",trueTopo);

  // Add lepton to jets (SR) or MET (VR)
  TVector2 missingETPrime =  *missingET;
  if(lepfromW!=1) missingETPrime = missingETPrime + TVector2(leptonTLVs[0].Px(),leptonTLVs[0].Py()) ;
  if(lepfromW!=2) missingETPrime = missingETPrime + TVector2(leptonTLVs[1].Px(),leptonTLVs[1].Py()) ;
  if(lepfromW!=3) missingETPrime = missingETPrime + TVector2(leptonTLVs[2].Px(),leptonTLVs[2].Py()) ;
  double MissingEtPrime = missingETPrime.Mod();

  // W pt and mT measurements
  double wpx = oneleptonTLV.Px()+missingET->Px();
  double wpy = oneleptonTLV.Py()+missingET->Py();
  double Wpt = std::sqrt(wpx*wpx+wpy*wpy);
  double met = std::sqrt(missingET->Px()*missingET->Px()+missingET->Py()*missingET->Py());
  double mt = std::sqrt(2.*oneleptonTLV.Pt() * met * (1. - (oneleptonTLV.Px()*missingET->Px() + oneleptonTLV.Py()*missingET->Py())/(oneleptonTLV.Pt()*met)));
  //lepptfromW = oneleptonTLV.Pt();

  // bad muons for MET cut: based on non isolated muons
  // FIXME do something special with isbadMETmuon when there are signal muons
  //if ( m_proxyUtils.isbadMETmuon(isolated_baseline_muons, MissingEt, *missingET) ) return true;
  //m_counter->increment(weight,incr++,"IsBadMETMuon",trueTopo);

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


  double mT2=-9;
  double mT2_noISR=-9;
  //if (good_jets.size()>=2) mT2 = m_proxyUtils.MT2(good_jets,missingETPrime);



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

  if(m_doSmallNtuple) {
    unsigned int runnum = RunNumber;
    if ( ! m_IsData  && ! m_IsTruth) runnum = mc_channel_number;

    std::vector<float> jetSmearSystW;

    // other cleaning tests
    unsigned int cleaning = 0;
    unsigned int power2 = 1;

    if(!m_IsTruth){

      // bad jet veto
      if ( !bad_jets.empty() ) cleaning += power2;
      power2 *= 2;//0x2

      // bad muon veto
      for ( size_t i = 0; i < isolated_baseline_muons.size(); i++) {
	if ( isolated_baseline_muons[i].passOVerlapRemoval() &&
	     isolated_baseline_muons[i].isBad() ) {
	  cleaning += power2;
	  break;
	}
      }
      power2 *= 2;//0x4

      // Cosmic muon cut
      if ( m_proxyUtils.CosmicMuon(isolated_baseline_muons) )  cleaning += power2;
      power2 *= 2;//0x8

      // bad Tile cut
      if ( m_proxyUtils.badTileVeto(good_jets,*missingET)) cleaning += power2;
      power2 *= 2;//0x10

      // Negative-cell cleaning cut (no longer used)
      power2 *= 2;//0x20

      // average timing of 2 leading jets
      if (fabs(time[0]) > 5) cleaning += power2;
      power2 *= 2;//0x40

      // FIXME why not in CRWT ?
      //bool chfTileVeto =  m_proxyUtils.chfTileVeto(good_jets);
      //if ( chfTileVeto ) cleaning += 4;

      bool chfVeto = m_proxyUtils.chfVeto(good_jets);
      if ( chfVeto )  cleaning += power2;
      power2 *= 2;//0x80

      bool * failMetCleaning = nullptr;
      if ( !store->retrieve<bool>(failMetCleaning,"failMetCleaning").isSuccess() ) throw std::runtime_error("could not retrieve failMetCleaning");
      if ( *failMetCleaning) cleaning+= power2;

    }
    else power2 *= 128;

    float dPhiBadTile = m_proxyUtils.dPhiBadTile(good_jets,*missingET);

    bool isNCBEvent = false;
    if ( m_IsData ) {
      bool* NCBEventFlag = 0;
      if ( !store->retrieve<bool>(NCBEventFlag,"NCBEventFlag").isSuccess() ) throw std::runtime_error("could not retrieve NCBEventFlag");
      isNCBEvent = *NCBEventFlag;
    }

    m_proxyUtils.FillNTVars(m_ntv, runnum, EventNumber, LumiBlockNumber, veto, weight, normWeight, *pileupWeights, genWeight,ttbarWeightHT,ttbarWeightPt2,ttbarAvgPt,WZweight, btag_weight, ctag_weight, b_jets.size(), c_jets.size(), MissingEtPrime, phi_met, missingET_TST->Mod(), missingET_TST->Phi(), Meff, meffincl, minDphi, RemainingminDPhi, good_jets, good_fat_jets, vD2_fat, visWmedium_fat, trueTopo, cleaning, time[0],jetSmearSystW,0,0.,0.,dPhiBadTile,isNCBEvent,m_IsTruth,baseline_taus,signal_taus);

    if ( systag == "" && !m_IsTruth ) {
      std::vector<float>* p_systweights = 0;
      if ( ! store->retrieve(p_systweights,"event_weights"+m_suffix).isSuccess() ) throw std::runtime_error("Could not retrieve event_weights"+m_suffix);
      m_ntv.systWeights = *p_systweights;

      std::vector<float>* p_btagSystweights = 0;
      if ( ! store->retrieve(p_btagSystweights,"btag_weights"+m_suffix).isSuccess() ) throw std::runtime_error("Could not retrieve btag_weights"+m_suffix);
      m_ntv.btagSystWeights = *p_btagSystweights;
    }

    m_proxyUtils.FillNTExtraVars(m_extrantv, MET_Track, MET_Track_phi, mT2, mT2_noISR, Ap);

    if (  m_fillTRJigsawVars ) m_proxyUtils.FillNTRJigsawVars(m_rjigsawntv, RJigsawVariables );


    std::vector<float> vReclJetMass ;
    std::vector<float> vReclJetPt;
    std::vector<float> vReclJetEta;
    std::vector<float> vReclJetPhi;
    std::vector<bool> visWtight ;
    std::vector<bool> visWmedium ;
    std::vector<bool> visZtight ;
    std::vector<bool> visZmedium ;
    std::vector<float> vD2;

    if( !m_IsTruth && m_fillReclusteringVars){
      m_proxyUtils.FillNTReclusteringVars(m_RTntv,good_jets,vReclJetMass,vReclJetPt,vReclJetEta,vReclJetPhi,vD2,visWmedium, visWtight, visZmedium, visZtight);
    }
    FillCR3LVars(m_cr3lntv, leptonTLVs, leptonCharges, lepfromW, InvMassLepPair, vptvarcone20, vptvarcone30, vtopoetcone20, m_IsTruth, leptonSignal, fake_weight, Zpt, Wpt, mt);
    m_tree->Fill();
  }
  return true;
}

void ZeroLeptonCR3L::finish()
{
  if ( m_DoSystematics ) {
    out() << m_counterRepository << std::endl;
  }
  else {
    out() << *m_counter << std::endl;
  }
}

void ZeroLeptonCR3L::FillCR3LVars(NTCR3LVars& cr3lvars, std::vector<TLorentzVector>& leptons, std::vector<int> lepsigns, int lepfromW, float InvMassLepPair, std::vector<float> vptvarcone20, std::vector<float> vptvarcone30, std::vector<float> vetcone20, bool m_IsTruth, std::vector<int> lepSignal, std::vector<float> fakeWeight, float Zpt, float Wpt, float mt)
{
  cr3lvars.Reset();
  cr3lvars.lep1sign = lepsigns.at(0);
  cr3lvars.lep1Pt  = (leptons.at(0)).Pt() * 0.001;
  cr3lvars.lep1Eta = (leptons.at(0)).Eta();
  cr3lvars.lep1Phi = (leptons.at(0)).Phi();
  cr3lvars.lep1Signal = lepSignal.at(0);
  if(!m_IsTruth){
    cr3lvars.lep1ptvarcone20 = vptvarcone20.at(0)*0.001;
    cr3lvars.lep1ptvarcone30 = vptvarcone30.at(0)*0.001;
    cr3lvars.lep1topoetcone20    = vetcone20.at(0)*0.001;
  }
  cr3lvars.lep2sign = lepsigns.at(1);
  cr3lvars.lep2Pt  = (leptons.at(1)).Pt() * 0.001;
  cr3lvars.lep2Eta = (leptons.at(1)).Eta();
  cr3lvars.lep2Phi = (leptons.at(1)).Phi();
  cr3lvars.lep2Signal = lepSignal.at(1);
  if(!m_IsTruth){
    cr3lvars.lep2ptvarcone20 = vptvarcone20.at(1)*0.001;
    cr3lvars.lep2ptvarcone30 = vptvarcone30.at(1)*0.001;
    cr3lvars.lep2topoetcone20    = vetcone20.at(1)*0.001;
  }
  cr3lvars.lep3sign = lepsigns.at(2);
  cr3lvars.lep3Pt  = (leptons.at(2)).Pt() * 0.001;
  cr3lvars.lep3Eta = (leptons.at(2)).Eta();
  cr3lvars.lep3Phi = (leptons.at(2)).Phi();
  cr3lvars.lep3Signal = lepSignal.at(2);
  if(!m_IsTruth){
    cr3lvars.lep3ptvarcone20 = vptvarcone20.at(2)*0.001;
    cr3lvars.lep3ptvarcone30 = vptvarcone30.at(2)*0.001;
    cr3lvars.lep3topoetcone20    = vetcone20.at(2)*0.001;
  }
  cr3lvars.lepfromW = lepfromW ;
  //cr3lvars.lepptfromW = lepptfromW * 0.001;

  cr3lvars.mll = InvMassLepPair * 0.001;
  cr3lvars.Zpt = Zpt * 0.001;
  cr3lvars.mt  = mt * 0.001;
  cr3lvars.Wpt = Wpt * 0.001;

  if(fakeWeight.size()>0) cr3lvars.fakeWeight = fakeWeight.at(0);
  if(fakeWeight.size()>1) cr3lvars.fakeWeightUp = fakeWeight.at(1);
  if(fakeWeight.size()>2) cr3lvars.fakeWeightDown = fakeWeight.at(2);
}



ClassImp(ZeroLeptonCR3L);

