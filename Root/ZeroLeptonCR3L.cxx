
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
    m_DoSystematics(false),
    m_period(INVALID),
    m_isVR(false),
    m_isMuonChannel(false),
    m_isElectronChannel(false),
    m_doRecl(false),
    m_suffix(""),
    m_suffixRecl(""),
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
  m_fillTRJigsawVars = config.get("fillTRJigsawVars",true);
  m_IsData = config.get("IsData",false);
  m_IsSignal = config.get("IsSignal",false);
  m_doRecl = config.get("doRecl",false);
  m_suffixRecl = config.get("suffixRecl","");
  m_IsTruth = config.get("IsTruth",false);
  m_DoSystematics = config.get("DoSystematics",false);

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

  m_derivationTag = derivationTagFromString(config.get("DerivationTag",""));
  if ( m_derivationTag == INVALID_Derivation ) throw(std::domain_error("ZeroLeptonSR: invalid derivation tag specified"));

  m_isVR = config.get("IsVR",false);
  if ( m_isVR ) m_stringRegion = "VR3L_SRAll";

  std::string cutfile = config.get("cutfile","None");
  if ( cutfile == "None" ) throw(std::domain_error("ZeroLeptonCR3L: invalid cut file specified"));
  m_cutVal.ReadCutValues(cutfile);

  m_suffix = config.get("suffix","");
  m_suffixRecl = config.get("suffixRecl","");
  m_physobjsFiller = new PhysObjProxyFiller(20000.f,10000.f,10000.f,m_suffix,m_doRecl,m_suffixRecl);
  m_physobjsFillerTruth = new PhysObjProxyFillerTruth(20000.f,20000.f,10000.f,m_suffix);
  m_proxyUtils = PhysObjProxyUtils(m_IsData);

  m_ZLUtils = ZeroLeptonUtils(m_IsData, m_derivationTag);
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

  if(! m_IsTruth){
    bool passEltrigger=false;
    bool passMutrigger=false;
    if((int)eventInfo->auxdata<char>("HLT_e24_lhmedium_iloose_L1EM18VH")==1 || (int)eventInfo->auxdata<char>("HLT_e60_lhmedium")==1) passEltrigger = true;
    if((int)eventInfo->auxdata<char>("HLT_mu20_iloose_L1MU15")==1 || (int)eventInfo->auxdata<char>("HLT_mu50")==1) passMutrigger = true;
    //if((int)eventInfo->auxdata<char>("HLT_e17_lhloose_L1EM15")==1 || (int)eventInfo->auxdata<char>("HLT_e17_loose_L1EM15")==1) passEltrigger = true;
    //if((int)eventInfo->auxdata<char>("HLT_mu14_iloose")==1 || (int)eventInfo->auxdata<char>("HLT_mu18")==1) passMutrigger = true;
    if( !(passEltrigger || passMutrigger) ) return true;
  }

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
  if(! m_IsTruth){
    if ( ! store->retrieve<TVector2>(missingET,"SUSYMET"+m_suffix+systag).isSuccess() ) throw std::runtime_error("could not retrieve SUSYMET"+m_suffix+systag);
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

  // primary vertex cut
  const xAOD::Vertex* primVertex = 0;
  if(! m_IsTruth){
    primVertex = ZeroLeptonUtils::GetPrimVtx(event);
    if ( !primVertex ||  !( primVertex->nTrackParticles() > 4) ) return true;
  }
  m_counter->increment(weight,incr++,"Vertex Cut",trueTopo);

  // FIXME fake lepton bkg estimate
  // do we need to re-implement that ?

  std::vector<TLorentzVector> leptonTLVs;
  std::vector<int> leptonCharges;

  int nW=1000;

  if ( m_isMuonChannel &&
       ( (!m_isVR && !m_IsTruth && (isolated_signal_muons.size()==3 || (isolated_signal_muons.size()==2 && isolated_signal_electrons.size()==1)) )
	 || ( m_isVR && !m_IsTruth && isolated_signal_muons.size()==2 && isolated_signal_electrons.size()==0 && (baseline_muons.size()==3 || baseline_electrons.size()==1) )
	 || (!m_isVR && m_IsTruth && (isolated_signal_muons_truth.size()==3 || (isolated_signal_muons_truth.size()==2 && isolated_signal_electrons_truth.size()==1)) )
	 || ( m_isVR && m_IsTruth &&  isolated_signal_muons_truth.size()==2 && isolated_signal_electrons_truth.size()==0 && (baseline_muons_truth.size()==3 || baseline_electrons_truth.size()==1) ))){

    if(!m_IsTruth){
      leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_muons[0]))));
      leptonCharges.push_back((int)(isolated_signal_muons[0].muon()->charge())*13);
      leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_muons[1]))));
      leptonCharges.push_back((int)(isolated_signal_muons[1].muon()->charge())*13);

      if(isolated_signal_muons.size()==3){
	leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_muons[2]))));
	leptonCharges.push_back((int)(isolated_signal_muons[2].muon()->charge())*13);
      } // 3 signal muons
      if(isolated_signal_muons.size()==2){
	if(!m_isVR){
	  leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_electrons[0]))));
	  leptonCharges.push_back((int)(isolated_signal_electrons[0].electron()->charge())*11);
	} // not validation region
	if( m_isVR){
	  if(baseline_muons.size()==3){
	    for(int i=0;i<3;i++){
	      if(baseline_muons[i].Pt()==isolated_signal_muons[0].Pt() || baseline_muons[i].Pt()==isolated_signal_muons[1].Pt())
		continue;
	      leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(baseline_muons[i]))));
	      leptonCharges.push_back((int)(baseline_muons[i].muon()->charge())*13);
	      nW = i ;
	    }
	  } // 3 baseline m
	  if(baseline_electrons.size()==1){
	    leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(baseline_electrons[0]))));
	    leptonCharges.push_back((int)(baseline_electrons[0].electron()->charge())*11);
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
	  if(baseline_muons_truth.size()==3){
            for(int i=0;i<3;i++){
              if(baseline_muons_truth[i].Pt()==isolated_signal_muons_truth[0].Pt() || baseline_muons_truth[i].Pt()==isolated_signal_muons_truth[1].Pt())
                continue;
              leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(baseline_muons_truth[i]))));
              leptonCharges.push_back((int)(baseline_muons_truth[i].muontruth()->charge())*13);
	      nW=i;
            }
          } // 3 baseline mu
          if(baseline_electrons_truth.size()==1){
            leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(baseline_electrons_truth[0]))));
            leptonCharges.push_back((int)(baseline_electrons_truth[0].eltruth()->charge())*11);
          } // 1 baseline e

	} // validation region
      } // 2 signal muons
    } // m_IsTruth
  }

  // ELECTRONS

  if ( m_isElectronChannel &&
       (  (!m_isVR && !m_IsTruth && (isolated_signal_electrons.size()==3 || (isolated_signal_electrons.size()==2 && isolated_signal_muons.size()==1)) )
	  || ( m_isVR && !m_IsTruth && isolated_signal_electrons.size()==2 && isolated_signal_muons.size()==0 && (baseline_electrons.size()==3 || baseline_muons.size()==1) )
	  || (!m_isVR && m_IsTruth && (isolated_signal_electrons_truth.size()==3 || (isolated_signal_electrons_truth.size()==2 && isolated_signal_muons_truth.size()==1)) )
	  || ( m_isVR && m_IsTruth &&  isolated_signal_electrons_truth.size()==2 && isolated_signal_muons_truth.size()==0 && (baseline_electrons_truth.size()==3 || baseline_muons_truth.size()==1) ))){

    if(!m_IsTruth){
      leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_electrons[0]))));
      leptonCharges.push_back((int)(isolated_signal_electrons[0].electron()->charge())*11);
      leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_electrons[1]))));
      leptonCharges.push_back((int)(isolated_signal_electrons[1].electron()->charge())*11);

      if(isolated_signal_electrons.size()==3){
        leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_electrons[2]))));
        leptonCharges.push_back((int)(isolated_signal_electrons[2].electron()->charge())*11);
      } // 3 signal e
      if(isolated_signal_electrons.size()==2){
        if(!m_isVR){
          leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_muons[0]))));
          leptonCharges.push_back((int)(isolated_signal_muons[0].muon()->charge())*13);
        } // not validation region
        if( m_isVR){
          if(baseline_electrons.size()==3){
            for(int i=0;i<3;i++){
              if(baseline_electrons[i].Pt()==isolated_signal_electrons[0].Pt() || baseline_electrons[i].Pt()==isolated_signal_electrons[1].Pt())
                continue;
              leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(baseline_electrons[i]))));
              leptonCharges.push_back((int)(baseline_electrons[i].electron()->charge())*11);
	      nW=i;
            }
          } // 3 baseline e
          if(baseline_muons.size()==1){
            leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(baseline_muons[0]))));
            leptonCharges.push_back((int)(baseline_muons[0].muon()->charge())*13);
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
          if(baseline_electrons_truth.size()==3){
            for(int i=0;i<3;i++){
              if(baseline_electrons_truth[i].Pt()==isolated_signal_electrons_truth[0].Pt() || baseline_electrons_truth[i].Pt()==isolated_signal_electrons_truth[1].Pt())
                continue;
              leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(baseline_electrons_truth[i]))));
              leptonCharges.push_back((int)(baseline_electrons_truth[i].eltruth()->charge())*11);
	      nW=i;
            }
          } // 3 baseline electrons
          if(baseline_muons_truth.size()==1){
            leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(baseline_muons_truth[0]))));
            leptonCharges.push_back((int)(baseline_muons_truth[0].muontruth()->charge())*13);
          } // 1 baseline muon
	} // validation region
      } // 2 signal e
    } // isTruth
  } // Electrons


  if ( leptonTLVs.empty() ) return true;

  int lepfromW = 0 ;
  float lepptfromW = 0 ;

  TLorentzVector dileptonTLV ;
  // CASE WITH THREE IDENTICAL LEPTONS - not needed for VR, as we have only two OS signal leptons
  if( !m_isVR && ((!m_IsTruth && (isolated_signal_electrons.size()==3 || isolated_signal_muons.size()==3))
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
      lepfromW = 3 ;
      lepptfromW = leptonTLVs[2].Pt();
      if ( leptonCharges[0]*leptonCharges[1] > 0 ) return true;
    }
    if(diff2<diff1 && diff2<diff3){
      dileptonTLV = leptonTLVs[0]+leptonTLVs[2];
      lepfromW = 2 ;
      lepptfromW = leptonTLVs[1].Pt();
      if ( leptonCharges[0]*leptonCharges[2] > 0 ) return true;
    }
    if(diff3<diff1 && diff3<diff2){
      dileptonTLV = leptonTLVs[1]+leptonTLVs[2];
      lepfromW = 1 ;
      lepptfromW = leptonTLVs[0].Pt();
      if ( leptonCharges[1]*leptonCharges[2] > 0 ) return true;
    }
  }
  else{
    dileptonTLV = leptonTLVs[0]+leptonTLVs[1];
    lepfromW = 3 ;
    lepptfromW = leptonTLVs[2].Pt();
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

      if(isolated_signal_electrons.size()==3 && !m_isVR){
	e1 = (isolated_signal_electrons[2].electron())->isolation(xAOD::Iso::topoetcone20);
	vtopoetcone20.push_back(e1);
	p1 = (isolated_signal_electrons[2].electron())->isolation(xAOD::Iso::ptvarcone20);
	vptvarcone20.push_back(p1);
	p2 = (isolated_signal_electrons[2].electron())->isolation(xAOD::Iso::ptvarcone30);
	vptvarcone30.push_back(p2);
      }
      if(isolated_signal_electrons.size()==2 && isolated_signal_muons.size()==1 && !m_isVR){
	e1 = (isolated_signal_muons[0].muon())->isolation(xAOD::Iso::topoetcone20);
	vtopoetcone20.push_back(e1);
	p1 = (isolated_signal_muons[0].muon())->isolation(xAOD::Iso::ptvarcone20);
	vptvarcone20.push_back(p1);
	p2 = (isolated_signal_muons[0].muon())->isolation(xAOD::Iso::ptvarcone30);
	vptvarcone30.push_back(p2);
      }
      if(baseline_electrons.size()==3 && m_isVR){
        e1 = (baseline_electrons[nW].electron())->isolation(xAOD::Iso::topoetcone20);
        vtopoetcone20.push_back(e1);
        p1 = (baseline_electrons[nW].electron())->isolation(xAOD::Iso::ptvarcone20);
        vptvarcone20.push_back(p1);
        p2 = (baseline_electrons[nW].electron())->isolation(xAOD::Iso::ptvarcone30);
        vptvarcone30.push_back(p2);
      }
      if(baseline_muons.size()==1 && m_isVR){
	e1 = (baseline_muons[0].muon())->isolation(xAOD::Iso::topoetcone20);
	vtopoetcone20.push_back(e1);
        p1 = (baseline_muons[0].muon())->isolation(xAOD::Iso::ptvarcone20);
        vptvarcone20.push_back(p1);
        p2 = (baseline_muons[0].muon())->isolation(xAOD::Iso::ptvarcone30);
        vptvarcone30.push_back(p2);
      }
    } // electron

    if(m_isMuonChannel && isolated_signal_muons.size()>=2){
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

      if(isolated_signal_muons.size()==3 && !m_isVR){
	e1 = (isolated_signal_muons[2].muon())->isolation(xAOD::Iso::topoetcone20);
	vtopoetcone20.push_back(e1);
	p1 = (isolated_signal_muons[2].muon())->isolation(xAOD::Iso::ptvarcone20);
	vptvarcone20.push_back(p1);
	p2 = (isolated_signal_muons[2].muon())->isolation(xAOD::Iso::ptvarcone30);
	vptvarcone30.push_back(p2);
      }
      if(isolated_signal_muons.size()==2 && isolated_signal_electrons.size()==1 && !m_isVR){
	e1 = (isolated_signal_electrons[0].electron())->isolation(xAOD::Iso::topoetcone20);
	vtopoetcone20.push_back(e1);
	p1 = (isolated_signal_electrons[0].electron())->isolation(xAOD::Iso::ptvarcone20);
	vptvarcone20.push_back(p1);
	p2 = (isolated_signal_electrons[0].electron())->isolation(xAOD::Iso::ptvarcone30);
	vptvarcone30.push_back(p2);
      }
      if(baseline_muons.size()==3 && m_isVR){
	e1 = (baseline_muons[nW].muon())->isolation(xAOD::Iso::topoetcone20);
	vtopoetcone20.push_back(e1);
	p1 = (baseline_muons[nW].muon())->isolation(xAOD::Iso::ptvarcone20);
	vptvarcone20.push_back(p1);
	p2 = (baseline_muons[nW].muon())->isolation(xAOD::Iso::ptvarcone30);
	vptvarcone30.push_back(p2);
      }
      if(baseline_electrons.size()==1 && m_isVR){
	e1 = (baseline_electrons[0].electron())->isolation(xAOD::Iso::topoetcone20);
	vtopoetcone20.push_back(e1);
	p1 = (baseline_electrons[0].electron())->isolation(xAOD::Iso::ptvarcone20);
	vptvarcone20.push_back(p1);
	p2 = (baseline_electrons[0].electron())->isolation(xAOD::Iso::ptvarcone30);
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
  if (InvMassLepPair < 66000 || InvMassLepPair > 116000) return true;
  m_counter->increment(weight,incr++,"M_ll Cut",trueTopo);

  // Add lepton to jets (SR) or MET (VR)
  TVector2 missingETPrime =  *missingET;
  if(lepfromW!=1) missingETPrime = missingETPrime + TVector2(leptonTLVs[0].Px(),leptonTLVs[0].Py()) ;
  if(lepfromW!=2) missingETPrime = missingETPrime + TVector2(leptonTLVs[1].Px(),leptonTLVs[1].Py()) ;
  if(lepfromW!=3) missingETPrime = missingETPrime + TVector2(leptonTLVs[2].Px(),leptonTLVs[2].Py()) ;
  double MissingEtPrime = missingETPrime.Mod();


  // bad muons for MET cut: based on non isolated muons
  // FIXME do something special with isbadMETmuon when there are signal muons
  //if ( m_proxyUtils.isbadMETmuon(baseline_muons, MissingEt, *missingET) ) return true;
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
      power2 *= 2;

      // bad muon veto
      for ( size_t i = 0; i < isolated_baseline_muons.size(); i++) {
	if ( isolated_baseline_muons[i].passOVerlapRemoval() &&
	     isolated_baseline_muons[i].isBad() ) {
	  cleaning += power2;
	  break;
	}
      }
      power2 *= 2;

      // Cosmic muon cut
      if ( m_proxyUtils.CosmicMuon(isolated_baseline_muons) )  cleaning += power2;
      power2 *= 2;

      // bad Tile cut
      if ( m_proxyUtils.badTileVeto(good_jets,*missingET)) cleaning += power2;
      power2 *= 2;

      // Negative-cell cleaning cut (no longer used)
      power2 *= 2;

      // average timing of 2 leading jets
      if (fabs(time[0]) > 5) cleaning += power2;
      power2 *= 2;

      // FIXME why not in CRWT ?
      //bool chfTileVeto =  m_proxyUtils.chfTileVeto(good_jets);
      //if ( chfTileVeto ) cleaning += 4;

      bool chfVeto = m_proxyUtils.chfVeto(good_jets);
      if ( chfVeto )  cleaning += power2;
      power2 *= 2;
    }
    else power2 *= 128;

    m_proxyUtils.FillNTVars(m_ntv, runnum, EventNumber, LumiBlockNumber, veto, weight, normWeight, *pileupWeights, genWeight,ttbarWeightHT,ttbarWeightPt2,ttbarAvgPt,WZweight, btag_weight, ctag_weight, b_jets.size(), c_jets.size(), MissingEtPrime, phi_met, Meff, meffincl, minDphi, RemainingminDPhi, good_jets, trueTopo, cleaning, time[0],jetSmearSystW,0,0.,0.,m_IsTruth,baseline_taus,signal_taus );

    if ( systag == "" && !m_IsTruth ) {
      std::vector<float>* p_systweights = 0;
      if ( ! store->retrieve(p_systweights,"event_weights"+m_suffix).isSuccess() ) throw std::runtime_error("Could not retrieve event_weights"+m_suffix);
      m_ntv.systWeights = *p_systweights;
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
    FillCR3LVars(m_cr3lntv, leptonTLVs, *missingET, leptonCharges, lepfromW, InvMassLepPair, vptvarcone20, vptvarcone30, vtopoetcone20, m_IsTruth, lepptfromW);
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

void ZeroLeptonCR3L::FillCR3LVars(NTCR3LVars& cr3lvars, std::vector<TLorentzVector>& leptons, const TVector2& metv, std::vector<int> lepsigns, int lepfromW, float InvMassLepPair, std::vector<float> vptvarcone20, std::vector<float> vptvarcone30, std::vector<float> vetcone20, bool m_IsTruth, float lepptfromW)
{
  cr3lvars.Reset();
  cr3lvars.lep1sign = lepsigns.at(0);
  cr3lvars.lep1Pt  = (leptons.at(0)).Pt() * 0.001;
  cr3lvars.lep1Eta = (leptons.at(0)).Eta();
  cr3lvars.lep1Phi = (leptons.at(0)).Phi();
  if(!m_IsTruth){
    cr3lvars.lep1ptvarcone20 = vptvarcone20.at(0);
    cr3lvars.lep1ptvarcone30 = vptvarcone30.at(0);
    cr3lvars.lep1topoetcone20    = vetcone20.at(0);
  }
  cr3lvars.lep2sign = lepsigns.at(1);
  cr3lvars.lep2Pt  = (leptons.at(1)).Pt() * 0.001;
  cr3lvars.lep2Eta = (leptons.at(1)).Eta();
  cr3lvars.lep2Phi = (leptons.at(1)).Phi();
  if(!m_IsTruth){
    cr3lvars.lep2ptvarcone20 = vptvarcone20.at(1);
    cr3lvars.lep2ptvarcone30 = vptvarcone30.at(1);
    cr3lvars.lep2topoetcone20    = vetcone20.at(1);
  }
  cr3lvars.lep3sign = lepsigns.at(2);
  cr3lvars.lep3Pt  = (leptons.at(2)).Pt() * 0.001;
  cr3lvars.lep3Eta = (leptons.at(2)).Eta();
  cr3lvars.lep3Phi = (leptons.at(2)).Phi();
  if(!m_IsTruth){
    cr3lvars.lep3ptvarcone20 = vptvarcone20.at(2);
    cr3lvars.lep3ptvarcone30 = vptvarcone30.at(2);
    cr3lvars.lep3topoetcone20    = vetcone20.at(2);
  }
  cr3lvars.lepfromW = lepfromW ;
  cr3lvars.lepptfromW = lepptfromW ;

  //double met = std::sqrt(metv.Px()*metv.Px()+metv.Py()*metv.Py());
  cr3lvars.mll = InvMassLepPair ; //(leptons.at(0)+leptons.at(1)).M() * 0.001;

  double zpx = leptons.at(0).Px()+leptons.at(1).Py();
  double zpy = leptons.at(0).Py()+leptons.at(1).Py();
  cr3lvars.Zpt = std::sqrt(zpx*zpx+zpy*zpy) * 0.001;

  double met = std::sqrt(metv.Px()*metv.Px()+metv.Py()*metv.Py());
  cr3lvars.mt = std::sqrt(2.*leptons.at(2).Pt() * met * (1. - (leptons.at(2).Px()*metv.Px() + leptons.at(2).Py()*metv.Py())/(leptons.at(2).Pt()*met))) * 0.001;

  double wpx = leptons.at(2).Px()+metv.Px();
  double wpy = leptons.at(2).Py()+metv.Py();
  cr3lvars.Wpt = std::sqrt(wpx*wpx+wpy*wpy) * 0.001;

}



ClassImp(ZeroLeptonCR3L);

