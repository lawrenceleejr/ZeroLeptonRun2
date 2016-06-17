#include "ZeroLeptonRun2/ZeroLeptonCRY.h"
#include "ZeroLeptonRun2/PhysObjProxies.h"
#include "ZeroLeptonRun2/CleaningHelper.h"

#include "ElectronPhotonSelectorTools/IAsgPhotonIsEMSelector.h"
#include "AsgTools/IAsgTool.h"
#include "PATCore/IAsgSelectionTool.h"
#include "ElectronPhotonSelectorTools/AsgElectronIsEMSelector.h"
#include "ElectronPhotonSelectorTools/AsgPhotonIsEMSelector.h"
#include "ElectronPhotonSelectorTools/AsgElectronLikelihoodTool.h"
#include "ElectronPhotonShowerShapeFudgeTool/ElectronPhotonShowerShapeFudgeTool.h"

#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TActiveStore.h"
#include "xAODRootAccess/TStore.h"
#include "xAODTracking/Vertex.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthParticle.h"
#include "PATInterfaces/SystematicSet.h"
#include "cafe/Processor.h"
#include "cafe/Controller.h"
#include "cafe/Config.h"
#include "TDirectory.h"
#include "TVector2.h"

#include <iostream>
#include <stdexcept>

ZeroLeptonCRY::ZeroLeptonCRY(const char *name)
  : cafe::Processor(name),
    m_tree(0),
    m_stringRegion("CRY_SRAll"),
    m_doSmallNtuple(true),
    m_fillTRJigsawVars(true),
    m_fillReclusteringVars(true),
    m_IsData(false),
    m_IsTruth(false),
    m_IsSignal(false),
    m_DoSystematics(false),
    m_period(INVALID),
    m_suffix(""),
    m_physobjsFiller(0),
    m_physobjsFillerTruth(0),
    m_cutVal(),
    m_proxyUtils(m_IsData),
    m_ZLUtils(m_IsData),
    m_counter(0),
    m_counterRepository("",false,0),
    m_photonSelIsEM(""),
    m_treeRepository()
{
  cafe::Config config(name);
  m_fillTRJigsawVars = config.get("fillTRJigsawVars",true);
  m_IsData = config.get("IsData",false);
  m_IsSignal = config.get("IsSignal",false);
  m_IsTruth = config.get("IsTruth",false);
  m_DoSystematics = config.get("DoSystematics",false);
  m_period = periodFromString(config.get("Period","p13tev"));
  if ( m_period == p7tev ) throw(std::domain_error("ZeroLeptonCRY does not support the 7tev run period"));
  if ( m_period == INVALID ) throw(std::domain_error("ZeroLeptonCRY: invalid run period specified"));

  std::string cutfile = config.get("cutfile","None");
  if ( cutfile == "None" ) throw(std::domain_error("ZeroLeptonCRY: invalid cut file specified"));
  m_cutVal.ReadCutValues(cutfile);


  m_suffix = config.get("suffix","");
  m_suffixSyst = "test";
  m_physobjsFiller = new PhysObjProxyFiller(20000.f,10000.f,10000.f,25000.f,m_suffix,m_doRecl,m_suffixRecl,m_suffixSyst);
  m_physobjsFillerTruth = new PhysObjProxyFillerTruth(20000.f,10000.f,10000.f,25000.f,m_suffix);
  m_proxyUtils = PhysObjProxyUtils(m_IsData);

  m_ZLUtils = ZeroLeptonUtils(m_IsData);

  string m_phId = "Tight";

  if (m_photonSelIsEM.empty()) {
    if ( asg::ToolStore::contains<AsgPhotonIsEMSelector>("PhotonSelIsEM_" + m_phId) ) {
      m_photonSelIsEM = asg::ToolStore::get<AsgPhotonIsEMSelector>("PhotonSelIsEM_" + m_phId);
    } else {
      AsgPhotonIsEMSelector* photonSelIsEM = new AsgPhotonIsEMSelector("PhotonSelIsEM_" + m_phId);
      if (m_phId == "Tight" ) {
	photonSelIsEM->setProperty("ConfigFile", "ElectronPhotonSelectorTools/offline/mc15_20150712/PhotonIsEMTightSelectorCutDefs.conf") ;
        photonSelIsEM->setProperty("isEMMask", static_cast<unsigned int> (egammaPID::PhotonTight) );
      }
      if (m_phId == "Medium" ) {
        photonSelIsEM->setProperty("ConfigFile", "ElectronPhotonSelectorTools/offline/mc15_20150712/PhotonIsEMMediumSelectorCutDefs.conf") ;
        photonSelIsEM->setProperty("isEMMask", static_cast<unsigned int> (egammaPID::PhotonMedium) );
      }
      if (m_phId == "Loose" ) {
        photonSelIsEM->setProperty("ConfigFile", "ElectronPhotonSelectorTools/offline/mc15_20150712/PhotonIsEMLooseSelectorCutDefs.conf") ;
        photonSelIsEM->setProperty("isEMMask", static_cast<unsigned int> (egammaPID::PhotonLoose) );
      }
      photonSelIsEM->initialize() ;
      m_photonSelIsEM = photonSelIsEM;
    }
  }


}

ZeroLeptonCRY::~ZeroLeptonCRY()
{
  if ( !m_DoSystematics && m_counter ) delete m_counter;
  if ( m_physobjsFiller ) delete m_physobjsFiller;
  if ( m_physobjsFillerTruth ) delete m_physobjsFillerTruth;
}

TTree* ZeroLeptonCRY::bookTree(const std::string& treename)
{
  const char* name(treename.c_str());
  TTree* tree = new TTree(name,"ZeroLepton final optimisation");
  tree->SetDirectory(getDirectory());
  bookNTVars(tree,m_ntv,false);
  if ( m_fillReclusteringVars ) bookNTReclusteringVars(tree,m_RTntv);
  bookNTExtraVars(tree,m_extrantv);
  if ( m_fillTRJigsawVars) bookNTRJigsawVars(tree,m_rjigsawntv);
  bookNTCRYVars(tree,m_cryntv);
  return tree;
}

TTree* ZeroLeptonCRY::getTree(const std::string& treename)
{
  std::map<std::string,TTree*>::const_iterator pos = m_treeRepository.find(treename);
  if ( pos == m_treeRepository.end() ) {
    pos = m_treeRepository.insert(std::make_pair(treename,bookTree(treename))).first;
  }
  return pos->second;
}

void ZeroLeptonCRY::begin()
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


bool ZeroLeptonCRY::processEvent(xAOD::TEvent& event)
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
      m_physobjsFiller->setSuffix(m_suffix);
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
  // FIXME : to be implemented for SM MC
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
  if ( ! m_IsData && (m_period == p8tev || m_period == p13tev2015 || m_period == p13tev2016) && !m_IsTruth) {
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
  bool passPhotonTrigger = false;
  if(! m_IsTruth){
    if ( m_IsData && m_period == p13tev2015 ) {
      if( (int)eventInfo->auxdata<char>("HLT_g120_loose")==1   ) passPhotonTrigger = true;
    } else if ( m_IsData && m_period == p13tev2016 ) {
      if( (int)eventInfo->auxdata<char>("HLT_g140_loose")==1   ) passPhotonTrigger = true;
    } else {
      if( (int)eventInfo->auxdata<char>("HLT_g120_loose")==1   ) passPhotonTrigger = true;
      if( (int)eventInfo->auxdata<char>("HLT_g140_loose")==1   ) passPhotonTrigger = true;
    }
  }

  if (passPhotonTrigger == false) return true;

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

  std::vector<float> time;
  m_proxyUtils.EnergyWeightedTime(good_jets,time);

  // isolated_xxx have overlap removed
  std::vector<ElectronProxy> baseline_electrons, isolated_baseline_electrons, isolated_signal_electrons;
  if(! m_IsTruth){
    m_physobjsFiller->FillElectronProxies(baseline_electrons, isolated_baseline_electrons, isolated_signal_electrons);
  }
  std::vector<ElectronTruthProxy> baseline_electrons_truth, isolated_baseline_electrons_truth, isolated_signal_electrons_truth;
  if(m_IsTruth){
    m_physobjsFillerTruth->FillElectronProxies(baseline_electrons_truth, isolated_baseline_electrons_truth, isolated_signal_electrons_truth);
  }

  // isolated_xxx have overlap removed
  std::vector<MuonProxy> baseline_muons, isolated_baseline_muons, isolated_signal_muons;
  if(! m_IsTruth){
    m_physobjsFiller->FillMuonProxies(baseline_muons, isolated_baseline_muons, isolated_signal_muons);
  }
  std::vector<MuonTruthProxy> baseline_muons_truth, isolated_baseline_muons_truth, isolated_signal_muons_truth;
  if(m_IsTruth){
    m_physobjsFillerTruth->FillMuonProxies(baseline_muons_truth, isolated_baseline_muons_truth, isolated_signal_muons_truth);
  }

  // isolated_xxx have overlap removed
  std::vector<PhotonProxy> baseline_photons, isolated_baseline_photons, isolated_signal_photons;
  if(! m_IsTruth){
    m_physobjsFiller->FillPhotonProxies(baseline_photons, isolated_baseline_photons, isolated_signal_photons);
  }
  std::vector<PhotonTruthProxy> baseline_photons_truth, isolated_baseline_photons_truth, isolated_signal_photons_truth;
  if(m_IsTruth){
    m_physobjsFillerTruth->FillPhotonProxies(baseline_photons_truth, isolated_baseline_photons_truth, isolated_signal_photons_truth);
  }

  // keep only leading photon
  int idph=0;
  for ( std::vector<PhotonProxy>::iterator it = isolated_baseline_photons.begin();
        it != isolated_baseline_photons.end(); ) {
    if ( idph>0 ) it = isolated_baseline_photons.erase(it);
    else it++;
    idph++;
  }

  //if(!m_photonSelIsEM.empty()){
  //  isEMvalue = m_photonSelIsEM->IsemValue();
  //}

  //std::vector<PhotonProxy> allphotons;
  //const xAOD::PhotonContainer* phContainer = 0;
  //float leadPhPt = 0.;
  //xAOD::PhotonContainer::const_iterator leadPh ;
  //std::vector<bool> vtight ;
  //std::vector<bool> vloose ;
  //std::vector<float> vtopoetcone20;
  //std::vector<float> vptvarcone20;
  //std::vector<float> vptcone20;
  //std::vector<float> vtopoetcone40;
  //std::vector<float> vptvarcone40;
  //std::vector<float> vptcone40;
  //std::vector<int>   vtruthType;
  //std::vector<int>   vtruthOrigin;
  //std::vector<int> visEMvalue;
  //std::vector<float> vpt;
  //std::vector<float> veta;

  //if(! m_IsTruth){
    //if ( !store->retrieve(phContainer,"SUSYPhotons"+m_suffix).isSuccess() ){
    //  throw std::runtime_error("Could not retrieve PhotonContainer with key SUSYPhotons"+m_suffix);
    //}
    //for ( auto phit = phContainer->begin(); phit != phContainer->end(); phit++){
      ////allphotons.push_back(PhotonProxy(*phit));
      //// Photon isolation /cvmfs/atlas.cern.ch/repo/sw/ASG/AnalysisBase/2.3.14/ElectronIsolationSelection/Root/IsolationSelectionTool.cxx
      //float topoetcone20=0;
      //float ptvarcone20=0;
      //float ptcone20=0;
      //float topoetcone40=0;
      //float ptvarcone40=0;
      //float ptcone40=0;
      //
      //(*phit)->isolationValue(topoetcone20,xAOD::Iso::topoetcone20);
      //(*phit)->isolationValue(ptvarcone20,xAOD::Iso::ptvarcone20);
      //(*phit)->isolationValue(ptcone20,xAOD::Iso::ptcone20);
      //(*phit)->isolationValue(topoetcone40,xAOD::Iso::topoetcone40);
      //(*phit)->isolationValue(ptvarcone40,xAOD::Iso::ptvarcone40);
      //(*phit)->isolationValue(ptcone40,xAOD::Iso::ptcone40);
      //
      //vptvarcone20.push_back(ptvarcone20);
      //vtopoetcone20.push_back(topoetcone20);
      //vptcone20.push_back(ptcone20);
      //vptvarcone40.push_back(ptvarcone40);
      //vtopoetcone40.push_back(topoetcone40);
      //vptcone40.push_back(ptcone40);

      //int isEMvalue = -1000;
      //if(!m_photonSelIsEM.empty()){
      //  isEMvalue = m_photonSelIsEM->IsemValue();
      //}
      //visEMvalue.push_back(isEMvalue);

      //if(!m_IsData){
      //  vtruthType.push_back((*phit)->auxdata<int>("truthType"));
      //  vtruthOrigin.push_back((*phit)->auxdata<int>("truthOrigin"));
      //}
      //else{
      //  vtruthType.push_back(-1000);
      //  vtruthOrigin.push_back(-1000);
      //}


      //
      //vpt.push_back((*phit)->pt());
      //veta.push_back((*phit)->eta());
      //
      // Photon quality
      //bool tight=false; // tight=1 passes the selection
      //if(!(*phit)->passSelection(tight,"Tight")){
      //	std::cout<<"WARNING: Tight decision is not available"<<std::endl;
      //}
      //bool loose=false;
      //if(!(*phit)->passSelection(loose,"Loose")){
      //	std::cout<<"WARNING: Loose decision is not available"<<std::endl;
      //}
      //vloose.push_back(loose);
      //vtight.push_back(tight);
      ////
      //// https://svnweb.cern.ch/trac/atlasoff/browser/PhysicsAnalysis/ElectronPhotonID/ElectronPhotonSelectorTools/trunk/ElectronPhotonSelectorTools/egammaPIDdefs.h
      ////unsigned int ph_isEMTight = (*phit)->auxdata<unsigned int>("isEMTight");
      ////visEMTight.push_back(ph_isEMTight);
      ////
      //if ( (*phit)->auxdecor<char>("baseline")==1 && (*phit)->pt() > leadPhPt ) {
      //  leadPhPt = (*phit)->pt();
      //  leadPh = phit;
      //}
  //}
  //}
  //
  //const xAOD::TruthParticleContainer* truthphotons = 0 ;
  //xAOD::TruthParticleContainer::const_iterator leadPhtruth;
  //if( m_IsTruth ){
  //  if ( !event.retrieve(truthphotons,"TruthPhotons").isSuccess() ){
  //    throw std::runtime_error("Could not retrieve truth particles with key TruthPhotons");
  //  }
  //  for ( auto phittruth = truthphotons->begin(); phittruth != truthphotons->end(); phittruth++){
  //    allphotons.push_back(PhotonProxy((*phittruth)->p4()));
  //
  //    if ( (*phittruth)->pt() > leadPhPt ) { // GERALDINE - ADD BASELINE DEFINITION
  //	leadPhPt = (*phittruth)->pt();
  //	leadPhtruth = phittruth;
  //    }
  //  }
  //}

  //if ( leadPhPt == 0. ) return true;
  //m_counter->increment(weight,incr++,"One photon",trueTopo);

  //out() << "Leading photon pt px py  " << leadPhPt ;
  //if ( leadPhPt > 0. ) out() << " " << (*leadPh)->p4().Px()<< " " << (*leadPh)->p4().Py();
  //out() <<std::endl;

  std::vector<TauProxy> baseline_taus, signal_taus;
  if(! m_IsTruth){
    m_physobjsFiller->FillTauProxies(baseline_taus, signal_taus);
  }

  // missing ET
  TVector2* missingET = 0;
  TVector2* missingET_TST = 0;
  TVector2 missingETCorr;
  if(! m_IsTruth){
    if ( ! store->retrieve<TVector2>(missingET,"SUSYMET"+m_suffix+systag).isSuccess() ) throw std::runtime_error("could not retrieve SUSYMET"+m_suffix+systag);
    if ( ! store->retrieve<TVector2>(missingET_TST,"SUSYMETTST"+m_suffix+systag).isSuccess() ) throw std::runtime_error("could not retrieve SUSYMET"+m_suffix+systag);
  }
  else {
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

  // at least one photon
  //if ( leadPhPt < m_cutVal.m_cutPhotonPtCRY ) return true;
  bool oneGamma = false;
  if( !isolated_baseline_photons.empty() && ! m_IsTruth )
    oneGamma = true;
  if( !isolated_baseline_photons_truth.empty() && m_IsTruth )
    oneGamma = true;
  if ( !oneGamma ) return true;
  m_counter->increment(weight,incr++,"1 photon",trueTopo);


  TLorentzVector photonTLV;
  if(! m_IsTruth) photonTLV.SetPtEtaPhiM(isolated_baseline_photons[0].photon()->pt(),isolated_baseline_photons[0].photon()->eta(),isolated_baseline_photons[0].photon()->phi(),0);
  if(  m_IsTruth) photonTLV = *(dynamic_cast<TLorentzVector*>(&(isolated_baseline_photons_truth[0])));

  int phSignal = 1;
  if(! m_IsTruth) phSignal = (isolated_baseline_photons[0].isSignal());
  int isEMvalue = -1000;


  //if(m_IsTruth) {
  missingETCorr.Set(missingET->Px()+photonTLV.Px(), missingET->Py()+photonTLV.Py());
  //}
  //else {
  //missingETCorr.Set(missingET->Px()+photonTLV.Px(), missingET->Py()+photonTLV.Py());
  //}
  double MissingEtCorr = missingETCorr.Mod();


  // photon isolation

  float topoetcone20 = 0 ;
  float ptvarcone20 = 0 ;
  float ptcone20 = 0 ;
  float topoetcone40 = 0 ;
  float ptvarcone40 = 0 ;
  float ptcone40 = 0 ;
  bool btight = false;
  bool bloose = false;
  int tight = 0 ;
  int loose = 0 ;
  int truthType = 0 ;
  int truthOrigin = 0 ;

  if(!m_IsTruth && isolated_baseline_photons.size()==1){
    topoetcone20    = (isolated_baseline_photons[0].photon())->isolation(xAOD::Iso::topoetcone20);
    ptvarcone20     = (isolated_baseline_photons[0].photon())->isolation(xAOD::Iso::ptvarcone20);
    ptcone20        = (isolated_baseline_photons[0].photon())->isolation(xAOD::Iso::ptcone20);
    topoetcone40    = (isolated_baseline_photons[0].photon())->isolation(xAOD::Iso::topoetcone40);
    ptvarcone40     = (isolated_baseline_photons[0].photon())->isolation(xAOD::Iso::ptvarcone40);
    ptcone40        = (isolated_baseline_photons[0].photon())->isolation(xAOD::Iso::ptcone40);

    if(!(isolated_baseline_photons[0].photon())->passSelection(btight,"Tight")){
      std::cout<<"WARNING: Tight decision is not available"<<std::endl;
    }
    if(!(isolated_baseline_photons[0].photon())->passSelection(bloose,"Loose")){
      std::cout<<"WARNING: Loose decision is not available"<<std::endl;
    }
    if(bloose) loose = 1;
    if(btight) tight = 1;

  }

  if(!m_IsData){
    if(!m_IsTruth){
      truthType   = (isolated_baseline_photons[0].photon())->auxdata<int>("truthType");
      truthOrigin = (isolated_baseline_photons[0].photon())->auxdata<int>("truthOrigin");
    }
    else{
      static const SG::AuxElement::ConstAccessor<unsigned int> acc_truthType("classifierParticleType");
      static const SG::AuxElement::ConstAccessor<unsigned int> acc_truthOrig("classifierParticleOrigin");
      truthType   = acc_truthType(*(isolated_signal_photons_truth[0].phtruth()));
      truthOrigin = acc_truthOrig(*(isolated_signal_photons_truth[0].phtruth()));
    }
  }else{
    truthType   = -1000;
    truthOrigin = -1000;
  }

  // Apply photon scale factor
  float phSF = 1.;
  if(! m_IsTruth) phSF = isolated_baseline_photons[0].photon()->auxdecor<float>("sf");
  if ( phSF != 0.f ) weight *= phSF;


  // MissingET cut
  if (!(MissingEtCorr > m_cutVal.m_cutEtMiss)) return true;
  m_counter->increment(weight,incr++,"MET cut",trueTopo);

  // Leading jet Pt cut
  if ( good_jets.empty() ) return true;
  if (!(good_jets[0].Pt() > m_cutVal.m_cutJetPt0)) return true;
  m_counter->increment(weight,incr++,"1 jet Pt > 130 GeV Selection",trueTopo);


  // Calculate variables for ntuple -----------------------------------------
  double phi_met = TMath::ATan2(missingETCorr.Y(),missingETCorr.X());
  double minDphi = m_proxyUtils.SmallestdPhi(good_jets,phi_met);
  double RemainingminDPhi = m_proxyUtils.SmallestRemainingdPhi(good_jets,phi_met);
  double Meff[6];

  for ( size_t i = 0; i<6; ++i ) {
    Meff[i] = m_proxyUtils.Meff(good_jets,
				std::max<size_t>(1,i+1),
				MissingEtCorr,
				m_cutVal.m_cutJetPt1,
				m_cutVal.m_cutJetPt4);
  }
  double meffincl = m_proxyUtils.Meff(good_jets,
				      good_jets.size(),
				      MissingEtCorr,
				      m_cutVal.m_cutJetPt4,
				      m_cutVal.m_cutJetPt4);

  // ttbar reweighting not available yet in SUSYOBJDef_xAOD
  float ttbarWeightHT = 1.;
  float ttbarWeightPt2 = 1.;
  float ttbarAvgPt = 0.;
  //float HT = meffincl -  MissingEtCorr;

  // Sherpa MassiveCB W/Z reweighting : not implemented yet in SUSYOBJDef_xAOD
  float WZweight = 1.;
  if( !m_IsData ){
    WZweight = eventInfo->auxdecor<float>("WZweight");
  }

  /* not used now
  double mT2=-9;
  if (good_jets.size()>=2) mT2 = m_proxyUtils.MT2(good_jets,missingETCorr);
  double mT2_noISR=-9;
  if (nonISR_jets.size()>=2) mT2_noISR = m_proxyUtils.MT2(nonISR_jets,missingETCorr);
  out() << " mT2 " << mT2 << " " << mT2_noISR << std::endl;
  */

  std::map<TString,float> RJigsawVariables;
  if ( m_fillTRJigsawVars) {
    m_proxyUtils.CalculateRJigsawVariables(good_jets,
					   missingETCorr.X(),
					   missingETCorr.Y(),
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
			      missingETCorr.X(),
			      missingETCorr.Y(),
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
    if ( ! m_IsData && ! m_IsTruth) runnum = mc_channel_number;

    std::vector<float> jetSmearSystW;

    if(!m_IsTruth){

      // bad jet veto
      if ( !bad_jets.empty() ) m_cleaningHelper.cleaning.at("badJetVeto") = true;

      // bad muon veto
      for ( size_t i = 0; i < isolated_baseline_muons.size(); i++) {
	if ( isolated_baseline_muons[i].passOVerlapRemoval() &&
	     isolated_baseline_muons[i].isBad() ) {
	  m_cleaningHelper.cleaning.at("badMuonVeto") = true;
	  break;
	}
      }

      // Cosmic muon cut
      if ( m_proxyUtils.CosmicMuon(isolated_baseline_muons) )  m_cleaningHelper.cleaning.at("cosmicMuonVeto") = true;

      // bad muons for MET cut: based on non isolated muons
      if ( m_proxyUtils.isbadMETmuon(baseline_muons, MissingEt, *missingET) )  m_cleaningHelper.cleaning.at("badMetMuonVeto") = true;

      // bad Tile cut
      if ( m_proxyUtils.badTileVeto(good_jets,*missingET)) m_cleaningHelper.cleaning.at("badTileVeto") = true;

      // average timing of 2 leading jets
      if (fabs(time[0]) > 5) m_cleaningHelper.cleaning.at("leadingJetTimingVeto") = true;

      if ( m_proxyUtils.chfVeto(good_jets) ) m_cleaningHelper.cleaning.at("chfVeto") = true;

      bool * failMetCleaning = nullptr;
      if ( !store->retrieve<bool>(failMetCleaning,"failMetCleaning").isSuccess() ) throw std::runtime_error("could not retrieve failMetCleaning");
      if ( *failMetCleaning) m_cleaningHelper.cleaning.at("metTSTCleaningVeto") = true;

    }

    unsigned long const cleaning = m_cleaningHelper.finalCleaning();
    float dPhiBadTile = m_proxyUtils.dPhiBadTile(good_jets,*missingET);

    bool isNCBEvent = false;
    if ( m_IsData ) {
      bool* NCBEventFlag = 0;
      if ( !store->retrieve<bool>(NCBEventFlag,"NCBEventFlag").isSuccess() ) throw std::runtime_error("could not retrieve NCBEventFlag");
      isNCBEvent = *NCBEventFlag;
    }

    m_proxyUtils.FillNTVars(m_ntv, runnum, EventNumber, LumiBlockNumber, veto, weight, normWeight, *pileupWeights, PRWHash, genWeight,ttbarWeightHT,ttbarWeightPt2,ttbarAvgPt,WZweight, b_jets.size(), c_jets.size(), MissingEtCorr, phi_met, missingET_TST->Mod(), missingET_TST->Phi(), Meff, meffincl, minDphi, RemainingminDPhi, good_jets, good_fat_jets, vD2_fat, visWmedium_fat, trueTopo, cleaning, time[0],jetSmearSystW,0, 0., dPhiBadTile,isNCBEvent,m_IsTruth,baseline_taus,signal_taus);

    if ( systag == ""  && !m_IsTruth ) {
      std::vector<float>* p_systweights = 0;
      if ( ! store->retrieve(p_systweights,"event_weights"+m_suffix).isSuccess() ) throw std::runtime_error("Could not retrieve event_weights"+m_suffix);
      m_ntv.systWeights = *p_systweights;

      std::vector<float>* p_btagSystweights = 0;
      if ( ! store->retrieve(p_btagSystweights,"btag_weights"+m_suffix).isSuccess() ) throw std::runtime_error("Could not retrieve btag_weights"+m_suffix);
      m_ntv.btagSystWeights = *p_btagSystweights;
    }


    m_proxyUtils.FillNTExtraVars(m_extrantv, MET_Track, MET_Track_phi, Ap);

    if ( m_fillTRJigsawVars) m_proxyUtils.FillNTRJigsawVars(m_rjigsawntv, RJigsawVariables );

    FillNTCRYVars(m_cryntv,photonTLV,*missingET,tight,loose,topoetcone20,ptvarcone20,
		  ptcone20, topoetcone40,ptvarcone40,ptcone40,
		  truthType,truthOrigin,isEMvalue,phSignal);

    std::vector<float> vReclJetMass ;
    std::vector<float> vReclJetPt;
    std::vector<float> vReclJetEta;
    std::vector<float> vReclJetPhi;
    std::vector<bool> visWtight ;
    std::vector<bool> visWmedium ;
    std::vector<bool> visZtight ;
    std::vector<bool> visZmedium ;
    std::vector<float> vD2;

    if(! m_IsTruth && m_fillReclusteringVars){
      m_proxyUtils.FillNTReclusteringVars(m_RTntv,good_jets,vReclJetMass,vReclJetPt,vReclJetEta,vReclJetPhi,vD2,visWmedium, visWtight, visZmedium, visZtight);
    }

    m_tree->Fill();
  }
  return true;
}

void ZeroLeptonCRY::finish()
{
  if ( m_DoSystematics ) {
    out() << m_counterRepository << std::endl;
  }
  else {
    out() << *m_counter << std::endl;
  }
}

void ZeroLeptonCRY::FillNTCRYVars(NTCRYVars& cryntv,
				  const TLorentzVector&  photon,
				  TVector2& origmisset, int tight, int loose,
				  float topoetcone20, float ptvarcone20,float ptcone20,
				  float topoetcone40, float ptvarcone40, float ptcone40,
				  int truthType, int truthOrigin,
				  int isEMvalue, int phSignal)
{
  cryntv.Reset();

  cryntv.phLoose        = loose ;
  cryntv.phTight        = tight ;
  cryntv.phTopoetcone20 = topoetcone20 * 0.001;
  cryntv.phPtvarcone20  = ptvarcone20 * 0.001;
  cryntv.phPtcone20     = ptcone20 * 0.001;
  cryntv.phTopoetcone40 = topoetcone40 * 0.001;
  cryntv.phPtvarcone40  = ptvarcone40 * 0.001 ;
  cryntv.phPtcone40     = ptcone40 * 0.001;
  cryntv.phTruthType    = truthType ;
  cryntv.phTruthOrigin  = truthOrigin ;
  cryntv.phisEMvalue    = isEMvalue ;
  cryntv.phSignal       = phSignal ;

  cryntv.origmet        = origmisset.Mod() * 0.001;
  cryntv.origmetPhi     = origmisset.Phi();

  cryntv.phPt           = photon.Pt() * 0.001;
  cryntv.phEta          = photon.Eta();
  cryntv.phPhi          = photon.Phi();

}

ClassImp(ZeroLeptonCRY);

