#include "ZeroLeptonRun2/ZeroLeptonCRY.h"
#include "ZeroLeptonRun2/PhysObjProxies.h"


#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TActiveStore.h"
#include "xAODRootAccess/TStore.h"
#include "xAODTracking/Vertex.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODEgamma/PhotonContainer.h"
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
    m_IsData(false),
    m_IsSignal(false),
    m_UseSystematics(false),
    m_period(INVALID),
    m_suffix(""),
    m_physobjsFiller(0),
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
  m_UseSystematics = config.get("UseSystematics",false);
  m_period = periodFromString(config.get("Period","p13tev"));
  if ( m_period == p7tev ) throw(std::domain_error("ZeroLeptonCRY does not support the 7tev run period"));
  if ( m_period == INVALID ) throw(std::domain_error("ZeroLeptonCRY: invalid run period specified"));

  m_derivationTag = derivationTagFromString(config.get("DerivationTag",""));
  if ( m_derivationTag == INVALID_Derivation ) throw(std::domain_error("ZeroLeptonSR: invalid derivation tag specified"));

  std::string cutfile = config.get("cutfile","None");
  if ( cutfile == "None" ) throw(std::domain_error("ZeroLeptonCRY: invalid cut file specified"));
  m_cutVal.ReadCutValues(cutfile);


  m_suffix = config.get("suffix","");
  m_physobjsFiller = new PhysObjProxyFiller(20000.f,10000.f,10000.f,m_suffix);
  m_proxyUtils = PhysObjProxyUtils(m_IsData);
  m_ZLUtils = ZeroLeptonUtils(m_IsData, m_derivationTag);
}

ZeroLeptonCRY::~ZeroLeptonCRY()
{
  if ( !m_UseSystematics && m_counter ) delete m_counter;
  if ( m_physobjsFiller ) delete m_physobjsFiller;
}

TTree* ZeroLeptonCRY::bookTree(const std::string& treename)
{
  const char* name(treename.c_str());
  TTree* tree = new TTree(name,"ZeroLepton final optimisation");
  tree->SetDirectory(getDirectory());
  bookNTVars(tree,m_ntv,false);
  bookNTReclusteringVars(tree,m_RTntv);
  bookNTExtraVars(tree,m_extrantv);
  //bookNTCRYVars(tree,m_cryntv);
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

  if ( m_UseSystematics ) {
    m_counterRepository = CounterRepository("ZeroLeptonCounter"+m_stringRegion,m_IsSignal,getDirectory());
  }
  else {
    m_counter = new Counter("ZeroLeptonCounter"+m_stringRegion,40,m_IsSignal);
    if(m_doSmallNtuple) m_tree = bookTree(sSR);
  }

}


bool ZeroLeptonCRY::processEvent(xAOD::TEvent& event)
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
  if ( !m_IsData ) {
    if ( !store->retrieve< std::vector<float> >(pileupWeights,"pileupWeights").isSuccess() ) throw std::runtime_error("could not retrieve pileupWeights");
    //out() << " pileup weight " << (*pileupWeights)[0] << std::endl;
    //weight *= (*pileupWeights)[0];
  }
  else {
    static std::vector<float> dummy(3,0.);
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
  if ( ! m_IsData && (m_period == p8tev || m_period == p13tev)) {
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
  m_physobjsFiller->FillJetProxies(good_jets,bad_jets,b_jets);
  std::vector<float> btag_weight(7,1.); // not implemented in SUSYTools
  std::vector<float> ctag_weight(7,1.); // not implemented in SUSYTools
  std::vector<float> time;
  m_proxyUtils.EnergyWeightedTime(good_jets,time);

  // isolated_xxx have overlap removed
  std::vector<ElectronProxy> baseline_electrons, isolated_baseline_electrons, isolated_signal_electrons;
  m_physobjsFiller->FillElectronProxies(baseline_electrons, isolated_baseline_electrons, isolated_signal_electrons);

  // isolated_xxx have overlap removed
  std::vector<MuonProxy> baseline_muons, isolated_baseline_muons, isolated_signal_muons;
  m_physobjsFiller->FillMuonProxies(baseline_muons, isolated_baseline_muons, isolated_signal_muons);

  const xAOD::PhotonContainer* phContainer = 0;
  if ( !store->retrieve(phContainer,"SUSYPhotons"+m_suffix).isSuccess() ){
    throw std::runtime_error("Could not retrieve PhotonContainer with key SUSYPhotons"+m_suffix);
  }
  xAOD::PhotonContainer::const_iterator leadPh = phContainer->end();
  float leadPhPt = 0.;
  for ( auto phit = phContainer->begin(); phit != phContainer->end(); phit++){
    if ( (*phit)->auxdecor<char>("baseline")==1 && (*phit)->pt() > leadPhPt ) {
      leadPhPt = (*phit)->pt();
      leadPh = phit;
    }
  }
  //out() << "Leading photon pt px py  " << leadPhPt ;
  //if ( leadPhPt > 0. ) out() << " " << (*leadPh)->p4().Px()<< " " << (*leadPh)->p4().Py();
  //out() <<std::endl;  

  // missing ET
  TVector2* missingET = new TVector2(0.,0.);
  TVector2* missingETCorr = missingET;
  if ( ! store->retrieve<TVector2>(missingET,"SUSYMET"+m_suffix).isSuccess() ) throw std::runtime_error("could not retrieve SUSYMET"+m_suffix);
  // the lead photon as if it was a Z->nunu
  if ( leadPhPt > 0. ) {
    //out() << "MET was " << missingET->Px() << " " << missingET->Py() << std::endl;
    missingETCorr->Set(missingET->Px()+(*leadPh)->p4().Px(), missingET->Py()+(*leadPh)->p4().Py());
    //out() << "MET  is " << missingET->Px() << " " << missingET->Py() << std::endl;
  } 
  double MissingEt = missingET->Mod();
  double MissingEtCorr = missingETCorr->Mod();


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
  m_ZLUtils.trackMET(event, MET_Track, MET_Track_phi);

  // primary vertex cut
  const xAOD::Vertex* primVertex = ZeroLeptonUtils::GetPrimVtx(event);
  if ( !primVertex ||  !( primVertex->nTrackParticles() > 4) ) return true;
  m_counter->increment(weight,incr++,"Vertex Cut",trueTopo);

  // Cosmic muon cut
  if ( m_proxyUtils.CosmicMuon(isolated_baseline_muons) ) return true;
  m_counter->increment(weight,incr++,"CosmicMuons",trueTopo);

  //FIXME isBadMuon not yet implement in SUSYObjDef_xAOD

  // bad muons for MET cut: based on non isolated muons
  if ( m_proxyUtils.isbadMETmuon(baseline_muons, MissingEt, *missingET) ) return true;
  m_counter->increment(weight,incr++,"IsBadMETMuon",trueTopo);

  // bad Tile cut
  if ( m_proxyUtils.badTileVeto(good_jets,*missingET)) return true;
  m_counter->increment(weight,incr++,"Bad Tile Veto",trueTopo);

  // Negative-cell cleaning cut
  bool HasNegCell = m_ZLUtils.NegCellCleaning(event,*missingET);
  //out() << " NegCell " << HasNegCell << std::endl;
  if ( HasNegCell ) return true;
  m_counter->increment(weight,incr++,"Negative-cell cleaning",trueTopo);


  // at least one photon
  if ( leadPhPt < m_cutVal.m_cutPhotonPtCRY ) return true;
  m_counter->increment(weight,incr++,">=1 photon",trueTopo);

 
  // MissingET cut
  if (!(MissingEtCorr > m_cutVal.m_cutEtMiss)) return true;
  m_counter->increment(weight,incr++,"MET cut",trueTopo);

  // Leading jet Pt cut
  if ( good_jets.empty() ) return true;
  if (!(good_jets[0].Pt() > m_cutVal.m_cutJetPt0)) return true; 
  m_counter->increment(weight,incr++,"1 jet Pt > 130 GeV Selection",trueTopo);



  // Calculate variables for ntuple -----------------------------------------
  double phi_met = TMath::ATan2(missingETCorr->Y(),missingETCorr->X());
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
  float HT = meffincl -  MissingEtCorr;

  // Sherpa MassiveCB W/Z reweighting : not implemented yet in SUSYOBJDef_xAOD
  float WZweight = 1.;


  double mT2=-9; 
  //if (good_jets.size()>=2) mT2 = m_proxyUtils.MT2(good_jets,*missingETCorr);
  double mT2_noISR=-9;
  //if (nonISR_jets.size()>=2) mT2_noISR = m_proxyUtils.MT2(nonISR_jets,*missingETCorr); 
  //out() << " mT2 " << mT2 << " " << mT2_noISR << std::endl; 


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
			      missingETCorr->X(),
			      missingETCorr->Y(),
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

    bool chfTileVeto =  m_proxyUtils.chfTileVeto(good_jets);
    if ( chfTileVeto ) cleaning += 4;

    bool chfVeto = m_proxyUtils.chfVeto(good_jets);
    if ( chfVeto ) cleaning += 8;

    m_proxyUtils.FillNTVars(m_ntv, runnum, EventNumber, veto, weight, normWeight, *pileupWeights, genWeight,ttbarWeightHT,ttbarWeightPt2,ttbarAvgPt,WZweight, btag_weight, ctag_weight, b_jets.size(), c_jets.size(), MissingEtCorr, phi_met, Meff, meffincl, minDphi, RemainingminDPhi, good_jets, trueTopo, cleaning, time[0],jetSmearSystW,0, 0., 0., 0., 0.);


    m_proxyUtils.FillNTExtraVars(m_extrantv, mT2,mT2_noISR,gaminvRp1 ,shatR ,mdeltaR ,cosptR ,gamma_R,dphi_BETA_R , dphi_leg1_leg2 , costhetaR ,dphi_BETA_Rp1_BETA_R,gamma_Rp1,costhetaRp1,Ap);
      
    m_proxyUtils.FillNTReclusteringVars(m_RTntv,good_jets);

    m_tree->Fill();
  }
  return true;
}

void ZeroLeptonCRY::finish()
{
  if ( m_UseSystematics ) {
    out() << m_counterRepository << std::endl;
  } 
  else {
    out() << *m_counter << std::endl;
  }
}

ClassImp(ZeroLeptonCRY);

