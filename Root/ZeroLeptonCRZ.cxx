
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
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthParticle.h"


#include <iostream>
#include <stdexcept>
#include <algorithm>

ZeroLeptonCRZ::ZeroLeptonCRZ(const char *name)
  : cafe::Processor(name), 
    m_tree(0), 
    m_stringRegion("CRZ_SRAll"), 
    m_doSmallNtuple(true),
    m_IsData(false),
    m_IsSignal(false),
    m_IsTruth(false),
    m_DoSystematics(false),
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
  m_DoSystematics = config.get("DoSystematics",false);

  m_period = periodFromString(config.get("Period","p13tev"));
  if ( m_period == p7tev ) throw(std::domain_error("ZeroLeptonCRZ does not support the 7tev run period"));
  if ( m_period == INVALID ) throw(std::domain_error("ZeroLeptonCRZ: invalid run period specified"));

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
  if ( !m_DoSystematics && m_counter ) delete m_counter;
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
  bookNTRJigsawVars(tree,m_rjigsawntv);
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

  if ( m_DoSystematics ) {
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
  for ( std::vector<ElectronProxy>::iterator it = isolated_signal_electrons.begin();
	it != isolated_signal_electrons.end(); ) {
    if ( it->Pt() < 25000. ) it = isolated_signal_electrons.erase(it);
    else it++;
  }
  for ( std::vector<ElectronTruthProxy>::iterator itt = isolated_signal_electrons_truth.begin();
        itt != isolated_signal_electrons_truth.end(); ) {
    if ( itt->Pt() < 25000. ) itt = isolated_signal_electrons_truth.erase(itt);
    else itt++;
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
  for ( std::vector<MuonTruthProxy>::iterator itt = isolated_signal_muons_truth.begin();
        itt != isolated_signal_muons_truth.end(); ) {
    if ( itt->Pt() < 25000. ) itt = isolated_signal_muons_truth.erase(itt);
    else itt++;
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
  if ( m_isMuonChannel && 
       ((!m_IsTruth && isolated_signal_muons.size()==2) || (m_IsTruth && isolated_signal_muons_truth.size()==2)) &&
       ((!m_IsTruth && isolated_baseline_electrons.size()==0) || (m_IsTruth && isolated_baseline_electrons_truth.size()==0)) ) {
    if(!m_IsTruth){
      leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_muons[0]))));
      leptonCharges.push_back((int)(isolated_signal_muons[0].muon()->charge())*13);
      leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_muons[1]))));
      leptonCharges.push_back((int)(isolated_signal_muons[1].muon()->charge())*13);
    }
    
    if(m_IsTruth){
      leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_muons_truth[0]))));
      leptonCharges.push_back((int)(isolated_signal_muons_truth[0].muontruth()->charge())*13);
      leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_muons_truth[1]))));
      leptonCharges.push_back((int)(isolated_signal_muons_truth[1].muontruth()->charge())*13);
    } 
    
  } 
  else if (  m_isElectronChannel && 
             ((!m_IsTruth && isolated_signal_electrons.size()==2) || (m_IsTruth && isolated_signal_electrons_truth.size()==2)) && 
	     ((!m_IsTruth && isolated_baseline_muons.size()==0) || (m_IsTruth && isolated_baseline_muons_truth.size()==0)) ) {
    if(!m_IsTruth){
      leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_electrons[0]))));
      leptonCharges.push_back((int)(isolated_signal_electrons[0].electron()->charge())*11);
      leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_electrons[1]))));
      leptonCharges.push_back((int)(isolated_signal_electrons[1].electron()->charge())*11);
    }
    
    if(m_IsTruth){
      leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_electrons_truth[0]))));
      leptonCharges.push_back((int)(isolated_signal_electrons_truth[0].eltruth()->charge())*11);
      leptonTLVs.push_back(*(dynamic_cast<TLorentzVector*>(&(isolated_signal_electrons_truth[1]))));
      leptonCharges.push_back((int)(isolated_signal_electrons_truth[1].eltruth()->charge())*11);
    }
    
  }
  if ( leptonTLVs.empty() ) return true;
  TLorentzVector dileptonTLV = leptonTLVs[0]+leptonTLVs[1];
  if ( leptonCharges[0]*leptonCharges[1] > 0 ) return true;

  /*
  float lepSF[2];
  lepSF[0] = 1.;
  lepSF[1] = 1.;
  if ( !m_IsData ) {
    if(!isolated_signal_electrons.empty()){
      isolated_signal_electrons[0].getSF(lepSF[0]);
      isolated_signal_electrons[1].getSF(lepSF[1]);
    }
    if(!isolated_signal_muons.empty()){
      isolated_signal_muons[0].getSF(lepSF[0]);
      isolated_signal_muons[1].getSF(lepSF[1]);
    }
  }
  weight *= lepSF[0];
  weight *= lepSF[1];
  */
  m_counter->increment(weight,incr++,"2 OS Signal Leptons",trueTopo);

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
  missingETPrime = missingETPrime +
    TVector2(leptonTLVs[0].Px(),leptonTLVs[0].Py())  + 
    TVector2(leptonTLVs[1].Px(),leptonTLVs[1].Py());
  double MissingEtPrime = missingETPrime.Mod();
  

  // bad muons for MET cut: based on non isolated muons
  // FIXME do something special with isbadMETmuon when there are signal muons
  //if ( m_proxyUtils.isbadMETmuon(baseline_muons, MissingEt, *missingET) ) return true;
  //m_counter->increment(weight,incr++,"IsBadMETMuon",trueTopo);


  // Negative-cell cleaning cut
  bool HasNegCell = 0 ; 
  if(! m_IsTruth){
    HasNegCell = m_ZLUtils.NegCellCleaning(event,*missingET);
  }
  //out() << " NegCell " << HasNegCell << std::endl;


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



  m_proxyUtils.RJigsawInit();
  
  std::map<TString,float> RJigsawVariables;

  m_proxyUtils.CalculateRJigsawVariables(good_jets, 
                                missingETPrime.X(),
                                missingETPrime.Y(),
                                RJigsawVariables);


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
      
      // Negative-cell cleaning cut
      if ( HasNegCell )  cleaning += power2;
      power2 *= 2;
      
      // leading jet timing
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

    m_proxyUtils.FillNTVars(m_ntv, runnum, EventNumber, LumiBlockNumber, veto, weight, normWeight, *pileupWeights, genWeight,ttbarWeightHT,ttbarWeightPt2,ttbarAvgPt,WZweight, btag_weight, ctag_weight, b_jets.size(), c_jets.size(), MissingEtPrime, phi_met, Meff, meffincl, minDphi, RemainingminDPhi, good_jets, trueTopo, cleaning, time[0],jetSmearSystW,0,0.,0.,m_IsTruth,baseline_taus,signal_taus);

    m_proxyUtils.FillNTExtraVars(m_extrantv, MET_Track, MET_Track_phi, mT2, mT2_noISR, gaminvRp1, shatR, mdeltaR, cosptR, gamma_R,dphi_BETA_R, dphi_leg1_leg2, costhetaR, dphi_BETA_Rp1_BETA_R, gamma_Rp1, costhetaRp1, Ap);

    m_proxyUtils.FillNTRJigsawVars(m_rjigsawntv, RJigsawVariables );
      

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
  if ( m_DoSystematics ) {
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
  crzvars.lep1Pt  = (leptons.at(0)).Pt() * 0.001;
  crzvars.lep1Eta = (leptons.at(0)).Eta();
  crzvars.lep1Phi = (leptons.at(0)).Phi();

  crzvars.lep2sign = lepsigns.at(1);
  crzvars.lep2Pt  = (leptons.at(1)).Pt() * 0.001;
  crzvars.lep2Eta = (leptons.at(1)).Eta();
  crzvars.lep2Phi = (leptons.at(1)).Phi();


  //double met = std::sqrt(metv.Px()*metv.Px()+metv.Py()*metv.Py());
  crzvars.mll = (leptons.at(0)+leptons.at(1)).M() * 0.001;

  double zpx = leptons.at(0).Px()+leptons.at(1).Py();
  double zpy = leptons.at(0).Py()+leptons.at(1).Py();
  crzvars.Zpt = std::sqrt(zpx*zpx+zpy*zpy) * 0.001;
}



ClassImp(ZeroLeptonCRZ);

