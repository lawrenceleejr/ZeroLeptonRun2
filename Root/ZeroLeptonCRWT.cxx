
#include "ZeroLeptonRun2/ZeroLeptonCRWT.h"
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

ZeroLeptonCRWT::ZeroLeptonCRWT(const char *name)
  : cafe::Processor(name), 
    m_tree(0), 
    m_stringRegion("CRWT_SRAll"), 
    m_doSmallNtuple(true),
    m_IsData(false),
    m_IsSignal(false),
    m_UseSystematics(false),
    m_period(INVALID),
    m_isVR(false),
    m_isMuonChannel(false),
    m_isElectronChannel(false),
    m_suffix(""),
    m_physobjsFiller(0),
    m_cutVal(),
    m_proxyUtils(m_IsData),
    m_ZLUtils(m_IsData),
    m_counter(0),
    m_counterRepository("",false,0),
    m_treeRepository()
{
  cafe::Config config(name);
  m_IsData = config.get("IsData",false);
  m_IsSignal = config.get("IsSignal",false);
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
  if ( m_period == p7tev ) throw(std::domain_error("ZeroLeptonCRWT does not support the 7tev run period"));
  if ( m_period == INVALID ) throw(std::domain_error("ZeroLeptonCRWT: invalid run period specified"));

  m_isVR = config.get("IsVR",false);
  if ( m_isVR ) m_stringRegion = "VRWT_SRAll";

  std::string cutfile = config.get("cutfile","None");
  if ( cutfile == "None" ) throw(std::domain_error("ZeroLeptonCRWT: invalid cut file specified"));
  m_cutVal.ReadCutValues(cutfile);

  m_suffix = config.get("suffix","");
  m_physobjsFiller = new PhysObjProxyFiller(20000.f,10000.f,10000.f,m_suffix);
  m_proxyUtils = PhysObjProxyUtils(m_IsData);
  m_ZLUtils = ZeroLeptonUtils(m_IsData);
}

ZeroLeptonCRWT::~ZeroLeptonCRWT()
{
  if ( !m_UseSystematics && m_counter ) delete m_counter;
  if ( m_physobjsFiller ) delete m_physobjsFiller;
}

TTree* ZeroLeptonCRWT::bookTree(const std::string& treename)
{
  const char* name(treename.c_str());
  TTree* tree = new TTree(name,"ZeroLepton final optimisation");
  tree->SetDirectory(getDirectory());
  bookNTVars(tree,m_ntv,false);
  bookNTReclusteringVars(tree,m_RTntv);
  bookNTCRWTVars(tree,m_crwtntv);
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

  if ( m_UseSystematics ) {
    m_counterRepository = CounterRepository("ZeroLeptonCounter"+m_stringRegion,m_IsSignal,getDirectory());
  }
  else {
    m_counter = new Counter("ZeroLeptonCounter"+m_stringRegion,40,m_IsSignal);
    if(m_doSmallNtuple) m_tree = bookTree(sSR);
  }

}


bool ZeroLeptonCRWT::processEvent(xAOD::TEvent& event)
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
  if ( ! m_IsData && m_period == p8tev ) {
    unsigned int* pveto = 0;
    if ( !store->retrieve<unsigned int>(pveto,"mc12VetoCode").isSuccess() ) throw std::runtime_error("could not retrieve mc12VetoCode");
    veto = *pveto;
    bool* mc12accept = 0;
    if ( !store->retrieve<bool>(mc12accept,"mc12Accept").isSuccess() ) throw std::runtime_error("could not retrieve mc12accept");
    if ( ! *mc12accept ) return true;
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

  // isolated_xxx have overlap removed
  std::vector<ElectronProxy> baseline_electrons, isolated_baseline_electrons, isolated_signal_electrons;
  m_physobjsFiller->FillElectronProxies(baseline_electrons, isolated_baseline_electrons, isolated_signal_electrons);
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
  m_physobjsFiller->FillMuonProxies(baseline_muons, isolated_baseline_muons, isolated_signal_muons);
  // keep only signal muons with Pt>25GeV
  for ( std::vector<MuonProxy>::iterator it = isolated_signal_muons.begin();
	it != isolated_signal_muons.end(); ) {
    if ( it->Pt() < 25000. ) it = isolated_signal_muons.erase(it);
    else it++;
  }
  // FIXME : trigger matching
  std::vector<MuonProxy> trigmatched_muons = isolated_signal_muons;

  // missing ET
  TVector2* missingET = new TVector2(0.,0.);
  if ( ! store->retrieve<TVector2>(missingET,"SUSYMET"+m_suffix).isSuccess() ) throw std::runtime_error("could not retrieve SUSYMET"+m_suffix);
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
  m_ZLUtils.trackMET(event, MET_Track, MET_Track_phi);

  // primary vertex cut
  const xAOD::Vertex* primVertex = ZeroLeptonUtils::GetPrimVtx(event);
  if ( !primVertex ||  !( primVertex->nTrackParticles() > 4) ) return true;
  m_counter->increment(weight,incr++,"Vertex Cut",trueTopo);

  //FIXME isBadMuon not yet implement in SUSYObjDef_xAOD

  // Cosmic muon cut
  if ( m_proxyUtils.CosmicMuon(isolated_baseline_muons) ) return true;
  m_counter->increment(weight,incr++,"CosmicMuons",trueTopo);

  // FIXME fake lepton bkg estimate
  // do we need to re-implement that ?

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
  TLorentzVector leptonTLV;
  int leptonCharge = 0;
  if ( m_isMuonChannel && 
       isolated_baseline_muons.size()==1 && 
       isolated_signal_muons.size()==1 && 
       trigmatched_muons.size()==1 && 
       signalmuonistrigmatched &&
       isolated_baseline_electrons.size()==0 ) {
    oneLepton = true;
    leptonTLV = *(dynamic_cast<TLorentzVector*>(&(isolated_signal_muons[0])));
    leptonCharge = (int)(isolated_signal_muons[0].muon()->charge());
  } 
  else if (  m_isElectronChannel && 
             isolated_baseline_electrons.size()==1 && 
	     isolated_signal_electrons.size()==1   &&
             trigmatched_electrons.size()==1 && 
	     signalelecistrigmatched && 
             isolated_baseline_muons.size()==0) {
    oneLepton = true;
    leptonTLV = *(dynamic_cast<TLorentzVector*>(&(isolated_signal_electrons[0])));
    leptonCharge = (int)(isolated_signal_electrons[0].electron()->charge());
  }
  if ( !oneLepton ) return true;
  m_counter->increment(weight,incr++,"1 Lepton",trueTopo);

  // FIXME: Apply Lepton scale factors
  
  // Add lepton to jets (SR) or MET (VR)
  TVector2 missingETPrime =  *missingET;
  if ( m_isVR ) {
    missingETPrime = missingETPrime + TVector2(leptonTLV.Px(),leptonTLV.Py());
  }
  else {
    good_jets.push_back(JetProxy(leptonTLV,true,true,true,false));
    std::sort(good_jets.begin(),good_jets.end(),PtOrder<JetProxy>);
  }
  double MissingEtPrime = missingETPrime.Mod();
  

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

  // MT cut 
  double mt = std::sqrt( 2.*leptonTLV.Pt()*MissingEt * 
			 (1.-(leptonTLV.Px()*missingET->Px() + leptonTLV.Py()*missingET->Py())/(leptonTLV.Pt()*MissingEt)) );
  if(!(mt >30000 && mt<100000)) return true; 
  m_counter->increment(weight,incr++,"MT cut",trueTopo);


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


  double mT2=-9; 
  //if (good_jets.size()>=2) mT2 = m_proxyUtils.MT2(good_jets,missingETPrime);
  double mT2_noISR=-9;
  //if (nonISR_jets.size()>=2) mT2_noISR = m_proxyUtils.MT2(nonISR_jets,missingETPrime); 
  //out() << " mT2 " << mT2 << " " << mT2_noISR << std::endl; 


  if(m_doSmallNtuple) { 
    unsigned int runnum = RunNumber;
    if ( ! m_IsData ) runnum = mc_channel_number;

    std::vector<float> jetSmearSystW;

    // other cleaning tests
    unsigned int cleaning = 0;
    if (fabs(time[0]) > 5) cleaning += 2;  //t2j use for all SRs

    // FIXME why not in CRWT ?
    //bool chfTileVeto =  m_proxyUtils.chfTileVeto(good_jets);
    //if ( chfTileVeto ) cleaning += 4;

    bool chfVeto = m_proxyUtils.chfVeto(good_jets);
    if ( chfVeto ) cleaning += 8;

    m_proxyUtils.FillNTVars(m_ntv, runnum, EventNumber, veto, weight, normWeight, *pileupWeights, genWeight,ttbarWeightHT,ttbarWeightPt2,ttbarAvgPt,WZweight, btag_weight, ctag_weight, b_jets.size(), c_jets.size(), MissingEtPrime, phi_met, Meff, meffincl, minDphi, RemainingminDPhi, good_jets, trueTopo, cleaning, time[0],jetSmearSystW,0, 0., 0.);
    m_proxyUtils.FillNTReclusteringVars(m_RTntv, good_jets);
    FillCRWTVars(m_crwtntv,leptonTLV,*missingET,leptonCharge);

    m_tree->Fill();
  }
  return true;
}

void ZeroLeptonCRWT::finish()
{
  if ( m_UseSystematics ) {
    out() << m_counterRepository << std::endl;
  } 
  else {
    out() << *m_counter << std::endl;
  }
}


void ZeroLeptonCRWT::FillCRWTVars(NTCRWTVars& crwtvars, const TLorentzVector& lepton, const TVector2& metv, int lepsign)
{
  crwtvars.Reset();
  crwtvars.lep1sign = lepsign;
  crwtvars.lep1Pt  = lepton.Pt();
  crwtvars.lep1Eta = lepton.Eta();
  crwtvars.lep1Phi = lepton.Phi();

  double met = std::sqrt(metv.Px()*metv.Px()+metv.Py()*metv.Py());
  crwtvars.mt = std::sqrt(2.*lepton.Pt() * met * (1. - (lepton.Px()*metv.Px() + lepton.Py()*metv.Py())/(lepton.Pt()*met)));

  double wpx = lepton.Px()+metv.Px();
  double wpy = lepton.Py()+metv.Py();
  crwtvars.Wpt = std::sqrt(wpx*wpx+wpy*wpy);
}

ClassImp(ZeroLeptonCRWT);

