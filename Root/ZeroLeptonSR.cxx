#include "ZeroLeptonRun2/ZeroLeptonSR.h"
#include "ZeroLeptonRun2/PhysObjProxies.h"


#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TActiveStore.h"
#include "xAODRootAccess/TStore.h"
#include "xAODTracking/Vertex.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODMissingET/MissingETContainer.h"
#include "PATInterfaces/SystematicSet.h"
#include "cafe/Processor.h"
#include "cafe/Controller.h"
#include "cafe/Config.h"
#include "TDirectory.h"
#include "TVector2.h"

#include <iostream>
#include <stdexcept>

ZeroLeptonSR::ZeroLeptonSR(const char *name)
  : cafe::Processor(name), 
    m_tree(0), 
    m_stringRegion("SRAll"), 
    m_doSmallNtuple(true),
    m_IsData(false),
    m_IsSignal(false),
    m_UseSystematics(false),
    m_period(INVALID),
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
  m_UseSystematics = config.get("UseSystematics",false);
  m_period = periodFromString(config.get("Period","p13tev"));
  if ( m_period == p7tev ) throw(std::domain_error("ZeroLeptonSR does not support the 7tev run period"));
  if ( m_period == INVALID ) throw(std::domain_error("ZeroLeptonSR: invalid run period specified"));

  std::string cutfile = config.get("cutfile","None");
  if ( cutfile == "None" ) throw(std::domain_error("ZeroLeptonSR: invalid cut file specified"));
  m_cutVal.ReadCutValues(cutfile);


  m_suffix = config.get("suffix","");
  m_physobjsFiller = new PhysObjProxyFiller(20000.f,10000.f,10000.f,m_suffix);
  m_proxyUtils = PhysObjProxyUtils(m_IsData);
  m_ZLUtils = ZeroLeptonUtils(m_IsData);
}

ZeroLeptonSR::~ZeroLeptonSR()
{
  if ( !m_UseSystematics && m_counter ) delete m_counter;
  if ( m_physobjsFiller ) delete m_physobjsFiller;
}

TTree* ZeroLeptonSR::bookTree(const std::string& treename)
{
  const char* name(treename.c_str());
  TTree* tree = new TTree(name,"ZeroLepton final optimisation");
  tree->SetDirectory(getDirectory());
  bookNTVars(tree,m_ntv,false);
  bookNTReclusteringVars(tree,m_RTntv);
  bookNTExtraVars(tree,m_extrantv);
  return tree;
}

TTree* ZeroLeptonSR::getTree(const std::string& treename)
{
  std::map<std::string,TTree*>::const_iterator pos = m_treeRepository.find(treename);
  if ( pos == m_treeRepository.end() ) {
    pos = m_treeRepository.insert(std::make_pair(treename,bookTree(treename))).first;
  }
  return pos->second;
}

void ZeroLeptonSR::begin()
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


bool ZeroLeptonSR::processEvent(xAOD::TEvent& event)
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
  if ( ! m_IsData && m_period == p8tev ) {
    unsigned int* pveto = 0;
    if ( !store->retrieve<unsigned int>(pveto,"mc12VetoCode").isSuccess() ) throw std::runtime_error("could not retrieve mc12VetoCode");
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

  // isolated_xxx have overlap removed
  std::vector<ElectronProxy> baseline_electrons, isolated_baseline_electrons, isolated_signal_electrons;
  m_physobjsFiller->FillElectronProxies(baseline_electrons, isolated_baseline_electrons, isolated_signal_electrons);

  // isolated_xxx have overlap removed
  std::vector<MuonProxy> baseline_muons, isolated_baseline_muons, isolated_signal_muons;
  m_physobjsFiller->FillMuonProxies(baseline_muons, isolated_baseline_muons, isolated_signal_muons);

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


  // 0 lepton
  if ( !isolated_baseline_electrons.empty() ) return true;
  if ( !isolated_baseline_muons.empty() ) return true;
  m_counter->increment(weight,incr++,"0 Lepton",trueTopo);
 
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
  if (!(MissingEt > m_cutVal.m_cutEtMiss)) return true;
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
  double phi_met = TMath::ATan2(missingET->Y(),missingET->X());
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
				       MissingEt,
				       m_cutVal.m_cutJetPt1,
				       m_cutVal.m_cutJetPt4);
    //out() << " Meff[" << i <<"] " << Meff[i] << std::endl;
  }
  double meffincl = m_proxyUtils.Meff(good_jets,
					     good_jets.size(),
					     MissingEt,
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
  m_proxyUtils.GetAlphaISRVar(good_jets,MissingEt,ISRvar_alpha);
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
  m_proxyUtils.GetISRJet(good_jets,isr_jet_indices,MissingEt,phi_met,"squark",false);  
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
  //if (good_jets.size()>=2) mT2 = m_proxyUtils.MT2(good_jets,*missingET);
  double mT2_noISR=-9;
  //if (nonISR_jets.size()>=2) mT2_noISR = m_proxyUtils.MT2(nonISR_jets,*missingET); 
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
			      missingET->X(),
			      missingET->Y(),
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

    const xAOD::MissingETContainer* metcontainer = 0;
    float metLHTOPOx = 0.f;
    float metLHTOPOy = 0.f;
    if( event.retrieve( metcontainer, "MET_LocHadTopo" ).isSuccess() ) {
      xAOD::MissingETContainer::const_iterator met_it = metcontainer->find("LocHadTopo");
      metLHTOPOx = (*met_it)->mpx();
      metLHTOPOy = (*met_it)->mpy();
    }

    m_proxyUtils.FillNTVars(m_ntv, runnum, EventNumber, veto, weight, normWeight, *pileupWeights, genWeight,ttbarWeightHT,ttbarWeightPt2,ttbarAvgPt,WZweight, btag_weight, ctag_weight, b_jets.size(), c_jets.size(), MissingEt, phi_met, Meff, meffincl, minDphi, RemainingminDPhi, good_jets, trueTopo, cleaning, time[0],jetSmearSystW,0, 0., 0., metLHTOPOx, metLHTOPOy);


    m_proxyUtils.FillNTExtraVars(m_extrantv, mT2,mT2_noISR,gaminvRp1 ,shatR ,mdeltaR ,cosptR ,gamma_R,dphi_BETA_R , dphi_leg1_leg2 , costhetaR ,dphi_BETA_Rp1_BETA_R,gamma_Rp1,costhetaRp1,Ap);
   
    m_proxyUtils.FillNTReclusteringVars(m_RTntv,good_jets);
    
    m_tree->Fill();
  }
  return true;
}

void ZeroLeptonSR::finish()
{
  if ( m_UseSystematics ) {
    out() << m_counterRepository << std::endl;
  } 
  else {
    out() << *m_counter << std::endl;
  }
}

ClassImp(ZeroLeptonSR);

