
#include "ZeroLeptonRun2/ZeroLeptonDataDrivenQCD.h"
#include "ZeroLeptonRun2/PhysObjProxies.h"
#include "cafe/Config.h"
#include "cafe/ParseRun.h"
#include "JetSmearing/JetMCSmearingTool.h"
#include "METUtilities/METMaker.h"
#include "xAODJet/JetContainer.h"
#include "xAODJet/JetAuxContainer.h"
#include "xAODMissingET/MissingET.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODMissingET/MissingETAuxContainer.h"
#include "xAODMissingET/MissingETComposition.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TActiveStore.h"
#include "xAODRootAccess/TStore.h"
#include "xAODCore/ShallowCopy.h"
#include "xAODBase/IParticleHelpers.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODTracking/Vertex.h"
#include "xAODMissingET/MissingETContainer.h"
#include "TRandom.h"
#include "TDirectory.h"
#include "TVector2.h"
#include "TFile.h"
#include "TH2F.h"
#include <stdexcept>

static SG::AuxElement::Decorator<char> dec_smeared("smearjet");
static SG::AuxElement::Decorator<char> dec_baseline("baseline");

ZeroLeptonDataDrivenQCD::ZeroLeptonDataDrivenQCD(const char *name): 
  cafe::Processor(name), 
  m_counter(0),
  m_saveSeedEvents(false),
  m_jetkey(),
  m_NumberSmear(10),
  m_SeedOption(1),
  m_QCDSeedSelector(),
  m_physobjsFiller(PhysObjProxyFiller(20000.f,10000.f,10000.f,25000.f,"",false,"","")),
  m_proxyUtils(true),
  m_ZLUtils(true),
  m_smearingTool(0),
  m_smearFnFile(0),
  m_smearFnFile2(0),
  m_GaussianCoreSmearingType("none"),
  m_containers(),
  m_cutVal(),
  m_processors()
{
  cafe::Config config(name);
  m_saveSeedEvents = config.get("SaveSeedEvents",m_saveSeedEvents);
  m_GaussianCoreSmearingType = config.get("GaussianCoreSmearingType","optimal");
  m_containers = config.getVString("Containers");
  m_jetkey = config.get("JetCollectionKey","SUSYJets");
  m_NumberSmear = config.get("NumberSmear", 10);
  m_SeedOption = config.get("SeedOption", 1);

  m_ZLUtils = ZeroLeptonUtils(true);//always true for JetSmearing

  InitialiseSmearing();

  m_prescaleTool = new PreScaleTool ();

  std::string cutfile = config.get("cutfile","None");
  if ( cutfile == "None" ) throw(std::domain_error("ZeroLeptonDataDrivenQCD: invalid cut file specified"));
  m_cutVal.ReadCutValues(cutfile);

  // list of child processors
  std::string run = config.get("Run","");
  if(run != "") {
    cafe::ParseRun parser;
    add(parser.parse(run));
  }

  m_metMaker = new met::METMaker("METMaker_ZeroLeptonRun2");
  m_metMaker->setProperty("ORCaloTaggedMuons", true);
  m_metMaker->initialize();
}

void ZeroLeptonDataDrivenQCD::begin()
{
 TDirectory* dir = getDirectory();
  if ( !dir ) {
    throw std::runtime_error("No root TDirectory defined in processor ZeroLeptonSR("+this->name()+")");
  }

  m_counter = new Counter("ZeroLeptonDataDrivenQCDCounter",40); 

  std::for_each(m_processors.begin(),m_processors.end(),
		std::mem_fun(&Processor::begin));
}

ZeroLeptonDataDrivenQCD::~ZeroLeptonDataDrivenQCD()
{
  if ( m_counter ) delete m_counter;  
  if (m_smearingTool) delete m_smearingTool;
  if (m_smearFnFile ) delete m_smearFnFile;
  if (m_smearFnFile2) delete m_smearFnFile2;

  for ( std::list<Processor*>::iterator it = m_processors.begin();
	it != m_processors.end();
	++it ) {
    delete *it;
  }
}

bool ZeroLeptonDataDrivenQCD::processEvent(xAOD::TEvent& event)
{
  float weight=1;

  // access the transient store
  xAOD::TStore* store = xAOD::TActiveStore::store();

  // counters
  int incr=0;
  m_counter->increment(1.,incr++,"NbOfEvents");
  m_counter->increment(weight,incr++,"runNumber");

  // eventInfo
  const xAOD::EventInfo* eventInfo = 0;
  if ( ! event.retrieve( eventInfo, "EventInfo").isSuccess() ) throw std::runtime_error("Could not retrieve EventInfo");
  uint32_t RunNumber = eventInfo->runNumber();
  unsigned long long EventNumber = eventInfo->eventNumber();

  if(eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ){
    throw std::logic_error("ZeroLeptonDataDrivenQCD can only be run on real data!");
  }

  bool* passGRL = 0;
  if ( !store->retrieve<bool>(passGRL,"passGRL").isSuccess() ) throw std::runtime_error("could not retrieve passGRL");
  if ( ! *passGRL ) return true;
  m_counter->increment(weight,incr++,"GRL");

  // primary vertex cut
  const xAOD::Vertex* primVertex = ZeroLeptonUtils::GetPrimVtx(event);
  if ( !primVertex ) return true;
  m_counter->increment(weight,incr++,"Vertex Cut");

  // These jets have overlap removed
  std::vector<JetProxy> good_jets, bad_jets, b_jets, c_jets;
  m_physobjsFiller.FillJetProxies(good_jets,bad_jets,b_jets);
  std::vector<float> btag_weight(7,1.); // not implemented in SUSYTools
  std::vector<float> ctag_weight(7,1.); // not implemented in SUSYTools

  // isolated_xxx have overlap removed
  std::vector<ElectronProxy> baseline_electrons, isolated_baseline_electrons, isolated_signal_electrons;
  m_physobjsFiller.FillElectronProxies(baseline_electrons, isolated_baseline_electrons, isolated_signal_electrons);

  // isolated_xxx have overlap removed
  std::vector<MuonProxy> baseline_muons, isolated_baseline_muons, isolated_signal_muons;
  m_physobjsFiller.FillMuonProxies(baseline_muons, isolated_baseline_muons, isolated_signal_muons);

  // missing ET
  TVector2* missingET = new TVector2(0.,0.);
  if ( ! store->retrieve<TVector2>(missingET,"SUSYMET").isSuccess() ) throw std::runtime_error("could not retrieve SUSYMET");
  double MissingEt = missingET->Mod();

  if ( good_jets.empty() ) return true;
  m_counter->increment(weight,incr++,"One jet");
  //FIXME need to get list of triggers and plateaux
  std::string single_jet_triggers [12] = {"HLT_j400","HLT_j360","HLT_j320","HLT_j260","HLT_j200","HLT_j175",
					"HLT_j150","HLT_j110","HLT_j85","HLT_j60","HLT_j25","HLT_j15"};
  std::vector<std::pair<std::string,bool > > triggersPass;
  for(const std::string &sjt : single_jet_triggers) {
    triggersPass.push_back( std::make_pair (sjt, (int)eventInfo->auxdata<char>(sjt)==1));
  }
  std::string best_trigger = m_prescaleTool->whichTrigger(triggersPass,good_jets[0].Pt());
  float triggerWeight = m_prescaleTool->getTriggerPrescale(best_trigger,RunNumber);
  if ( triggerWeight<1.0 ) return true;

  weight *= triggerWeight;
  m_counter->increment(weight,incr++,"Trigger");

  // LAr, Tile, reco problems in data
  bool* badDetectorQuality = 0 ; 
  if ( !store->retrieve<bool>(badDetectorQuality,"badDetectorQuality").isSuccess() ) throw std::runtime_error("could not retrieve badDetectorQuality");
  if ( *badDetectorQuality ) return true;
  m_counter->increment(weight,incr++,"Detector cleaning");

  // 0 lepton
  if ( !isolated_baseline_electrons.empty() ) return true;
  if ( !isolated_baseline_muons.empty() ) return true;
  m_counter->increment(weight,incr++,"0 Lepton");

  // apply seed selection cut?

  const xAOD::MissingETContainer *originalmet = 0;
  if ( ! store->retrieve(originalmet, "MET_ZL").isSuccess() ) {
    throw std::runtime_error("Unable to retrieve MET_ZL");
  }
  xAOD::MissingETContainer::const_iterator mettst_refjet_it = originalmet->find("RefJet");
  if (mettst_refjet_it==originalmet->end()) out() << "No RefJet inside MET container" << std::endl;
  double sumet = (*mettst_refjet_it)->sumet();
  double met_significance = (MissingEt-8000.)/sqrt(sumet*1e3);
  int Nbjets = b_jets.size();

  if( met_significance > (0.5 + 0.1*Nbjets) ) return true;

  double average_jet_pt = 0.;
  for(size_t i=0; i<good_jets.size(); i++) {
    average_jet_pt += good_jets[i].Pt();
  }
  average_jet_pt /= (double)good_jets.size();

  float met_over_pt = MissingEt/average_jet_pt;
  if( m_SeedOption==2 ){
    if( met_over_pt>0.2 ) return true;
  }

  m_counter->increment(weight,incr++,"seed selection");

  if ( m_saveSeedEvents ) {
    write_xAOD_event();
    return true;
  }

  // ---------------------   Apply smearing & recompute MET

  eventInfo->auxdecor<float>("triggerWeight") = triggerWeight;
  eventInfo->auxdecor<float>("met_over_pt") = met_over_pt;

  // get Jet container with calibrated jets
  xAOD::JetContainer* jets = 0;
  if ( !store->retrieve(jets, "SUSYJets").isSuccess() ) {
    throw std::runtime_error("Could not retrieve JetContainer with key SUSYJets");
  }
  xAOD::JetContainer::iterator jet_itr = jets->begin();
  xAOD::JetContainer::iterator jet_end = jets->end();
  for(; jet_itr != jet_end; ++jet_itr ) {
    if( (*jet_itr)->auxdata<char>("passOR")==1 && (*jet_itr)->auxdata<char>("baseline")==1 ){
      dec_smeared(**jet_itr) = true;
    }else{
      dec_smeared(**jet_itr) = false;
    }
  }
  
  // get MET recomputed after jet calibration
  const xAOD::MissingETContainer* metcontainer = 0;
  if( !event.retrieve( metcontainer, "MET_Core_AntiKt4EMTopo" ).isSuccess() ) {
    throw std::runtime_error("Cuold not retrieve METContainer with key MET_Core_AntiKt4EMTopo");
  }
  const xAOD::MissingETAssociationMap *metMap = 0;
  if ( !event.retrieve( metMap, "METAssoc_AntiKt4EMTopo") ) {
    throw std::runtime_error("Could not retrieve METAssociationMap with key METAssoc_AntiKt4EMTopo");
  }

  gRandom->SetSeed(RunNumber+EventNumber);

  std::vector<std::unique_ptr<SmearData> > smrMc;
  m_smearingTool->DoSmearing(smrMc,*jets);

  xAOD::JetContainer* smearedjets = new xAOD::JetContainer;
  xAOD::AuxContainerBase* smearedjetsaux = new xAOD::AuxContainerBase;
  smearedjets->setStore( smearedjetsaux );
  if ( ! store->record( smearedjets, m_jetkey+"_smeared").isSuccess() ){
    throw std::runtime_error("Could not store "+m_jetkey+"_smeared");
  }
  if ( ! store->record( smearedjetsaux, m_jetkey+"_smearedAux.").isSuccess() ){
    throw std::runtime_error("Could not store "+m_jetkey+"_smearedAux.");
  } 

  TVector2* missingET_smeared = new TVector2(0., 0.);
  if ( ! store->record(missingET_smeared,"SUSYMET_smeared").isSuccess() ) throw std::runtime_error("Could not store SUSYMET with tag SUSYMET_smeared");

  xAOD::MissingETContainer* SMRmetContainer = new xAOD::MissingETContainer;
  xAOD::MissingETAuxContainer* SMRmetContainerAux = new xAOD::MissingETAuxContainer;
  SMRmetContainer->setStore(SMRmetContainerAux);
  if ( ! store->record(SMRmetContainer,"SMRMet").isSuccess() ) throw std::runtime_error("Could not store SMRMet");
  if ( ! store->record(SMRmetContainerAux, "SMRMetAux.").isSuccess() ) throw std::runtime_error("Could not store SMRMETAux");

  //----------------------------------------   Electrons
  const xAOD::ElectronContainer* electrons = 0;
  const xAOD::ShallowAuxContainer* electrons_aux = 0;
  if ( ! store->retrieve(electrons, "SUSYElectrons").isSuccess() ) throw std::runtime_error("Could not retrieve SUSYElectrons");
  if ( ! store->retrieve(electrons_aux, "SUSYElectronsAux.").isSuccess() ) throw std::runtime_error("Could not retrieve SUSYElectronsAux");

  //----------------------------------------   Muons
  const xAOD::MuonContainer* muons = 0;
  const xAOD::ShallowAuxContainer* muons_aux = 0;
  if ( ! store->retrieve(muons, "SUSYMuons").isSuccess() ) throw std::runtime_error("Could not retrieve SUSYMuons");
  if ( ! store->retrieve(muons_aux, "SUSYMuonsAux.").isSuccess() ) throw std::runtime_error("Could not retrieve SUSYMuonsAux");

  //----------------------------------------   Photons
  const xAOD::PhotonContainer* photons = 0;
  xAOD::ShallowAuxContainer* photons_aux = 0;
  if ( ! store->retrieve(photons, "SUSYPhotons").isSuccess() ) throw std::runtime_error("Could not retrieve SUSYPhotons");
  if ( ! store->retrieve(photons_aux, "SUSYPhotonsAux.").isSuccess() ) throw std::runtime_error("Could not retrieve SUSYPhotonsAux");

  std::vector<std::unique_ptr<SmearData> >::iterator firstSmear = smrMc.begin();
  std::vector<std::unique_ptr<SmearData> >::iterator lastSmear = smrMc.end();
  for(; firstSmear != lastSmear; ++firstSmear){
    smearedjets->clear();
    xAOD::JetContainer::iterator firstJet = (*firstSmear)->jetContainer->begin();
    xAOD::JetContainer::iterator lastJet = (*firstSmear)->jetContainer->end();
    for(; firstJet != lastJet; ++firstJet){
      if( (*firstJet)->auxdata< char >("smearjet") == true ){
	xAOD::Jet *jet = new xAOD::Jet;
	smearedjets->push_back( jet );
	*jet = **firstJet;
      }
    }
    SMRmetContainer->clear();
    GetMET(*SMRmetContainer, metcontainer, metMap, (*firstSmear)->jetContainer, electrons, muons, photons, 0, true);
    
    xAOD::MissingETContainer::const_iterator met_it = SMRmetContainer->find("Final");
    if ( met_it == SMRmetContainer->end() ) throw std::runtime_error("Could not find Final MET after running  GetMET");
    missingET_smeared->Set((*met_it)->mpx(), (*met_it)->mpy());
   
    for(std::list<cafe::Processor*>::iterator it = m_processors.begin();
	it != m_processors.end();
	++it) {
      (*it)->incEventCount();
      if(!(*it)->processEvent(event)) return false;
    }
  }

  //delete jets;
  //   delete jets_aux;
  /*delete muons;
  delete muons_aux;
  delete electrons;
  delete electrons_aux;
  delete photons;
  delete photons_aux;*/
  //delete rebuiltmetc;
  //    delete rebuiltmetcAux;

  return true;
}


void ZeroLeptonDataDrivenQCD::GetResponseMaps(TFile* smearFnFile,TH2F*& nominal, TH2F*& tailHigh, TH2F*& tailLow)
{
  std::string histNameBase = "responseEJES_p2419SUSY11";
  std::string histName;
  TH2F* tempMap = 0;

  histName = histNameBase; //+std::string("optimal");
  tempMap = (TH2F*)smearFnFile->Get(histName.c_str());
  if ( !tempMap )  throw std::invalid_argument("ZeroLeptonDataDrivenQCD: could not find histogram "+histName);
  nominal = (TH2F*)tempMap->Clone("map");

  /** not implemented now
  histName = histNameBase+std::string("high");
  tempMap = (TH2F*)smearFnFile->Get(histName.c_str());
  if ( !tempMap )  throw std::invalid_argument("ZeroLeptonDataDrivenQCD: could not find histogram "+histName);
  tailHigh = (TH2F*)tempMap->Clone("map");

  histName = histNameBase+std::string("low");
  tempMap = (TH2F*)smearFnFile->Get(histName.c_str());
  if ( !tempMap )  throw std::invalid_argument("ZeroLeptonDataDrivenQCD: could not find histogram "+histName);
  tailLow = (TH2F*)tempMap->Clone("map");
  **/
}

void ZeroLeptonDataDrivenQCD::InitialiseSmearing()
{
  char* ROOTCOREDIR = std::getenv("ROOTCOREBIN");

  std::string filebveto = std::string(ROOTCOREDIR) + "/data/JetSmearing/MC15/R_map2015_bveto_OP77_EJES_p2419SUSY11_MarchHADD.root";
  std::string filebtag = std::string(ROOTCOREDIR) + "/data/JetSmearing/MC15/R_map2015_btag_OP77_EJES_p2419SUSY11_MarchHADD.root";

  m_smearFnFile  = new TFile(filebveto.c_str());
  if ( !m_smearFnFile || m_smearFnFile->IsZombie() ) {
    throw std::invalid_argument("ZeroLeptonDataDrivenQCD: could not open root file "+filebveto);    
  }
  m_smearFnFile2 = new TFile(filebtag.c_str());
  if ( !m_smearFnFile2 || m_smearFnFile2->IsZombie() ) {
    throw std::invalid_argument("ZeroLeptonDataDrivenQCD: could not open root file "+filebtag);    
  }

  TH2F *bveto_nominal=0, *bveto_tailHigh=0, *bveto_tailLow=0;
  GetResponseMaps(m_smearFnFile,bveto_nominal,bveto_tailHigh,bveto_tailLow);

  TH2F *btag_nominal=0, *btag_tailHigh=0, *btag_tailLow=0;
  GetResponseMaps(m_smearFnFile2,btag_nominal,btag_tailHigh,btag_tailLow);

  m_smearingTool  = new SUSY::JetMCSmearingTool("ZLJetMCSmearingTool");
  m_smearingTool->SetResponseMaps(bveto_nominal,btag_nominal);
  /** not implemeted now
  m_smearingTool->AddSystematicResponseMaps(bveto_tailHigh,btag_tailHigh);
  m_smearingTool->AddSystematicResponseMaps(bveto_tailLow,btag_tailLow);
  **/
  m_smearingTool->setProperty("NumberOfSmearedEvents",m_NumberSmear);
  m_smearingTool->setProperty("DoBJetSmearing",true);
  m_smearingTool->setProperty("UseTailWeights",false);
  m_smearingTool->setProperty("DoGaussianCoreSmearing",false);
  m_smearingTool->setProperty("GaussianCoreSmearingType",m_GaussianCoreSmearingType);
  m_smearingTool->setProperty("ShiftMeanType", "none");
  m_smearingTool->setProperty("DoPhiSmearing",false);
  
  m_smearingTool->initialize();
}

void ZeroLeptonDataDrivenQCD::finish()
{
  std::for_each(m_processors.begin(),m_processors.end(),
		std::mem_fun(&Processor::finish));
  out() << *m_counter << std::endl;
}

void ZeroLeptonDataDrivenQCD::inputFileOpened(TFile *file) 
{
  for(std::list<Processor*>::iterator it =m_processors.begin();
      it !=m_processors.end();
      ++it) {
    (*it)->inputFileOpened(file);
  }
}

void ZeroLeptonDataDrivenQCD::inputFileClosing(TFile *file)
{
  for(std::list<Processor*>::iterator it =m_processors.begin();
      it !=m_processors.end();
      ++it) {
    (*it)->inputFileClosing(file);
  }
}

bool ZeroLeptonDataDrivenQCD::add(const std::list<cafe::Processor*>& procs)
{
  for (std::list<Processor*>::const_iterator it = procs.begin();
       it != procs.end();
       ++it) {
    add(*it);
  }
  return true;
}

bool ZeroLeptonDataDrivenQCD::add(cafe::Processor *proc)
{
  if(proc) {
    //proc->setParent(this);
    setParent(proc,this);
    if(debug() > 0) {
      proc->setDebug(debug());
    }
    out() << "ZeroLeptonDataDrivenQCD[" << name() << "]: Adding " << proc->fullName() << std::endl;
   m_processors.push_back(proc);
    return true;
  } else {
    return false;
  }
}


bool ZeroLeptonDataDrivenQCD::GetMET(xAOD::MissingETContainer& met,
				     const xAOD::MissingETContainer* metcore,
				     const xAOD::MissingETAssociationMap* metMap,		
				     const xAOD::JetContainer* jet,
				     const xAOD::ElectronContainer* elec,
				     const xAOD::MuonContainer* muon,
				     const xAOD::PhotonContainer* gamma,
				     const xAOD::TauJetContainer* taujet,
				     bool doTST, bool doJVTCut,
				     const xAOD::IParticleContainer* invis)
{
  std::string m_eleTerm = "RefEle";
  std::string m_gammaTerm = "RefGamma";
  std::string m_tauTerm = "RefTau";
  std::string m_jetTerm = "RefJet";
  std::string m_muonTerm = "Muons";
  std::string m_outMETTerm = "Final";

  std::string softTerm = "SoftClus";
  if (doTST) {
    softTerm = "PVSoftTrk";
  }

  metMap->resetObjSelectionFlags();

  // allow creation of proxy MET by flagging objects for "neutrino/ification" as already selected
  if (invis) {
    m_metMaker->markInvisible(invis, metMap);
  }

  if (elec) {
    ConstDataVector<xAOD::ElectronContainer> metelectron(SG::VIEW_ELEMENTS);
    for (const auto& el : *elec) {
      if (dec_baseline(*el)) {
        bool veto(false);
        if (invis) {
          for (const auto& ipart : *invis) {
            if (ipart == el) {veto = true; break;}
          }
        }
        if (!veto) metelectron.push_back(el);
      }
    }
    m_metMaker->rebuildMET(m_eleTerm, xAOD::Type::Electron, &met, metelectron.asDataVector(), metMap);
  }

  if (gamma) {
    ConstDataVector<xAOD::PhotonContainer> metgamma(SG::VIEW_ELEMENTS);
    for (const auto& ph : *gamma) {
      if (dec_baseline(*ph)) {
        bool veto(false);
        if (invis) {
          for (const auto& ipart : *invis) {
            if (ipart == ph) {veto = true; break;}
          }
        }
        if (!veto) metgamma.push_back(ph);
      }
    }
    m_metMaker->rebuildMET(m_gammaTerm, xAOD::Type::Photon, &met, metgamma.asDataVector(), metMap);
  }

  if (taujet) {
    ConstDataVector<xAOD::TauJetContainer> mettau(SG::VIEW_ELEMENTS);
    for (const auto& tau : *taujet) {
      if (dec_baseline(*tau)) {
        bool veto(false);
        if (invis) {
          for (const auto& ipart : *invis) {
            if (ipart == tau) {veto = true; break;}
          }
        }
        if (!veto) mettau.push_back(tau);
      }
    }
    m_metMaker->rebuildMET(m_tauTerm, xAOD::Type::Tau, &met, mettau.asDataVector(), metMap);
  }

  if (muon) {
    ConstDataVector<xAOD::MuonContainer> metmuon(SG::VIEW_ELEMENTS);
    for (const auto& mu : *muon) {
      bool veto(false);
      if (dec_baseline(*mu)) {
        if (invis) {
          for (const auto& ipart : *invis) {
            if (ipart == mu) {veto = true; break;}
          }
        }
        if (!veto) metmuon.push_back(mu);
      }
    }
    m_metMaker->rebuildMET(m_muonTerm, xAOD::Type::Muon, &met, metmuon.asDataVector(), metMap);
  }

  if (!jet) {
    out() << "WARNING: " << "Invalid jet container specified for MET rebuilding!"<< std::endl;
    return StatusCode::SUCCESS;
  }
  m_metMaker->rebuildJetMET(m_jetTerm, softTerm, &met, jet, metcore, metMap, doJVTCut);

  /* MET systematics is needed?
  if (!m_IsData) {
    if ( m_metSystTool->applyCorrection(*met[softTerm]) != CP::CorrectionCode::Ok ) {
      out() << "WARNING: " << "GetMET: Failed to apply MET soft term systematics."<< std::endl;
    }
    out() << "VERBOSE: " << "New soft term value: " << met[softTerm]->met() << std::endl;
  }
  */
  m_metMaker->buildMETSum(m_outMETTerm, &met, met[softTerm]->source());

  //out() << "VERBOSE: " <<  "Rebuilt MET: Missing Et (x,y): (" << met[m_outMETTerm]->mpx() << "," <<  met[m_outMETTerm]->mpy() << ")"<< std::endl;

  return true;
}


ClassImp(ZeroLeptonDataDrivenQCD);

