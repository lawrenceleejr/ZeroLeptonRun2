
#include "ZeroLeptonRun2/ZeroLeptonDataDrivenQCD.h"
#include "ZeroLeptonRun2/PhysObjProxies.h"
#include "cafe/Config.h"
#include "cafe/ParseRun.h"
#include "JetSmearing/JetMCSmearingTool.h"
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

ZeroLeptonDataDrivenQCD::ZeroLeptonDataDrivenQCD(const char *name): 
  cafe::Processor(name), 
  m_counter(0),
  m_seed(0),
  m_saveSeedEvents(false),
  m_jetkey(),
  m_QCDSeedSelector(),
  m_physobjsFiller(PhysObjProxyFiller(20000.f,10000.f,10000.f,"",false,"")),
  m_proxyUtils(true),
  m_ZLUtils(true, NotADerivation),
  m_smearingTool(0),
  m_smearFnFile(0),
  m_smearFnFile2(0),
  m_GaussianCoreSmearingType("none"),
  m_containers(),
  m_cutVal(),
  m_processors(),
  m_derivationTag(INVALID_Derivation)
{
  cafe::Config config(name);
  m_seed = config.get("RandomSeed",m_seed);
  m_saveSeedEvents = config.get("SaveSeedEvents",m_saveSeedEvents);
  m_GaussianCoreSmearingType = config.get("GaussianCoreSmearingType","optimal");
  m_containers = config.getVString("Containers");
  m_jetkey = config.get("JetCollectionKey","xxxx");

  m_derivationTag = derivationTagFromString(config.get("DerivationTag",""));
  if ( m_derivationTag == INVALID_Derivation ) throw(std::domain_error("ZeroLeptonSR: invalid derivation tag specified"));

  m_ZLUtils = ZeroLeptonUtils(true, m_derivationTag);

  gRandom->SetSeed(gRandom->GetSeed()+m_seed);
  InitialiseSmearing();

  std::string cutfile = config.get("cutfile","None");
  if ( cutfile == "None" ) throw(std::domain_error("ZeroLeptonDataDrivenQCD: invalid cut file specified"));
  m_cutVal.ReadCutValues(cutfile);

  // list of child processors
  std::string run = config.get("Run","");
  if(run != "") {
    cafe::ParseRun parser;
    add(parser.parse(run));
  }

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
  //uint32_t RunNumber = eventInfo->runNumber();
  //unsigned long long EventNumber = eventInfo->eventNumber();

  if(eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ){
    throw std::logic_error("ZeroLeptonDataDrivenQCD can only be run on real data!");
  }

  bool* passGRL = 0;
  if ( !store->retrieve<bool>(passGRL,"passGRL").isSuccess() ) throw std::runtime_error("could not retrieve passGRL");
  if ( ! *passGRL ) return true;
  m_counter->increment(weight,incr++,"GRL");

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
  //FIXME no trigger information ! 
  float triggerWeight = 1.0;
  //if (!m_QCDSeedSelector.passTrigger_2012Data(r->RunNumber,r->EF_j460_a4tchad,r->EF_j360_a4tchad,r->EF_j280_a4tchad,r->EF_j220_a4tchad,r->EF_j180_a4tchad,r->EF_j145_a4tchad,r->EF_j110_a4tchad,r->EF_j80_a4tchad,r->EF_j55_a4tchad,m_util->numJets(), m_util->jet_Pt(isolated_good_jets_indices[0]), triggerWeight)) return;
  out() << " lead jet pT " <<  good_jets[0].Pt() << std::endl;
  if ( good_jets[0].Pt() < 410000. ) return true;

  weight *= triggerWeight;
  m_counter->increment(weight,incr++,"Trigger");

  // bad jet veto
  if ( !bad_jets.empty() ) return true;
  m_counter->increment(weight,incr++,"JetCleaning");

  // LAr, Tile, reco problems in data
  bool* badDetectorQuality = 0 ; 
  if ( !store->retrieve<bool>(badDetectorQuality,"badDetectorQuality").isSuccess() ) throw std::runtime_error("could not retrieve badDetectorQuality");
  if ( *badDetectorQuality ) return true;
  m_counter->increment(weight,incr++,"Detector cleaning");

  // primary vertex cut
  const xAOD::Vertex* primVertex = ZeroLeptonUtils::GetPrimVtx(event);
  if ( !primVertex ||  !( primVertex->nTrackParticles() > 4) ) return true;
  m_counter->increment(weight,incr++,"Vertex Cut");

  // Cosmic muon cut
  if ( m_proxyUtils.CosmicMuon(isolated_baseline_muons) ) return true;
  m_counter->increment(weight,incr++,"CosmicMuons");

  //FIXME isBadMuon not yet implement in SUSYObjDef_xAOD

  // bad muons for MET cut: based on non isolated muons
  if ( m_proxyUtils.isbadMETmuon(baseline_muons, MissingEt, *missingET) ) return true;
  m_counter->increment(weight,incr++,"IsBadMETMuon");

  // bad Tile cut
  if ( m_proxyUtils.badTileVeto(good_jets,*missingET)) return true;
  m_counter->increment(weight,incr++,"Bad Tile Veto");

  // no longer used
  m_counter->increment(weight,incr++,"Negative-cell cleaning");

  // 0 lepton
  if ( !isolated_baseline_electrons.empty() ) return true;
  if ( !isolated_baseline_muons.empty() ) return true;
  m_counter->increment(weight,incr++,"0 Lepton");
 
  //FIXME: in the Run 1 code but not in SR
  if ( m_proxyUtils.chfVeto(good_jets) ) return true;
  m_counter->increment(weight,incr++,"ChfCut");
  if ( m_proxyUtils.chfTileVeto(good_jets) ) return true;
  m_counter->increment(weight,incr++,"ChfCutForTile");

  if ( m_saveSeedEvents ) {
    write_xAOD_event();
    return true;
  }

  // ---------------------   Apply smearing & recompute MET


  // get Jet container with calibrated jets
  const xAOD::JetContainer* jets = 0;
  if ( !store->retrieve(jets, "SUSYJets").isSuccess() ) {
    throw std::runtime_error("Could not retrieve JetContainer with key SUSYJets");
  }

  // get MET recomputed after jet calibration
  const xAOD::MissingETContainer* metcontainer = 0;
  if( !store->retrieve( metcontainer, "MET_ZL" ).isSuccess() ) {
    throw std::runtime_error("Could not retrieve MissingETContainer with key MET_ZL");
  }
  const xAOD::MissingETComponentMap* metMap = 0;
  if( !event.retrieve(metMap,"METMap_RefFinal").isSuccess() ) {
    throw std::runtime_error("Could not retrieveMissingETComponentMap  with key METMap_ZL");
  }
  xAOD::MissingETContainer::const_iterator met_it = metcontainer->find("Final");
  //out() << " MET before smearing " <<  (*met_it)->mpx() << " " << (*met_it)->mpy() << std::endl;

  std::vector<SmearData> smrMc;
  m_smearingTool->DoSmearing(smrMc,*jets);


  // setup containers which will contain the smeared quantities
  xAOD::MissingETContainer* outCont = new xAOD::MissingETContainer() ;
  xAOD::MissingETAuxContainer* outAuxCont = new xAOD::MissingETAuxContainer();
  outCont->setStore(outAuxCont);
  if ( ! store->record(outCont,"MET_RefFinal_smeared").isSuccess() ) throw std::runtime_error("Could not register MissingETContainer with tag MET_RefFinal_smeared");
  if ( ! store->record(outAuxCont,"MET_AuxRefFinal_smeared").isSuccess() ) throw std::runtime_error("Could not register MissingETAuxContainer with tag MET_AuxRefFinal_smeared");

  const xAOD::JetContainer* originaljets = 0;
  if ( !event.retrieve(originaljets, m_jetkey).isSuccess() ) {
    throw std::runtime_error("Could not retrieve JetContainer with key "+m_jetkey);
  }
  std::pair< xAOD::JetContainer*, xAOD::ShallowAuxContainer* > smearedjets = xAOD::shallowCopyContainer( *originaljets );
  if ( ! store->record(smearedjets.first,"SUSYJets_smeared").isSuccess() ) {
    throw std::runtime_error("Could not store SUSYJets_smeared");
  }
  if ( ! store->record(smearedjets.second,"SUSYJets_smearedAux.").isSuccess() ) {
    throw std::runtime_error("Could not store SUSYJets_smearedAux.");
  }
  // fix links to containers to please METrebuilder
  if ( !xAOD::setOriginalObjectLink(*originaljets, *smearedjets.first) ) throw std::runtime_error("Could not set original links in jet container copy");

  if ( jets->size() != originaljets->size() ) throw std::logic_error("SUSY jets container should have the same size as original jet container");

  for ( std::size_t iexp = 0; iexp < smrMc.size(); ++iexp){
    xAOD::JetContainer::iterator targetit = smearedjets.first->begin();
    xAOD::JetContainer::const_iterator fromit = smrMc[iexp].jets.begin();
    for ( ; fromit != smrMc[iexp].jets.end(); ++fromit, ++targetit ) {
      (*targetit)->setJetP4((*fromit)->jetP4());
    }

    for(std::list<cafe::Processor*>::iterator it = m_processors.begin();
	it != m_processors.end();
	++it) {
      (*it)->incEventCount();
      if(!(*it)->processEvent(event)) return false;
    }

  }

  return true;
}


void ZeroLeptonDataDrivenQCD::GetResponseMaps(TFile* smearFnFile,TH2F*& nominal, TH2F*& tailHigh, TH2F*& tailLow)
{
  std::string histNameBase = "responseWE010_p1328_doublegaus_";
  std::string histName;
  TH2F* tempMap = 0;

  histName = histNameBase+std::string("optimal");
  tempMap = (TH2F*)smearFnFile->Get(histName.c_str());
  if ( !tempMap )  throw std::invalid_argument("ZeroLeptonDataDrivenQCD: could not find histogram "+histName);
  nominal = (TH2F*)tempMap->Clone("map");

  histName = histNameBase+std::string("high");
  tempMap = (TH2F*)smearFnFile->Get(histName.c_str());
  if ( !tempMap )  throw std::invalid_argument("ZeroLeptonDataDrivenQCD: could not find histogram "+histName);
  tailHigh = (TH2F*)tempMap->Clone("map");

  histName = histNameBase+std::string("low");
  tempMap = (TH2F*)smearFnFile->Get(histName.c_str());
  if ( !tempMap )  throw std::invalid_argument("ZeroLeptonDataDrivenQCD: could not find histogram "+histName);
  tailLow = (TH2F*)tempMap->Clone("map");
}

void ZeroLeptonDataDrivenQCD::InitialiseSmearing()
{

  std::string histNameBase = "responseWE010_p1328_doublegaus_";

  char* ROOTCOREDIR = std::getenv("ROOTCOREBIN");

  std::string filebveto = std::string(ROOTCOREDIR) + "/data/JetSmearing/July13/R_map2012_bveto_WE010_p1328_doublegaus_Responses_dphi.root";
  std::string filebtag = std::string(ROOTCOREDIR) + "/data/JetSmearing/July13/R_map2012_btag_WE010_p1328_doublegaus_Responses_dphi.root";

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
  m_smearingTool->AddSystematicResponseMaps(bveto_tailHigh,btag_tailHigh);
  m_smearingTool->AddSystematicResponseMaps(bveto_tailLow,btag_tailLow);
  unsigned int numsmear = 10;
  m_smearingTool->setProperty("NumberOfSmearedEvents",numsmear);
  m_smearingTool->setProperty("DoBJetSmearing",true);
  m_smearingTool->setProperty("UseTailWeights",true);
  m_smearingTool->setProperty("DoPhiSmearing",true);
  m_smearingTool->setProperty("DoGaussianCoreSmearing",true);

  if ( m_GaussianCoreSmearingType == "none" ) {
    m_smearingTool->setProperty("GaussianCoreSmearingType",SUSY::none);
  }
  else if ( m_GaussianCoreSmearingType == "optimal" ) {
    m_smearingTool->setProperty("GaussianCoreSmearingType",SUSY::optimal);
  }
  else if ( m_GaussianCoreSmearingType == "high" ) {
    m_smearingTool->setProperty("GaussianCoreSmearingType",SUSY::high);
  }
  else if ( m_GaussianCoreSmearingType == "low" ) {
    m_smearingTool->setProperty("GaussianCoreSmearingType",SUSY::low);
  }
  else  {
    throw std::logic_error("Could not recognize gaussing smearing type "+m_GaussianCoreSmearingType);
  }
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

ClassImp(ZeroLeptonDataDrivenQCD);

