
#include "ZeroLeptonRun2/BuildSUSYObjects.h"

#include "xAODRootAccess/tools/ReturnCheck.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TActiveStore.h"
#include "xAODRootAccess/TStore.h"
#include "xAODCore/ShallowCopy.h"
#include "xAODBase/IParticleHelpers.h"

#include "xAODJet/JetContainer.h"
#include "xAODJet/JetAuxContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODTau/TauJetContainer.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODMissingET/MissingETAuxContainer.h"
#include "PATInterfaces/SystematicSet.h"

#include "SUSYTools/SUSYObjDef_xAOD.h"
#include "cafe/Config.h"

#include "TVector2.h"

#include <stdexcept>

BuildSUSYObjects::BuildSUSYObjects(const char *name)
  : cafe::Processor(name),
    m_SUSYObjTool(0),
    m_IsData(false),
    m_IsAtlfast(false),
    m_UseSmearedJets(false),
    m_UseSystematics(false),
    m_PhotonInOR(false),
    m_jetkey(),
    m_suffix(),
    m_period(INVALID)
{
  cafe::Config config(name);
  m_IsData = config.get("IsData",false);
  m_IsAtlfast = config.get("IsAtlfast",false);
  m_jetkey = config.get("JetCollectionKey","xxxx");
  m_suffix = config.get("suffix","");
  m_UseSmearedJets = config.get("UseSmearedJets",false);
  m_UseSystematics = config.get("UseSystematics",false);
  if ( m_UseSmearedJets && m_UseSystematics ) throw std::logic_error("Cannot use jet smearing and systematics variations at the same time");
  m_PhotonInOR = config.get("PhotonInOR",false);
  m_period = periodFromString(config.get("Period","p13tev"));
  if ( m_period == p7tev ) throw(std::domain_error("BuildSUSYObjects does not support the 7tev run period"));


  m_SUSYObjTool = new ST::SUSYObjDef_xAOD("ZLST");
  m_SUSYObjTool->msg().setLevel( MSG::WARNING);
  m_SUSYObjTool->setProperty("IsData",(int)m_IsData);
  m_SUSYObjTool->setProperty("IsAtlfast",(int)m_IsAtlfast);
  m_SUSYObjTool->setProperty("METTauTerm","");
  if ( m_period == p8tev ) m_SUSYObjTool->setProperty("Is8TeV", true);
  else m_SUSYObjTool->setProperty("Is8TeV", false);

  if ( !m_SUSYObjTool->SUSYToolsInit().isSuccess() ) throw std::runtime_error("Could not initialise SUSYOBjDef ! ]SUSYToolsInit()]");

  if ( !m_SUSYObjTool->initialize().isSuccess() ) throw std::runtime_error("Could not initialise SUSYOBjDef !");
}

BuildSUSYObjects::~BuildSUSYObjects()
{
  if ( m_SUSYObjTool ) delete m_SUSYObjTool;
}

bool BuildSUSYObjects::processEvent(xAOD::TEvent& event)
{
  // Need to lump in all object treatment because of overlap removal and
  // xAOD Physics object forward declaration problems
  xAOD::TStore* store = xAOD::TActiveStore::store();

  if ( m_UseSystematics ) {
    m_SUSYObjTool->resetSystematics();
    // remove output from previous systematics
    store->remove("SUSYMET"+m_suffix).ignore();
    store->remove("SUSYJets").ignore();
    store->remove("SUSYJetsAux.").ignore();
    store->remove("SUSYMuons").ignore();
    store->remove("SUSYMuonsAux.").ignore();
    store->remove("SUSYElectrons").ignore();
    store->remove("SUSYElectronsAux.").ignore();
    store->remove("SUSYPhotons").ignore();
    store->remove("SUSYPhotonsAux.").ignore();
    //store->remove("SUSYTaus").ignore();
    //store->remove("SUSYTausAux.").ignore();
    store->remove("MET_MyRefFinal").ignore();
    store->remove("MET_MyRefFinalAux.").ignore();

    const CP::SystematicSet* currentSyst = 0;
    if ( !store->retrieve(currentSyst,"CurrentSystematicSet").isSuccess() ) throw std::runtime_error("Could not retrieve CurrentSystematicSet");

    if ( m_SUSYObjTool->applySystematicVariation(*currentSyst) != CP::SystematicCode::Ok) {
      out() << "Could not apply systematics " << currentSyst->name() << std::endl;
      return false; // stop processing in current algorithm list
    }

  }
  else if ( m_UseSmearedJets  ){
    // remove output from previous smearing iterations
    store->remove("SUSYMET"+m_suffix).ignore();
    store->remove("MET_MyRefFinal").ignore();
    store->remove("MET_MyRefFinalAux.").ignore();
  }
  //store->print();



  // Jets
  const xAOD::JetContainer* inputjets = 0;
  xAOD::JetContainer* outputjets = 0;
  if ( m_UseSmearedJets ) {
    if ( !store->retrieve(inputjets, "SUSYJets"+m_suffix).isSuccess() ) {
      throw std::runtime_error("Could not retrieve JetContainer with key SUSYJets"+m_suffix);
    }
    outputjets = const_cast<xAOD::JetContainer*>(inputjets);
  }
  else {
    if ( !event.retrieve(inputjets, m_jetkey).isSuccess() ) {
      throw std::runtime_error("Could not retrieve JetContainer with key "+m_jetkey);
    }
    std::pair< xAOD::JetContainer*, xAOD::ShallowAuxContainer* > susyjets = xAOD::shallowCopyContainer( *inputjets );
    if ( ! store->record(susyjets.first,"SUSYJets"+m_suffix).isSuccess() ) {
      throw std::runtime_error("Could not store SUSYJets"+m_suffix);
    }
    if ( ! store->record(susyjets.second,"SUSYJets"+m_suffix+"Aux.").isSuccess() ) {
      throw std::runtime_error("Could not store SUSYJets"+m_suffix+"Aux.");
    }
    outputjets = susyjets.first;
  }
  if ( !xAOD::setOriginalObjectLink(*inputjets, *outputjets) ) throw std::runtime_error("Could not set original links in jet container copy");

  // calibrate and fill properties
  xAOD::JetContainer::iterator jet_itr = outputjets->begin();
  xAOD::JetContainer::iterator jet_end = outputjets->end();
  for( ; jet_itr != jet_end; ++jet_itr ) {    
    if ( m_UseSmearedJets ) { 
      // no calibration since it was done in a first call to BuildSUSYObjects
      (*jet_itr)->auxdecor<char>("baseline") = 1;
    }
    else {
      if ( ! m_SUSYObjTool->FillJet(**jet_itr, 20000., 10.).isSuccess() ) throw std::runtime_error("Error in FillJet");    

    }
    if ( m_period == p8tev ) {
      m_SUSYObjTool->IsBJet(**jet_itr,false,1.85);
    }
    else if ( m_period == p13tev ) {
      m_SUSYObjTool->IsBJet(**jet_itr,true,0.7892);
    }

    //out() << "pt " << (*jet_itr)->pt() << " baseline " << (int)((*jet_itr)->auxdecor<char>("baseline")) << " bad " << (int)((*jet_itr)->auxdecor<char>("bad")) << " bjet " << (int)((*jet_itr)->auxdecor<char>("bjet")) <<  " container " << (*jet_itr)->container() << std::endl;


  }

  // Muons
  std::pair< xAOD::MuonContainer*, xAOD::ShallowAuxContainer* > susymuons = std::make_pair<xAOD::MuonContainer*, xAOD::ShallowAuxContainer* >(NULL,NULL);
  if ( m_UseSmearedJets ) {
    if ( !store->retrieve(susymuons.first ,"SUSYMuons").isSuccess() ) {
      throw std::runtime_error("Could not retrieve MuonContainer with key SUSYMuons");
    }
    
  }
  else {
    const xAOD::MuonContainer* muons = 0;
    if ( !event.retrieve( muons,"Muons").isSuccess() ) {
      throw std::runtime_error("Could not retrieve MuonContainer with key Muons");
    }
    susymuons = xAOD::shallowCopyContainer( *muons );

    // calibrate and fill properties
    xAOD::MuonContainer::iterator mu_itr = (susymuons.first)->begin();
    xAOD::MuonContainer::iterator mu_end = (susymuons.first)->end();
    if ( !xAOD::setOriginalObjectLink(*muons, *susymuons.first) ) throw std::runtime_error("Could not set original links in Muon container copy");
    for( ; mu_itr != mu_end; ++mu_itr ) {
      if ( ! m_SUSYObjTool->FillMuon(**mu_itr).isSuccess() ) throw std::runtime_error("Error in FillMuon");
      m_SUSYObjTool->IsSignalMuon(**mu_itr);
      m_SUSYObjTool->IsCosmicMuon(**mu_itr);

      // kill non baseline muon by setting 4-vector to small value
      if ( ((*mu_itr)->muonType() != xAOD::Muon::Combined &&
	   (*mu_itr)->muonType() != xAOD::Muon::MuonStandAlone && 
	    (*mu_itr)->muonType() != xAOD::Muon::SegmentTagged) ||
	   (*mu_itr)->auxdecor<char>("baseline") != 1)
      {
	(*mu_itr)->setP4(1.,(*mu_itr)->eta(),(*mu_itr)->phi()); 
      }
      //out() << " Muon " << (*mu_itr)->pt() << " " << (*mu_itr)->eta()
      //    << " " << (*mu_itr)->phi() << std::endl;
    }
    
    if ( ! store->record(susymuons.first,"SUSYMuons"+m_suffix).isSuccess() ) {
      throw std::runtime_error("Could not store SUSYMuons"+m_suffix);
    }
    if ( ! store->record(susymuons.second,"SUSYMuons"+m_suffix+"Aux.").isSuccess()) {
      throw std::runtime_error("Could not store SUSYMuons"+m_suffix+"Aux.");
    }
  }

  // Electrons
  std::pair< xAOD::ElectronContainer*, xAOD::ShallowAuxContainer* > susyelectrons = std::make_pair< xAOD::ElectronContainer*, xAOD::ShallowAuxContainer* >(NULL,NULL);
  if ( m_UseSmearedJets ) {
    if ( !store->retrieve(susyelectrons.first, "SUSYElectrons").isSuccess() ){
      throw std::runtime_error("Could not retrieve ElectronContainer with key SUSYElectrons");
    }
  }
  else {
    const xAOD::ElectronContainer* electrons = 0;
    if ( !event.retrieve(electrons, "ElectronCollection").isSuccess() ){
      throw std::runtime_error("Could not retrieve ElectronContainer with key ElectronCollection");
    }
    susyelectrons = xAOD::shallowCopyContainer(*electrons);
    
    // calibrate and fill properties
    xAOD::ElectronContainer::iterator el_itr = susyelectrons.first->begin();
    xAOD::ElectronContainer::iterator el_end = susyelectrons.first->end();
    if ( !xAOD::setOriginalObjectLink(*electrons, *susyelectrons.first) ) throw std::runtime_error("Could not set original links in electron container copy");
    for( ; el_itr != el_end; ++el_itr ) {
      if ( ! m_SUSYObjTool->FillElectron(**el_itr,10000.,2.47).isSuccess() ) throw std::runtime_error("Error in FillElectron");
      m_SUSYObjTool->IsSignalElectron(**el_itr);
    }

    if ( ! store->record(susyelectrons.first,"SUSYElectrons"+m_suffix).isSuccess() ) {
      throw std::runtime_error("Could not store SUSYElectrons"+m_suffix);
    }
    if ( ! store->record(susyelectrons.second,"SUSYElectrons"+m_suffix+"Aux.").isSuccess()) {
      throw std::runtime_error("Could not store SUSYElectrons"+m_suffix+"Aux.");
    }
  }

  // Photons
  std::pair< xAOD::PhotonContainer*, xAOD::ShallowAuxContainer* > susyphotons = std::make_pair< xAOD::PhotonContainer*, xAOD::ShallowAuxContainer* >(NULL,NULL);
  if ( m_UseSmearedJets ) {
    if ( !store->retrieve(susyphotons.first,"SUSYPhotons").isSuccess() ){
      throw std::runtime_error("Could not retrieve PhotonContainer with key SUSYPhotons");
    }
  }
  else {
    const xAOD::PhotonContainer* photons = 0;
    if ( !event.retrieve(photons,"PhotonCollection").isSuccess() ){
      throw std::runtime_error("Could not retrieve PhotonContainer with key PhotonCollection");
    }

    susyphotons = xAOD::shallowCopyContainer(*photons);
    
    // calibrate and fill properties
    xAOD::PhotonContainer::iterator ph_itr = susyphotons.first->begin();
    xAOD::PhotonContainer::iterator ph_end = susyphotons.first->end();
    if ( !xAOD::setOriginalObjectLink(*photons, *susyphotons.first) ) throw std::runtime_error("Could not set original links in photon container copy");
        
    for( ; ph_itr != ph_end; ++ph_itr ) {
      if ( ! m_SUSYObjTool->FillPhoton(**ph_itr).isSuccess() ) throw std::runtime_error("Error in FillPhoton");
    }

    if ( ! store->record(susyphotons.first,"SUSYPhotons"+m_suffix).isSuccess() ) {
      throw std::runtime_error("Could not store SUSYPhotons"+m_suffix);
    }
    if ( ! store->record(susyphotons.second,"SUSYPhotons"+m_suffix+"Aux.").isSuccess()) {
      throw std::runtime_error("Could not store SUSYPhotons"+m_suffix+"Aux.");
    }
  }

  /*
  // Taus
  std::pair< xAOD::TauJetContainer*, xAOD::ShallowAuxContainer* > susytaus = std::make_pair< xAOD::TauJetContainer*, xAOD::ShallowAuxContainer* >(NULL,NULL);
  if ( m_UseSmearedJets ) {
    if (!store->retrieve(susytaus.first,"SUSYTaus").isSuccess()){
      throw std::runtime_error("Could not retrieve TauJetContainer with key SUSYTaus");
    }
  }
  else {
    const xAOD::TauJetContainer* taus = 0;
    if (!event.retrieve(taus,"TauRecContainer").isSuccess()){
      throw std::runtime_error("Could not retrieve TauJetContainer with key TauRecContainer");
    }

    susytaus = xAOD::shallowCopyContainer(*taus);
    xAOD::TauJetContainer::iterator tau_itr = susytaus.first->begin();
    xAOD::TauJetContainer::iterator tau_end = susytaus.first->end();
    for( ; tau_itr != tau_end; ++tau_itr ) {
      if ( ! m_SUSYObjTool->FillTau( **tau_itr).isSuccess() ) throw std::runtime_error("Error in FillTau");
    }
    
    if ( ! store->record(susytaus.first,"SUSYTaus"+m_suffix).isSuccess() ) {
      throw std::runtime_error("Could not store SUSYTaus"+m_suffix);
    }
    if ( ! store->record(susytaus.second,"SUSYTaus"+m_suffix+"Aux.").isSuccess()) {
      throw std::runtime_error("Could not store SUSYTaus"+m_suffix+"Aux.");
    }
  }
  */

  // Overlap removal
  if ( m_PhotonInOR ) {
    if ( ! m_SUSYObjTool->OverlapRemoval(susyelectrons.first, susymuons.first, outputjets, susyphotons.first, false, 0.2, 0.4, 0.4, 0.01, 0.05, 0.2, 0.4).isSuccess() ) throw std::runtime_error("Error in OverlapRemoval");
  }
  else {
    if ( ! m_SUSYObjTool->OverlapRemoval(susyelectrons.first, susymuons.first, outputjets, false, 0.2, 0.4, 0.4, 0.01, 0.05).isSuccess() ) throw std::runtime_error("Error in OverlapRemoval");
  }

  // Missing ET
  TVector2* MissingET = new TVector2(0.,0.);
  xAOD::MissingETContainer* rebuiltmetc = new xAOD::MissingETContainer();
  xAOD::MissingETAuxContainer* rebuiltmetcAux = new xAOD::MissingETAuxContainer();
  rebuiltmetc->setStore(rebuiltmetcAux);
  if ( ! store->record(rebuiltmetc,"MET_MyRefFinal"+m_suffix).isSuccess() ) {
    throw std::runtime_error("Unable to store MissingETContainer with tag MET_MyRefFinal"+m_suffix);
  }
  if ( ! store->record(rebuiltmetcAux,"MET_MyRefFinal"+m_suffix+"Aux.").isSuccess() ) {
    throw std::runtime_error("Unable to store MissingETAuxContainer with tag MET_MyRefFinal"+m_suffix+"Aux");
  }

  if ( ! m_SUSYObjTool->GetMET(*rebuiltmetc,
			       outputjets,
			       susyelectrons.first,
			       susymuons.first,
			       susyphotons.first,
			       0,
			       false).isSuccess() 
       ) throw std::runtime_error("Error in GetMET");
  xAOD::MissingETContainer::const_iterator met_it = rebuiltmetc->find("Final");
  if ( met_it == rebuiltmetc->end() ) throw std::runtime_error("Could not find Final MET after running  GetMET");
  MissingET->Set((*met_it)->mpx(), (*met_it)->mpy());
  if ( ! store->record(MissingET,"SUSYMET"+m_suffix).isSuccess() ) {
    throw std::runtime_error("Could not store SUSYMET");
  }

  //out() << " MET after smearing " <<  MissingET->X() << " " << MissingET->Y()  << std::endl;


  return true;
}


ClassImp(BuildSUSYObjects);

