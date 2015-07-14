
#include "ZeroLeptonRun2/BuildSUSYObjects.h"
#include "ZeroLeptonRun2/ZeroLeptonUtils.h"

#include "xAODRootAccess/tools/ReturnCheck.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TActiveStore.h"
#include "xAODRootAccess/TStore.h"
#include "xAODCore/ShallowCopy.h"
#include "xAODBase/IParticleHelpers.h"
#include "AthContainers/AuxElement.h"

#include "xAODEventInfo/EventInfo.h"
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
    m_taukey(),
    m_suffix(),
    m_period(INVALID),
    m_derivationTag(INVALID_Derivation),
    m_JESNuisanceParameterSet(0),
    m_ECKey(""),
    m_PCKey("")
{
  cafe::Config config(name);
  m_IsData = config.get("IsData",false);
  m_IsAtlfast = config.get("IsAtlfast",false);
  m_jetkey = config.get("JetContainerKey","xxxx");
  m_taukey = config.get("TauContainerKey","xxxx");
  m_suffix = config.get("suffix","");
  m_ECKey = config.get("ElectronContainerKey","ElectronCollection");
  m_PCKey = config.get("PhotonContainerKey","PhotonCollection");
  m_UseSmearedJets = config.get("UseSmearedJets",false);
  m_UseSystematics = config.get("UseSystematics",false);
  if ( m_UseSmearedJets && m_UseSystematics ) throw std::logic_error("Cannot use jet smearing and systematics variations at the same time");
  m_PhotonInOR = config.get("PhotonInOR",false);
  m_period = periodFromString(config.get("Period","p13tev"));
  if ( m_period == p7tev ) throw(std::domain_error("BuildSUSYObjects does not support the 7tev run period"));
  if ( m_period == p8tev ) throw(std::domain_error("Due to interface changes in SUSYTools BuildSUSYObjects no longer supports the 8tev run period"));

  m_derivationTag = derivationTagFromString(config.get("DerivationTag",""));
  if ( m_derivationTag == INVALID_Derivation ) throw(std::domain_error("ZeroLeptonSR: invalid derivation tag specified"));

  m_JESNuisanceParameterSet = config.get("JESNuisanceParameterSet",0);

}

BuildSUSYObjects::~BuildSUSYObjects()
{
  //FIXME crash in SUSYObjDef_xAOD destructor
  //if ( m_SUSYObjTool ) delete m_SUSYObjTool;
}

void BuildSUSYObjects::initSUSYTools()
{
  // SUSYObjDef_xAOD needs a TEvent to be defined, moved from constructor to begin()
  m_SUSYObjTool = new ST::SUSYObjDef_xAOD( m_PhotonInOR ? "ZLST_GAMMA" : "ZLST");
  m_SUSYObjTool->msg().setLevel( MSG::WARNING);
  ST::SettingDataSource datasource = m_IsData ? ST::Data : (m_IsAtlfast ? ST::AtlfastII : ST::FullSim);
  m_SUSYObjTool->setProperty("DataSource",datasource).ignore();
  m_SUSYObjTool->setProperty("METTauTerm","").ignore();
  if ( m_period == p8tev ) m_SUSYObjTool->setProperty("Is8TeV", true).ignore();
  else m_SUSYObjTool->setProperty("Is8TeV", false).ignore();

  xAOD::JetInput::Type jetType =  ZeroLeptonUtils::JetTypeFromString(m_jetkey);
  if ( jetType == xAOD::JetInput::Uncategorized ) throw(std::domain_error("ZeroLeptonSR: could not identify JetType"));
  m_SUSYObjTool->setProperty("JetInputType", jetType).ignore();

  m_SUSYObjTool->setProperty("JESNuisanceParameterSet",m_JESNuisanceParameterSet).ignore();
  m_SUSYObjTool->setProperty("DoJetAreaCalib",true).ignore();
  m_SUSYObjTool->setProperty("DoJetGSCCalib",true).ignore();
  if ( m_derivationTag == p1872 ) {
    m_SUSYObjTool->setProperty("METInputCont","MET_RefFinalFix").ignore();
    m_SUSYObjTool->setProperty("METInputMap","METMap_RefFinalFix").ignore();
  }
  //m_SUSYObjTool->setProperty("IsoWP","Gradient").ignore();

  // set our own tau selection
  TauAnalysisTools::TauSelectionTool* tauSelTool;
  TauAnalysisTools::TauEfficiencyCorrectionsTool* tauEffTool;
  tauSelTool = new TauAnalysisTools::TauSelectionTool( m_PhotonInOR ? "TauSelectionTool_GAMMA" : "TauSelectionTool");
  tauSelTool->msg().setLevel( MSG::WARNING);
  //tauSelTool->msg().setLevel( MSG::VERBOSE);
  tauSelTool->setProperty("PtMin", 20. ).ignore(); // pt in GeV
  std::vector<double> vAbsEtaRegion = {0, 1.37, 1.52, 2.5};
  tauSelTool->setProperty("AbsEtaRegion", vAbsEtaRegion).ignore();
  tauSelTool->setProperty("AbsCharge", 1.).ignore();
  std::vector<size_t> vNTracks = {1, 3};
  tauSelTool->setProperty( "NTracks", vNTracks).ignore();
  tauSelTool->setProperty( "JetIDWP",  int(TauAnalysisTools::JETIDBDTMEDIUM)).ignore();

  tauSelTool->setProperty(
    "SelectionCuts",
    (int) TauAnalysisTools::SelectionCuts(TauAnalysisTools::CutPt |
					  TauAnalysisTools::CutAbsEta    |
					  TauAnalysisTools::CutAbsCharge |
					  TauAnalysisTools::CutNTrack    |
					  TauAnalysisTools::CutJetIDWP )
			  ).ignore();

  tauSelTool->initialize().ignore();
  m_tauSelTool = tauSelTool;

  tauEffTool = new TauAnalysisTools::TauEfficiencyCorrectionsTool("TauEfficiencyCorrectionsTool",tauSelTool);
  tauEffTool->msg().setLevel( MSG::WARNING);
  //tauEffTool->msg().setLevel( MSG::VERBOSE);
  tauEffTool->setProperty("EfficiencyCorrectionType", (int)TauAnalysisTools::SFJetID).ignore();
  tauEffTool->setProperty("SysDirection", 1).ignore();
  tauEffTool->initialize().ignore();
  m_tauEffTool = tauEffTool;

  m_SUSYObjTool->setProperty("TauSelectionTool",m_tauSelTool).ignore();
  m_SUSYObjTool->setProperty("TauEfficiencyCorrectionsTool",m_tauEffTool).ignore();

  const CP::SystematicSet& recommendedSystematics = m_tauEffTool->recommendedSystematics();
  std::vector<CP::SystematicSet> systSetList;
  systSetList.reserve(recommendedSystematics.size()*2); // allow for continuous systematics
  systSetList.push_back(CP::SystematicSet());
  for(const auto& syst : recommendedSystematics){
    systSetList.push_back(CP::SystematicSet());
    systSetList.back().insert(syst);
  }
  m_tauEffSystSetList = systSetList;


  if ( !m_SUSYObjTool->SUSYToolsInit().isSuccess() ) throw std::runtime_error("Could not initialise SUSYOBjDef ! ]SUSYToolsInit()]");

  if ( !m_SUSYObjTool->initialize().isSuccess() ) throw std::runtime_error("Could not initialise SUSYOBjDef !");
}


bool BuildSUSYObjects::processEvent(xAOD::TEvent& event)
{

  // SUSYTools initialisation must be delayed until we have a TEvent associated
  // with a file due to xAODConfigTool 
  if ( !m_SUSYObjTool ) initSUSYTools();

  // Need to lump in all object treatment because of overlap removal and
  // xAOD Physics object forward declaration problems
  xAOD::TStore* store = xAOD::TActiveStore::store();

  if ( m_UseSystematics ) {
    m_SUSYObjTool->resetSystematics();
    // remove output from previous systematics
    store->remove("SUSYMET"+m_suffix).ignore();
    store->remove("SUSYJets"+m_suffix).ignore();
    store->remove("SUSYJets"+m_suffix+"Aux.").ignore();
    store->remove("SUSYMuons"+m_suffix).ignore();
    store->remove("SUSYMuons"+m_suffix+"Aux.").ignore();
    store->remove("SUSYElectrons"+m_suffix).ignore();
    store->remove("SUSYElectrons"+m_suffix+"Aux.").ignore();
    store->remove("SUSYPhotons"+m_suffix).ignore();
    store->remove("SUSYPhotons"+m_suffix+"Aux.").ignore();
    store->remove("SUSYTaus").ignore();
    store->remove("SUSYTausAux.").ignore();
    store->remove("MET_ZL"+m_suffix).ignore();
    store->remove("MET_ZL"+m_suffix+"Aux.").ignore();

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
    store->remove("MET_ZL").ignore();
    store->remove("MET_ZLAux.").ignore();
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
      if ( ! m_SUSYObjTool->FillJet(**jet_itr, true).isSuccess() ) throw std::runtime_error("Error in FillJet");    
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
      m_SUSYObjTool->IsBadMuon(**mu_itr);

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
    
      float muSF = 0;
      muSF = (float) m_SUSYObjTool->GetTotalMuonSF(*susymuons.first);
      //float testSF = (float) m_SUSYObjTool->GetSignalMuonSF(**mu_itr);
      //std::cout << "MUON WEIGHT : " << muSF << "  " << (*mu_itr)->pt() << "  " << testSF << std::endl;
      (*mu_itr)->auxdecor<float>("sf") = muSF ; 
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
    if ( !event.retrieve(electrons, m_ECKey).isSuccess() ){
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

      float elSF=0;
      elSF = (float) m_SUSYObjTool->GetTotalElectronSF(*susyelectrons.first);
      (*el_itr)->auxdecor<float>("sf") = elSF ; 
      //float testSF = (float) m_SUSYObjTool->GetSignalElecSF(**el_itr);
      //std::cout << "ELECTRON WEIGHT : " << elSF << "  " << (*el_itr)->pt() << "  " << testSF << std::endl;
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
    if ( !event.retrieve(photons,m_PCKey).isSuccess() ){
      throw std::runtime_error("Could not retrieve PhotonContainer with key PhotonCollection");
    }

    susyphotons = xAOD::shallowCopyContainer(*photons);
    
    // calibrate and fill properties
    xAOD::PhotonContainer::iterator ph_itr = susyphotons.first->begin();
    xAOD::PhotonContainer::iterator ph_end = susyphotons.first->end();
    if ( !xAOD::setOriginalObjectLink(*photons, *susyphotons.first) ) throw std::runtime_error("Could not set original links in photon container copy");
        
    for( ; ph_itr != ph_end; ++ph_itr ) {
      if ( ! m_SUSYObjTool->FillPhoton(**ph_itr).isSuccess() ) throw std::runtime_error("Error in FillPhoton");
      m_SUSYObjTool->IsSignalPhoton(**ph_itr,25000.);
      //if ((*ph_itr)->pt()>10000.) out() << "Photon : pt " << (*ph_itr)->pt() << " baseline " << (int)((*ph_itr)->auxdecor<char>("baseline")) << " signal " << (int)((*ph_itr)->auxdecor<char>("signal")) << std::endl;
    }

    if ( ! store->record(susyphotons.first,"SUSYPhotons"+m_suffix).isSuccess() ) {
      throw std::runtime_error("Could not store SUSYPhotons"+m_suffix);
    }
    if ( ! store->record(susyphotons.second,"SUSYPhotons"+m_suffix+"Aux.").isSuccess()) {
      throw std::runtime_error("Could not store SUSYPhotons"+m_suffix+"Aux.");
    }
  }

  // Taus
  std::pair< xAOD::TauJetContainer*, xAOD::ShallowAuxContainer* > susytaus = std::make_pair< xAOD::TauJetContainer*, xAOD::ShallowAuxContainer* >(NULL,NULL);
  if ( m_UseSmearedJets ) {
    if (!store->retrieve(susytaus.first,"SUSYTaus").isSuccess()){
      throw std::runtime_error("Could not retrieve TauJetContainer with key SUSYTaus");
    }
  }
  else {
    const xAOD::TauJetContainer* taus = 0;
    if (!event.retrieve(taus,m_taukey).isSuccess()){
      throw std::runtime_error("Could not retrieve TauJetContainer with key TauRecContainer");
    }
    
    xAOD::TauJet::Decorator<float> dec_SFJetID("SFJetID");
    xAOD::TauJet::Decorator<float> dec_SFJetIDStatUp("SFJetIDStatUp");
    xAOD::TauJet::Decorator<float> dec_SFJetIDStatDown("SFJetIDStatDown");
    xAOD::TauJet::Decorator<float> dec_SFJetIDSystUp("SFJetIDSystUp");
    xAOD::TauJet::Decorator<float> dec_SFJetIDSystDown("SFJetIDSystDown");
    susytaus = xAOD::shallowCopyContainer(*taus);
    xAOD::TauJetContainer::iterator tau_itr = susytaus.first->begin();
    xAOD::TauJetContainer::iterator tau_end = susytaus.first->end();
    for( ; tau_itr != tau_end; ++tau_itr ) {
      
      if ( ! m_SUSYObjTool->FillTau( **tau_itr).isSuccess() ) throw std::runtime_error("Error in FillTau");
      if( (*tau_itr)->auxdecor<char>("baseline")==1 ){
	for ( const auto& syst : m_tauEffSystSetList ){
	  // one by one apply systematic variation 
	  if (m_tauEffTool->applySystematicVariation(syst) != CP::SystematicCode::Ok){
	    throw std::runtime_error("Could not configure for systematic variatoin" );
	  }else{
	    m_tauEffTool->applyEfficiencyScaleFactor( **tau_itr );
	    std::string systName = syst.name();
	    float sf = (float)( *tau_itr )->auxdata< double >("TauScaleFactorJetID");
	    //std::cout<<"TauScaleFactorJetID syst:"<<systName<<" "<<sf<<std::endl;
	    if( systName == "" )                               dec_SFJetID( **tau_itr )         = sf;
	    else if( systName == "TAUS_EFF_JETID_STAT__1up"   ) dec_SFJetIDStatUp( **tau_itr )   = sf;
	    else if( systName == "TAUS_EFF_JETID_STAT__1down" ) dec_SFJetIDStatDown( **tau_itr ) = sf;
	    else if( systName == "TAUS_EFF_JETID_SYST__1up"   ) dec_SFJetIDSystUp( **tau_itr )    = sf;
	    else if( systName == "TAUS_EFF_JETID_SYST__1down" ) dec_SFJetIDSystDown( **tau_itr )  = sf;
	  }
	  CP::SystematicSet defaultSet;
	  if(m_tauEffTool->applySystematicVariation(defaultSet) != CP::SystematicCode::Ok){
	    throw std::runtime_error("Could not configure TauEfficiencyCorrectionsTool for default systematic setting");
	  }
	}
      }
    }
    
    if ( ! store->record(susytaus.first,"SUSYTaus"+m_suffix).isSuccess() ) {
      throw std::runtime_error("Could not store SUSYTaus"+m_suffix);
    }
    if ( ! store->record(susytaus.second,"SUSYTaus"+m_suffix+"Aux.").isSuccess()) {
      throw std::runtime_error("Could not store SUSYTaus"+m_suffix+"Aux.");
    }
  }


  // Overlap removal
  if ( m_PhotonInOR ) {
    if ( ! m_SUSYObjTool->OverlapRemoval(susyelectrons.first, susymuons.first, outputjets, susyphotons.first, false, false, false, 0.2, 0.4, 0.4, 0.01, 0.05, 0.4, 0.4, 0.4).isSuccess() ) throw std::runtime_error("Error in OverlapRemoval");
  }
  else {
    if ( ! m_SUSYObjTool->OverlapRemoval(susyelectrons.first, susymuons.first, outputjets, false, false, false, 0.2, 0.4, 0.4, 0.01, 0.05).isSuccess() ) throw std::runtime_error("Error in OverlapRemoval");
  }

  // signal and btag jet now depend on OR, so loop again on jets
  for ( const auto& jet_itr : *outputjets ) {
    if ( !m_UseSmearedJets ) { 
      m_SUSYObjTool->IsSignalJet(*jet_itr, 20000.,10., -1e+99); // no JVT cut
    }
    if ( m_period == p8tev ) {
      //FIXME m_SUSYObjTool->IsBJet(**jet_itr,false,1.85);
    }
    else if ( m_period == p13tev ) {
      m_SUSYObjTool->IsBJet(*jet_itr);
    }

  }


  // Missing ET
  TVector2* MissingET = new TVector2(0.,0.);
  xAOD::MissingETContainer* rebuiltmetc = new xAOD::MissingETContainer();
  xAOD::MissingETAuxContainer* rebuiltmetcAux = new xAOD::MissingETAuxContainer();
  rebuiltmetc->setStore(rebuiltmetcAux);
  if ( ! store->record(rebuiltmetc,"MET_ZL"+m_suffix).isSuccess() ) {
    throw std::runtime_error("Unable to store MissingETContainer with tag MET_ZL"+m_suffix);
  }
  if ( ! store->record(rebuiltmetcAux,"MET_ZL"+m_suffix+"Aux.").isSuccess() ) {
    throw std::runtime_error("Unable to store MissingETAuxContainer with tag MET_ZL"+m_suffix+"Aux");
  }

  if ( ! m_SUSYObjTool->GetMET(*rebuiltmetc,
			       outputjets,
			       susyelectrons.first,
			       susymuons.first,
			       susyphotons.first,
			       0,
			       false,
			       false,
			       0).isSuccess() 
       ) throw std::runtime_error("Error in GetMET");
  xAOD::MissingETContainer::const_iterator met_it = rebuiltmetc->find("Final");
  if ( met_it == rebuiltmetc->end() ) throw std::runtime_error("Could not find Final MET after running  GetMET");
  MissingET->Set((*met_it)->mpx(), (*met_it)->mpy());
  if ( ! store->record(MissingET,"SUSYMET"+m_suffix).isSuccess() ) {
    throw std::runtime_error("Could not store SUSYMET");
  }

  //out() << " MET after smearing " <<  MissingET->X() << " " << MissingET->Y()  << std::endl;

  if ( m_suffix == "" ) fillTriggerInfo(event);

  return true;
}

void BuildSUSYObjects::fillTriggerInfo(xAOD::TEvent& event) const
{
#ifndef ZLDC14
  static std::vector<std::string> trigNames = {
    "L1_XE50",
    "L1_XE70",
    "HLT_xe70",
    "HLT_xe70_pueta",
    "HLT_xe100",
    "HLT_xe100_pueta",
    "HLT_e28_tight_iloose",
    "HLT_e60_medium",
    "HLT_mu26_imedium",
    "HLT_mu50",
    "HLT_j30_xe10_razor170",
    "HLT_xe70_tc_em",
    "HLT_xe70_tc_lcw",
    "HLT_xe70_mht",
    "HLT_xe70_pufit",
    "HLT_xe100_tc_em",
    "HLT_xe100_tc_lcw",
    "HLT_xe100_mht",
    "HLT_xe100_pufit",
    "HLT_3j175",
    "HLT_4j85",
    "HLT_5j85",
    "HLT_6j25",
    "HLT_6j45_0eta240",
    "HLT_6j55_0eta240_L14J20",
    "HLT_7j45",
    "L1_2J15",
    "HLT_2j55_bloose",
    "HLT_j80_xe80",
    "HLT_e24_lhmedium_iloose_L1EM18VH",
    "HLT_e60_lhmedium",
    "HLT_mu20_iloose_L1MU15",
    "HLT_mu40",
    "HLT_mu50",
    "HLT_g120_loose",
    "HLT_g120_lhloose",
    "HLT_mu18",
    "HLT_e17_lhloose_L1EM15",
    "HLT_e17_loose_L1EM15",
    "HLT_mu14_iloose"    
  };

  const xAOD::EventInfo* eventInfo = 0;
  if ( ! event.retrieve( eventInfo, "EventInfo").isSuccess() ) throw std::runtime_error("BuildSUSYObjects: Could not retrieve EventInfo");

  std::bitset<32> triggers;
  for ( size_t i = 0; i < trigNames.size(); i++ ) {
    char pass =  m_SUSYObjTool->IsTrigPassed(trigNames[i]);
    triggers[i] = pass;
    eventInfo->auxdecor<char>(trigNames[i]) = pass;
  }

  xAOD::TStore* store = xAOD::TActiveStore::store();
  unsigned long* triggerSet = new unsigned long;
  *triggerSet =  triggers.to_ulong();
  if ( ! store->record(triggerSet,"triggerbits").isSuccess() ) {
    throw std::runtime_error("Could not store trigger bits");
  }
  //std::cout << " trigger " << triggers << " " << *triggerSet << std::endl;
#endif
}



ClassImp(BuildSUSYObjects);

