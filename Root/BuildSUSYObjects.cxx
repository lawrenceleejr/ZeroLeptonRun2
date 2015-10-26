
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
#include "TRegexp.h"

#include <stdexcept>
#include <vector>

BuildSUSYObjects::BuildSUSYObjects(const char *name)
  : cafe::Processor(name),
    m_SUSYObjTool(0),
    m_IsData(false),
    m_IsAtlfast(false),
    m_UseSmearedJets(false),
    m_DoSystematics(false),
    m_PhotonInOR(false),
    m_jetkey(),
    m_taukey(),
    m_suffix(),
    m_period(INVALID),
    m_derivationTag(INVALID_Derivation),
    m_JESNuisanceParameterSet(0),
    m_ECKey(""),
    m_PCKey(""),
    m_SystInfoList(),
    m_SystMatch()
{
  cafe::Config config(name);
  m_IsData = config.get("IsData",false);
  m_Is25ns = config.get("Is25ns",true);
  m_IsAtlfast = config.get("IsAtlfast",false);
  m_jetkey = config.get("JetContainerKey","xxxx");
  m_taukey = config.get("TauContainerKey","xxxx");
  m_suffix = config.get("suffix","");
  m_ECKey = config.get("ElectronContainerKey","ElectronCollection");
  m_PCKey = config.get("PhotonContainerKey","PhotonCollection");
  m_UseSmearedJets = config.get("UseSmearedJets",false);
  m_DoSystematics = config.get("DoSystematics",false);
  if ( m_UseSmearedJets && m_DoSystematics ) throw std::logic_error("Cannot use jet smearing and systematics variations at the same time");
  m_PhotonInOR = config.get("PhotonInOR",false);
  m_period = periodFromString(config.get("Period","p13tev"));
  if ( m_period == p7tev ) throw(std::domain_error("BuildSUSYObjects does not support the 7tev run period"));
  if ( m_period == p8tev ) throw(std::domain_error("Due to interface changes in SUSYTools BuildSUSYObjects no longer supports the 8tev run period"));

  m_derivationTag = derivationTagFromString(config.get("DerivationTag",""));
  if ( m_derivationTag == INVALID_Derivation ) throw(std::domain_error("ZeroLeptonSR: invalid derivation tag specified"));

  m_JESNuisanceParameterSet = config.get("JESNuisanceParameterSet",0);

  m_SystMatch = config.getVString("SystMatch");
}

BuildSUSYObjects::~BuildSUSYObjects()
{
  //FIXME crash in SUSYObjDef_xAOD destructor
  //if ( m_SUSYObjTool ) delete m_SUSYObjTool;
}

void  BuildSUSYObjects::inputFileOpened(TFile *file)
{
  // SUSYTools initialisation must be delayed until we have a TEvent associated
  // with a file due to xAODConfigTool 
  if ( !m_SUSYObjTool ) {
    initSUSYTools();
  }
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
  if ( m_PhotonInOR ) { m_SUSYObjTool->setProperty("DoPhotonOR",true); }
  m_SUSYObjTool->setProperty("PhotonIsoWP","Cone20").ignore();
  m_SUSYObjTool->setProperty("PhotonId","Loose").ignore();

  // set our own tau selection

  // truth tau matching is needed for FillTau
  // https://svnweb.cern.ch/trac/atlasoff/browser/PhysicsAnalysis/TauID/TauAnalysisTools/trunk/doc/README-TauTruthMatchingTool.rst
  if (  !m_IsData && m_tauTruthMatchTool.empty() ) {
    TauAnalysisTools::TauTruthMatchingTool* tauTruthMatchTool = new TauAnalysisTools::TauTruthMatchingTool("TauTruthMatchingTool");
    tauTruthMatchTool->setProperty("WriteTruthTaus",true);
    tauTruthMatchTool->initialize().ignore();
    m_tauTruthMatchTool = tauTruthMatchTool;
  }

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

  tauEffTool = new TauAnalysisTools::TauEfficiencyCorrectionsTool("TauEfficiencyCorrectionsTool"+m_suffix,tauSelTool);
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


  //if ( !m_SUSYObjTool->SUSYToolsInit().isSuccess() ) throw std::runtime_error("Could not initialise SUSYOBjDef ! ]SUSYToolsInit()]");

  if ( !m_SUSYObjTool->initialize().isSuccess() ) throw std::runtime_error("Could not initialise SUSYOBjDef !");

  if ( m_DoSystematics ) {
    std::vector<ST::SystInfo> sysInfos = m_SUSYObjTool->getSystInfoList();
    if (  m_SystMatch.empty() ) {
      m_SystInfoList = sysInfos;
    }
    else {
      for ( const auto & sys : sysInfos ) {
	const CP::SystematicSet& systSet = sys.systset;
	std::string name = systSet.name(); 
	bool matched = false;
	if ( name == "" ) {
	  matched = true;
	}
	else {
	  // check if name matches wildcard expression
	  for ( const auto& wildcardexp : m_SystMatch ){
	    TRegexp re = TRegexp(wildcardexp.c_str(),kTRUE);
	    Ssiz_t l;
	    if ( re.Index(name,&l) >= 0 ) {
	      matched = true;
	      break;
	    }
	  }
	}
	if ( matched ) {
	  m_SystInfoList.push_back(sys);
	}
      }
    }
  }
  else {
    // fill with default "no systematics"
    ST::SystInfo infodef;
    infodef.affectsKinematics = false;
    infodef.affectsWeights = false;
    infodef.affectsType = ST::Unknown;
      m_SystInfoList.push_back(infodef);
  }
  for ( const auto& sys : m_SystInfoList ) {
    std::cout << "systematics: " << sys.systset.name()  << std::endl;
  }
}


bool BuildSUSYObjects::processEvent(xAOD::TEvent& event)
{

  // for muon trigger SF
  if ( ! m_SUSYObjTool->setRunNumber(267639).isSuccess() ) throw std::runtime_error("Could not set reference run number in SUSYTools !");


  // active storage to put the physics object collections
  xAOD::TStore* store = xAOD::TActiveStore::store();

  // event info
  const xAOD::EventInfo* eventInfo = 0;
  if ( ! event.retrieve( eventInfo, "EventInfo").isSuccess() ) throw std::runtime_error("BuildSUSYObjects: Could not retrieve EventInfo");


  // make to use no systematics
  m_SUSYObjTool->resetSystematics();

  //----------------------------------------   Taus
  std::pair< xAOD::TauJetContainer*, xAOD::ShallowAuxContainer* > susytaus = std::make_pair< xAOD::TauJetContainer*, xAOD::ShallowAuxContainer* >(NULL,NULL);
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
  if ( !m_IsData ) m_tauTruthMatchTool->initializeEvent();
  for( ; tau_itr != tau_end; ++tau_itr ) {
    if ( !m_IsData ) m_tauTruthMatchTool->getTruth( **tau_itr);
    if ( ! m_SUSYObjTool->FillTau( **tau_itr).isSuccess() ) throw std::runtime_error("Error in FillTau");
    if(  !m_IsData && (*tau_itr)->auxdecor<char>("baseline")==1 ){
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



  // loop over systematics variations
  std::vector<float>* event_weights = new std::vector<float>;
  std::vector<std::string>* event_weights_names = new std::vector<std::string>;
  std::vector<CP::SystematicSet>* sys_variations_kinematics = new std::vector<CP::SystematicSet>;
  for  ( const auto & sys : m_SystInfoList ) {
    const CP::SystematicSet& systSet = sys.systset;
    std::string tag = "";
    if ( ! systSet.empty() ) tag = "_"+systSet.name()+"_";

    if ( m_SUSYObjTool->applySystematicVariation(systSet) != CP::SystematicCode::Ok) {
      throw std::runtime_error("Could not apply systematics "+systSet.name());
    }

    bool storeVariation = systSet.name()=="" || sys.affectsKinematics;

    //----------------------------------------   Jets
    xAOD::JetContainer* jets = 0;
    xAOD::ShallowAuxContainer* jets_aux = 0;
    if (! m_SUSYObjTool->GetJets(jets,jets_aux,false).isSuccess() ) {
      throw std::runtime_error("Could not retrieve Jets");
    }
    if ( storeVariation ) {
      if ( ! store->record(jets,"SUSYJets"+m_suffix+tag).isSuccess() ) {
	throw std::runtime_error("Could not store SUSYJets"+m_suffix+tag);
      }
      if ( ! store->record(jets_aux,"SUSYJets"+m_suffix+tag+"Aux.").isSuccess() ) {
	throw std::runtime_error("Could not store SUSYJets"+m_suffix+tag+"Aux.");
      }
    }

    //----------------------------------------   Muons
    xAOD::MuonContainer* muons = 0;
    xAOD::ShallowAuxContainer* muons_aux = 0;
    if (! m_SUSYObjTool->GetMuons(muons,muons_aux,false).isSuccess() ) {
      throw std::runtime_error("Could not retrieve Muons");
    }
    if ( storeVariation ) {
      if ( ! store->record(muons,"SUSYMuons"+m_suffix+tag).isSuccess() ) {
	throw std::runtime_error("Could not store SUSYMuons"+m_suffix+tag);
      }
      if ( ! store->record(muons_aux,"SUSYMuons"+m_suffix+tag+"Aux.").isSuccess() ) {
	throw std::runtime_error("Could not store SUSYMuons"+m_suffix+tag+"Aux.");
      }
    }
    //out() <<  "Muons"+m_suffix+tag+" muons " << std::endl;
    for ( const auto& mu : *muons ) {
      // declare calo and forward muons non baseline so they don't get used in MET
      if ( mu->muonType() != xAOD::Muon::Combined &&
	   mu->muonType() != xAOD::Muon::MuonStandAlone && 
	   mu->muonType() != xAOD::Muon::SegmentTagged ) {
	mu->auxdecor<char>("baseline") = 0;
      }
      m_SUSYObjTool->IsSignalMuon(*mu,10000.);
      m_SUSYObjTool->IsCosmicMuon(*mu);
      m_SUSYObjTool->IsBadMuon(*mu);
      /*
      out() << " Muon " << mu->pt() << " " << mu->eta()
	    << " " << mu->phi() 
	    << " bad " <<  (int) mu->auxdata<char>("bad") 
	    << " baseline " <<  (int) mu->auxdata<char>("baseline") 
	    << " signal " <<  (int) mu->auxdata<char>("signal") 
	    <<std::endl;
      */
    }


     //----------------------------------------   Electrons
    xAOD::ElectronContainer* electrons = 0;
    xAOD::ShallowAuxContainer* electrons_aux = 0;
    if (! m_SUSYObjTool->GetElectrons(electrons,electrons_aux,false).isSuccess() ) {
      throw std::runtime_error("Could not retrieve Electrons");
    }
    if ( storeVariation ) {
      if ( ! store->record(electrons,"SUSYElectrons"+m_suffix+tag).isSuccess() ) {
	throw std::runtime_error("Could not store SUSYElectrons"+m_suffix+tag);
      }
      if ( ! store->record(electrons_aux,"SUSYElectrons"+m_suffix+tag+"Aux.").isSuccess() ) {
	throw std::runtime_error("Could not store SUSYElectrons"+m_suffix+tag+"Aux.");
      }
    }
    //out() <<  "Electrons"+m_suffix+tag+" electrons " << std::endl;
    for ( const auto& el : *electrons ) {
      m_SUSYObjTool->IsSignalElectron(*el,10000.);
      /*
      out() << " Electron " << el->pt() << " " << el->eta()
	    << " " << el->phi() 
	    << " baseline " <<  (int) el->auxdata<char>("baseline") 
	    << " signal " <<  (int) el->auxdata<char>("signal") 
	    <<std::endl;
      */
    }


     //----------------------------------------   Photons
    xAOD::PhotonContainer* photons = 0;
    xAOD::ShallowAuxContainer* photons_aux = 0;
    if (! m_SUSYObjTool->GetPhotons(photons,photons_aux,false).isSuccess() ) {
      throw std::runtime_error("Could not retrieve Photons");
    }
    if ( storeVariation ) {
      if ( ! store->record(photons,"SUSYPhotons"+m_suffix+tag).isSuccess() ) {
	throw std::runtime_error("Could not store SUSYPhotons"+m_suffix+tag);
      }
      if ( ! store->record(photons_aux,"SUSYPhotons"+m_suffix+tag+"Aux.").isSuccess() ) {
	throw std::runtime_error("Could not store SUSYPhotons"+m_suffix+tag+"Aux.");
      }
    }
    //out() <<  "Photons"+m_suffix+tag+" photons " << std::endl;
    float phSF = 1;
    for ( const auto& ph : *photons ) {
      m_SUSYObjTool->IsSignalPhoton(*ph);
      if ( !m_IsData && ph->auxdata<char>("signal") != 0 ) {
	float sf = m_SUSYObjTool->GetSignalPhotonSF(*ph);
	ph->auxdecor<float>("sf") = sf;
	phSF *= sf;
      }
      else
      {
	ph->auxdecor<float>("sf") = 1.;
      }
      /*
      out() << " Photon " << ph->pt() << " " << ph->eta()
	    << " " << ph->phi() 
	    << " baseline " <<  (int) ph->auxdata<char>("baseline") 
	    << " signal " <<  (int) ph->auxdata<char>("signal") 
	    <<std::endl;
      */
    }
    eventInfo->auxdecor<float>("phSF") = phSF ; 

    // Overlap removal
    if ( m_PhotonInOR ) {
      if ( ! m_SUSYObjTool->OverlapRemoval(electrons, muons, jets, photons).isSuccess() ) throw std::runtime_error("Error in OverlapRemoval");
    }
    else {
      if ( ! m_SUSYObjTool->OverlapRemoval(electrons, muons, jets).isSuccess() ) throw std::runtime_error("Error in OverlapRemoval");
    }

    // signal and btag jet now depend on OR, so loop again on jets
    //out() <<  "SUSYJets"+m_suffix+tag+" jets " << std::endl;
    for ( const auto& jet : *jets ) {
      if ( !m_UseSmearedJets ) { 
	m_SUSYObjTool->IsBadJet(*jet); // no JVT cut
	// currently does not add any decoration
	bool issignal = m_SUSYObjTool->IsSignalJet(*jet);
	if ( issignal ) {
	  jet->auxdecor<char>("signal") = 1;
	}
	else {
	  jet->auxdecor<char>("signal") = 0;
	}
      }
      if ( m_period == p8tev ) {
	//FIXME m_SUSYObjTool->IsBJet(**jet_itr,false,1.85);
      }
      else if ( m_period == p13tev ) {
	m_SUSYObjTool->IsBJet(*jet);
      }
      /*
      out() << " Jet " << jet->pt() << " " << jet->eta()
	    << " " << jet->phi() 
	    << " bad " <<  (int) jet->auxdata<char>("bad") 
	    << " baseline " <<  (int) jet->auxdata<char>("baseline") 
	    << " signal " <<  (int) jet->auxdata<char>("signal") 
	    <<std::endl;
      */
    }

    // GetTotalMuonSF also test for OR
    float muSF = 1.f;
    if ( !m_IsData ) muSF = (float) m_SUSYObjTool->GetTotalMuonSF(*muons,true,true,"HLT_mu20_iloose_L1MU15_OR_HLT_mu50");
    eventInfo->auxdecor<float>("muSF") = muSF ; 

    // idem for GetTotalElectronSF
    float elSF = 1.f;
    if ( !m_IsData ) elSF = (float) m_SUSYObjTool->GetTotalElectronSF(*electrons);
    eventInfo->auxdecor<float>("elSF") = elSF ; 



    xAOD::MissingETContainer* rebuiltmetc = new xAOD::MissingETContainer();
    xAOD::MissingETAuxContainer* rebuiltmetcAux = new xAOD::MissingETAuxContainer();
    rebuiltmetc->setStore(rebuiltmetcAux);
    if ( storeVariation ) {
      if ( ! store->record(rebuiltmetc,"MET_ZL"+m_suffix+tag).isSuccess() ) {
	throw std::runtime_error("Unable to store MissingETContainer with tag MET_ZL"+m_suffix+tag);
      }
      if ( ! store->record(rebuiltmetcAux,"MET_ZL"+m_suffix+tag+"Aux.").isSuccess() ) {
	throw std::runtime_error("Unable to store MissingETAuxContainer with tag MET_ZL"+m_suffix+tag+"Aux");
      }
    }
    
    if ( ! m_SUSYObjTool->GetMET(*rebuiltmetc,
				 jets,
				 electrons,
				 muons,
				 photons,
				 0,
				 true,
				 true,
				 0).isSuccess() 
	 ) throw std::runtime_error("Error in GetMET");
    TVector2* MissingET = new TVector2(0.,0.);
    xAOD::MissingETContainer::const_iterator met_it = rebuiltmetc->find("Final");
    if ( met_it == rebuiltmetc->end() ) throw std::runtime_error("Could not find Final MET after running  GetMET");
    MissingET->Set((*met_it)->mpx(), (*met_it)->mpy());
    if ( ! store->record(MissingET,"SUSYMET"+m_suffix+tag).isSuccess() ) {
      throw std::runtime_error("Could not store SUSYMET with tag SUSYMET"+m_suffix+tag);
    }

    if ( sys.affectsKinematics || systSet.name()=="" ) {
      sys_variations_kinematics->push_back(systSet);
    }
    if ( sys.affectsWeights || systSet.name()=="" ) {
      event_weights->push_back(elSF*phSF*muSF);
      event_weights_names->push_back(systSet.name());
    }

    // delete shallow copies if not stored
    if ( !storeVariation ) {
      delete jets;
      delete jets_aux;
      delete muons;
      delete muons_aux;
      delete electrons;
      delete electrons_aux;
      delete photons;
      delete photons_aux;
      delete rebuiltmetc;
      delete rebuiltmetcAux;
    }

  }

  // store variations that affect kinematics -> need to rerun SR+CR
  if ( ! store->record(sys_variations_kinematics,"sys_variations_kinematics"+m_suffix).isSuccess() ) {
    throw std::runtime_error("Could not store sys_variations_kinematics"+m_suffix);
  }
  // store event weights variations
  if ( ! store->record(event_weights_names,"event_weights_names"+m_suffix).isSuccess() ) {
    throw std::runtime_error("Could not store event_weights_names"+m_suffix);
  }
  if ( ! store->record(event_weights,"event_weights"+m_suffix).isSuccess() ) {
    throw std::runtime_error("Could not store event_weights"+m_suffix);
  }




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
    "HLT_e24_lhmedium_iloose_L1EM20VH",
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

