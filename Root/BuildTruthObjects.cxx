#include "ZeroLeptonRun2/BuildTruthObjects.h"

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
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthParticle.h"
#include "xAODTau/TauJetContainer.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODMissingET/MissingETAuxContainer.h"
#include "PATInterfaces/SystematicSet.h"
#include "MCTruthClassifier/MCTruthClassifierDefs.h"
#include "FourMomUtils/xAODP4Helpers.h"

#include "cafe/Config.h"

#include "TVector2.h"

#include <stdexcept>

static const SG::AuxElement::Decorator<char> dec_passOR("passOR");
static const SG::AuxElement::ConstAccessor<unsigned int> acc_truthType("classifierParticleType");

BuildTruthObjects::BuildTruthObjects(const char *name)
  : cafe::Processor(name),
    m_IsData(false),
    m_PhotonInOR(false),
    m_jetkey(),
    m_suffix()
{
  cafe::Config config(name);
  m_IsData = config.get("IsData",false);
  m_jetkey = config.get("JetCollectionKey","xxxx");
  m_PhotonInOR = config.get("PhotonInOR",false);
  m_suffix = config.get("suffix","");


}

BuildTruthObjects::~BuildTruthObjects()
{

}

bool BuildTruthObjects::processEvent(xAOD::TEvent& event)
{

  xAOD::TStore* store = xAOD::TActiveStore::store();

  // Jets
  const xAOD::JetContainer* inputjets = 0;
  xAOD::JetContainer* outputjets = 0;
  if ( !event.retrieve(inputjets, m_jetkey).isSuccess() ) {
    throw std::runtime_error("Could not retrieve JetContainer with key "+m_jetkey);
  }
  std::pair< xAOD::JetContainer*, xAOD::ShallowAuxContainer* > truthjets = xAOD::shallowCopyContainer( *inputjets );
  // if ( ! store->record<xAOD::JetContainer>(truthjets.first,"TruthJets"+m_suffix).isSuccess() ) { // BEFORE
  if ( ! store->record(truthjets.first,"myTruthJets"+m_suffix).isSuccess() ) {
    throw std::runtime_error("Could not store myTruthJets"+m_suffix);
  }
  //  if ( ! store->record<xAOD::ShallowAuxContainer>(truthjets.second,"TruthJets"+m_suffix+"Aux.").isSuccess() ) { //BEFORE
  if ( ! store->record(truthjets.second,"myTruthJets"+m_suffix+"Aux.").isSuccess() ) {
    throw std::runtime_error("Could not store myTruthJets"+m_suffix+"Aux.");
  }
  outputjets = truthjets.first;

  if ( !xAOD::setOriginalObjectLink(*inputjets, *outputjets) ) throw std::runtime_error("Could not set original links in jet container copy");

  // Muons
  const xAOD::TruthParticleContainer* truthmuons = 0;
  if ( !event.retrieve( truthmuons, "TruthMuons").isSuccess() ) {
    throw std::runtime_error("Could not retrieve truth particles with key TruthMuons");
  }

  std::pair< xAOD::TruthParticleContainer*, xAOD::ShallowAuxContainer* > truthmuons2 = xAOD::shallowCopyContainer( *truthmuons );
  //  if ( ! store->record<xAOD::TruthParticleContainer>(truthmuons2.first,"TruthMuons").isSuccess() ) { //BEFORE
  if ( ! store->record(truthmuons2.first,"myTruthMuons"+m_suffix).isSuccess() ) {
    throw std::runtime_error("Could not store myTruthMuons"+m_suffix);
  }
  //if ( ! store->record<xAOD::ShallowAuxContainer>(truthmuons2.second,"TruthMuonsAux.").isSuccess()) { //BEFORE
  if ( ! store->record(truthmuons2.second,"myTruthMuons"+m_suffix+"Aux.").isSuccess()) {
    throw std::runtime_error("Could not store myTruthMuons"+m_suffix+"Aux.");
  }

  // Electrons
  const xAOD::TruthParticleContainer* truthelectrons = 0 ;
  if ( !event.retrieve(truthelectrons, "TruthElectrons").isSuccess() ){
    throw std::runtime_error("Could not retrieve truth particles with key TruthElectrons");
  }

  std::pair< xAOD::TruthParticleContainer*, xAOD::ShallowAuxContainer* > truthelectrons2 = xAOD::shallowCopyContainer( *truthelectrons );
  // if ( ! store->record<xAOD::TruthParticleContainer>(truthelectrons2.first,"TruthElectrons"+m_suffix).isSuccess() ) { //BEFORE
  if ( ! store->record(truthelectrons2.first,"myTruthElectrons"+m_suffix).isSuccess() ) {
    throw std::runtime_error("Could not store myTruthElectrons"+m_suffix);
  }
  //if ( ! store->record<xAOD::ShallowAuxContainer>(truthelectrons2.second,"TruthElectrons"+m_suffix+"Aux.").isSuccess() ) { //BEFORE
  if ( ! store->record(truthelectrons2.second,"myTruthElectrons"+m_suffix+"Aux.").isSuccess() ) {
    throw std::runtime_error("Could not store myTruthElectrons"+m_suffix+"Aux.");
  }


  // Photons
  const xAOD::TruthParticleContainer* truthphotons = 0 ;
  if ( !event.retrieve(truthphotons,"TruthPhotons").isSuccess() ){
    throw std::runtime_error("Could not retrieve truth particles with key TruthPhotons");
  }

  std::pair< xAOD::TruthParticleContainer*, xAOD::ShallowAuxContainer* > truthphotons2 = xAOD::shallowCopyContainer( *truthphotons );

  if ( ! store->record(truthphotons2.first,"myTruthPhotons"+m_suffix).isSuccess() ) {
    throw std::runtime_error("Could not store myTruthPhotons"+m_suffix);
  }
  if ( ! store->record(truthphotons2.second,"myTruthPhotons"+m_suffix+"Aux.").isSuccess()) {
    throw std::runtime_error("Could not store myTruthPhotons"+m_suffix+"Aux.");
  }



  // Taus
  const xAOD::TruthParticleContainer* truthtaus = 0 ;

  if (!event.retrieve(truthtaus,"TruthTaus").isSuccess()){
    throw std::runtime_error("Could not retrieve truth particles with key TruthTaus");
  }

  std::pair< xAOD::TruthParticleContainer*, xAOD::ShallowAuxContainer* > truthtaus2 = xAOD::shallowCopyContainer( *truthtaus );

  if ( ! store->record(truthtaus2.first,"myTruthTaus"+m_suffix).isSuccess() ) {
    throw std::runtime_error("Could not store myTruthTaus"+m_suffix);
  }
  if ( ! store->record(truthtaus2.second,"myTruthTaus"+m_suffix+"Aux.").isSuccess()) {
    throw std::runtime_error("Could not store myTruthTaus"+m_suffix+"Aux.");
  }



  // Overlap removal

  if ( m_PhotonInOR ) {
    if (! OverlapRemoval(truthelectrons2.first, truthmuons2.first, outputjets, truthphotons2.first, false,
			 0.2 /*dRejet*/, 0.4/*dRjetmu*/, 0.4/*dRjete*/, 0.01/*dRemu*/
			 ,0.05/*dRee*/, 0.2/*dRphjet*/, 0.4/*dReph*/, 0.000/*dRmuph*/) ) throw std::runtime_error("Error in OverlapRemoval");
  }
  else{
  if (! OverlapRemoval(truthelectrons2.first, truthmuons2.first, outputjets, false, 0.2, 0.4, 0.4, 0.01, 0.05) ) throw std::runtime_error("Error in OverlapRemoval");
  }


  // Missing ET
  const xAOD::MissingETContainer* truthCont = 0 ;
  if ( !event.retrieve(truthCont, "MET_Truth").isSuccess() ){
    throw std::runtime_error("Could not retrieve truth met with key MET_Truth");
  }

  TVector2* MissingET = new TVector2(0.,0.);
  const xAOD::MissingET* met = truthCont->front();

  MissingET->Set(met->mpx(), met->mpy());

  if ( ! store->record(MissingET,"TruthMET"+m_suffix).isSuccess() ) {
    throw std::runtime_error("Could not store TruthMET"+m_suffix);
  }

return true;
}


bool BuildTruthObjects::OverlapRemoval(const xAOD::TruthParticleContainer *electrons, const xAOD::TruthParticleContainer *muons, const xAOD::JetContainer *jets, bool doHarmonization, double dRejet, double dRjetmu, double dRjete, double dRemu, double dRee)
{
  int Njetin=0;
  int Nelin=0;
  int Nmuin=0;

  for(const auto& jet : *jets) {
    bool jet_sel = jet->pt() > 20000 && std::abs(jet->eta()) < 2.8;
    dec_passOR(*jet) = jet_sel;
    if(jet_sel) Njetin++;
  }
  for(const auto& mu : *muons) {
    bool mu_sel = mu->pt() > 10000 && std::abs(mu->eta()) < 2.4 && acc_truthType(*mu)==MCTruthPartClassifier::IsoMuon;
    dec_passOR( *mu ) = mu_sel;
    if(mu_sel) Nmuin++;
  }
  for(const auto& el : *electrons) {
    bool el_sel = el->pt() > 10000 && std::abs(el->eta()) < 2.47 && acc_truthType(*el)==MCTruthPartClassifier::IsoElectron;
    dec_passOR( *el ) = el_sel;
    if(el_sel) Nelin++;
    else continue;

    //#######################

    for(const auto& jet : *jets) {
      if( !dec_passOR(*jet) ) continue;
      if (xAOD::P4Helpers::isInDeltaR(*el, *jet, dRejet)) {
	//std::cout << " Rejecting jet at (eta,phi)=(" << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10) <<  jet->eta() <<"," << jet->phi() <<") "
	//	  << " due to electron at (eta,phi)=(" << el->eta() <<"," << el->phi() <<")" << std::endl ;
	dec_passOR( *jet ) = false;
      }
    }
  } // END loop over electrons
  // Remove electrons and muons overlapping with jets
  for(const auto& el : *electrons) {
    if( !dec_passOR(*el) ) continue;

    for(const auto& jet : *jets) {
      if ( !dec_passOR( *jet ) ) continue;

      if (xAOD::P4Helpers::isInDeltaR(*el, *jet, dRjete)) {
	dec_passOR(*el) = false;
      }
    }
  }

  for(const auto& mu : *muons) {
    if( !dec_passOR(*mu) ) continue;

    for(const auto& jet : *jets) {
      if ( !dec_passOR( *jet ) ) continue;

      if (xAOD::P4Helpers::isInDeltaR(*mu, *jet, dRjetmu)) {
	dec_passOR( *mu ) = false;
      }
    }
  }

  // Remove electrons and muons overlapping with each other
  for(const auto& el : *electrons) {
    if( !dec_passOR(*el) ) continue;

    for(const auto& mu : *muons) {
      if ( !dec_passOR( *mu ) ) continue;

      if (xAOD::P4Helpers::isInDeltaR(*el, *mu, dRemu)) {
	dec_passOR( *el ) = false;
	dec_passOR( *mu ) = false;
      }
    }
  }
  // Remove electrons overlapping with each other
  for(const auto& el : *electrons) {
    if( !dec_passOR(*el) ) continue ;

    for(const auto& el2 : *electrons) {
      if(el == el2) continue;
      if ( !dec_passOR( *el2 ) ) continue;

      if (xAOD::P4Helpers::isInDeltaR(*el, *el2, dRee)) {
	if(el->pt() < el2->pt()){
	  //std::cout << " Rejecting electron at (eta,phi)=(" << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10)  << el->eta() <<"," << el->phi() <<") "
	  //	    << " and muon at (eta,phi)=(" << (*el2_itr)->eta() <<"," << (*el2_itr)->phi() <<")"<< std::endl ;
	  dec_passOR(*el) = false;
	}else{
	  //std::cout << " Rejecting electron at (eta,phi)=(" << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10) << (*el2_itr)->eta() <<"," << (*el2_itr)->phi() <<") "
	  //	    << " and muon at (eta,phi)=(" << el->eta() <<"," << el->phi() <<")"<< std::endl ;
	  dec_passOR( *el2 ) = false;
	}
      }
    }
  }
  // Count number of objects after overlap removal
  int Nel=0;
  for(const auto& el : *electrons) {
    if(dec_passOR( *el )) Nel++;
  }

  int Nmu=0;
  for(const auto& mu : *muons) {
    if(dec_passOR( *mu )) Nmu++;
  }

  int Njet=0;
  for(const auto& jet : *jets) {
    if(dec_passOR( *jet )) Njet++;
  }
  //std::cout << " Before overlap removal: Nel=" << Nelin <<", Nmu="<< Nmuin <<", Njet=" << Njetin<< std::endl ;
  //std::cout << " After  overlap removal: Nel=" << Nel <<", Nmu="<< Nmu <<", Njet=" << Njet<< std::endl ;
  return true;

}

bool BuildTruthObjects::OverlapRemoval(const xAOD::TruthParticleContainer *electrons, const xAOD::TruthParticleContainer *muons, const xAOD::JetContainer *jets, const xAOD::TruthParticleContainer *photons, const bool doHarmonization, const double dRejet, const double dRjetmu, const double dRjete, double dRemu, double dRee, double dRphjet, double dReph, double dRmuph){

  int Njetin=0;
  int Nelin=0;
  int Nmuin=0;
  int Nphin=0;

  for(const auto& jet : *jets) {
    bool jet_sel = jet->pt() > 20000 && std::abs(jet->eta()) < 2.8;
    dec_passOR(*jet) = jet_sel;
    if(jet_sel) Njetin++;
  }
  for(const auto& mu : *muons) {
    bool mu_sel = mu->pt() > 10000 && std::abs(mu->eta()) < 2.4 && acc_truthType(*mu)==MCTruthPartClassifier::IsoMuon;
    dec_passOR( *mu ) = mu_sel;
    if(mu_sel) Nmuin++;
  }
  for(const auto& el : *electrons) {
    bool el_sel = el->pt() > 10000 && std::abs(el->eta()) < 2.47 && acc_truthType(*el)==MCTruthPartClassifier::IsoElectron;
    dec_passOR( *el ) = el_sel;
    if(el_sel) Nelin++;
    else continue;

    for(const auto& jet : *jets) {
      if( !dec_passOR(*jet) ) continue;

      if (xAOD::P4Helpers::isInDeltaR(*el, *jet, dRejet)) {
	//std::cout << " Rejecting jet at (eta,phi)=(" << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10)  << jet->eta() <<"," << jet->phi() <<") "
	//	  << " due to electron at (eta,phi)=(" << el->eta() <<"," << el->phi() <<")"<< std::endl ;
	dec_passOR(*jet) = false;
      }
    }
  } // END loop over electrons

  //std::cout << " Before overlap removal: Nel=" << Nelin <<", Nmu="<< Nmuin <<", Njet=" << Njetin<< std::endl ;

  // remove jets overlapping with (baseline/signal) photons
  for(const auto& ph : *photons) {
    bool ph_sel = ph->pt() > 25000 && std::abs(ph->eta()) < 2.37 && acc_truthType(*ph)==MCTruthPartClassifier::IsoPhoton;
    dec_passOR( *ph ) = ph_sel;
    if(ph_sel) Nphin++;
    else  continue;

    for(const auto& jet : *jets) {
      if( !dec_passOR(*jet) ) continue ;
      if (xAOD::P4Helpers::isInDeltaR(*ph, *jet, dRphjet)) {
	// std::cout << " Rejecting jet at (eta,phi)=(" << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10) << jet->eta() <<"," << jet->phi() <<") "
	// 	  << " due to photon at (eta,phi)=(" << ph->eta() <<"," << ph->phi() <<")"<< std::endl ;
	dec_passOR(*jet) = false;
      }
    }
  }// END loop over photons

  // Remove electrons and muons overlapping with jets and photons
  for(const auto& el : *electrons) {
    if( !dec_passOR(*el) ) continue ;

    for(const auto& jet : *jets) {
      if( !dec_passOR(*jet) ) continue;

      if (xAOD::P4Helpers::isInDeltaR(*el, *jet, dRjete)) {
	//std::cout <<  " Rejecting electron at (eta,phi)=(" << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10)  << el->eta() <<"," << el->phi() <<") "
	//	  << " due to jet at (eta,phi)=(" << jet->eta() <<"," << jet->phi() <<")"<< std::endl ;
	dec_passOR(*el) = false;
      }
    }
  }


  for(const auto& mu : *muons) {
    if( !dec_passOR(*mu) ) continue;

    for(const auto& jet : *jets) {
      if( !dec_passOR(*jet) ) continue;

      if (xAOD::P4Helpers::isInDeltaR(*mu, *jet, dRjetmu)) {
	dec_passOR(*mu) = false;
      }
    }
  }



  // Remove electrons and muons overlapping with each other
  for(const auto& el : *electrons) {

    if( !dec_passOR(*el) ) continue ;//if( !dec_passOR(*el) ) continue;

    for(const auto& mu : *muons) {
      if( !dec_passOR(*mu) ) continue ; //if ( !dec_passOR( *mu ) ) continue;

      if (xAOD::P4Helpers::isInDeltaR(*el, *mu, dRemu)) {
	//std::cout << " Rejecting both electron at (eta,phi)=(" << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10)  << el->eta() <<"," << el->phi() <<") "
	//	  << " and muon at (eta,phi)=(" << mu->eta() <<"," << mu->phi() <<")"<< std::endl ;
	dec_passOR(*el) = false;
	dec_passOR(*mu) = false;
      }
    }
  }


  // Remove electrons overlapping with each other
  for(const auto& el : *electrons) {

    if( !dec_passOR(*el) ) continue ;  //if( !dec_passOR(*el) ) continue;

    for(const auto& el2 : *electrons ) {

      if(el == el2) continue;
      if( !dec_passOR(*el2) ) continue;

      if (xAOD::P4Helpers::isInDeltaR(*el, *el2, dRee)) {
	if(el->pt() < el2->pt()){
	  dec_passOR( *el ) = false;
	  //std::cout << " Rejecting electron at (eta,phi)=(" << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10) << el->eta() <<"," << el->phi() <<") "
	  //	    << " and muon at (eta,phi)=(" << (*el2_itr)->eta() <<"," << (*el2_itr)->phi() <<")"<< std::endl ;
	}else{
	  dec_passOR( *el2 ) = false;
	  //std::cout << " Rejecting electron at (eta,phi)=(" << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10)  << (*el2_itr)->eta() <<"," << (*el2_itr)->phi() <<") "
	  //	    << " and muon at (eta,phi)=(" << el->eta() <<"," << el->phi() <<")"<< std::endl ;
	}
      }
    }
  }


  // Remove photons if overlapping with electrons

  for(const auto& ph : *photons) {

    if( !dec_passOR(*ph) ) continue ; //if( !dec_passOR( *ph ) )
    continue;

    for(const auto& el : *electrons) {
      if( !dec_passOR(*el) ) continue ; //if( !dec_passOR(*el) ) continue;

      if (xAOD::P4Helpers::isInDeltaR(*el, *ph, dReph)) {
	//std::cout << " Rejecting photon at (eta,phi)=(" << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10)  << ph->eta() <<"," << ph->phi() <<") "
	//	  << " due to electron at (eta,phi)=(" << el->eta() <<"," << el->phi() <<")"<< std::endl ;
	dec_passOR( *ph ) = false;
      }
    }
  }


  // Remove photons if overlapping with muons
  for(const auto& ph : *photons) {

    if( !dec_passOR(*ph) ) continue ; //if( !dec_passOR( *ph ) )
    continue;

    for(const auto& mu : *muons) {
      if( !dec_passOR(*mu) ) continue ; //if( !dec_passOR(*mu) ) continue;
      if (xAOD::P4Helpers::isInDeltaR(*mu, *ph, dRmuph)) {
	//std::cout << " Rejecting photon at (eta,phi)=(" << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10) << ph->eta() <<"," << ph->phi() <<") "
	//	  << " due to muon at (eta,phi)=(" << mu->eta() <<"," << mu->phi() <<")"<< std::endl ;
	dec_passOR( *ph ) = false;
      }
    }
  }


  // Count number of objects after overlap removal
  int Nel=0;
  for(const auto& el : *electrons) {
    if(dec_passOR(*el) == 1) Nel++ ;//if(dec_passOR( *el )) Nel++;
  }

  int Nmu=0;
  for(const auto& mu : *muons) {
    if(dec_passOR(*mu) == 1) Nmu++ ;//if(dec_passOR( *mu )) Nmu++;
  }

  int Njet=0;
  for(const auto& jet : *jets) {
    if(dec_passOR(*jet) == 1) Njet++ ;//if(dec_passOR( *jet )) Njet++;
  }

  int Nph=0;
  for(const auto& ph : *photons) {
    if(dec_passOR(*ph) == 1) Nph++ ;//if(dec_passOR( *ph )) Nph++;
  }

  //std::cout << " After overlap removal: Nel=" << Nel <<", Nmu="<< Nmu <<", Njet=" << Njet << ", Nph=" << Nph << std::endl;
  return true;

}




ClassImp(BuildTruthObjects);

