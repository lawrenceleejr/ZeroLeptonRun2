

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


#include "cafe/Config.h"

#include "TVector2.h"

#include <stdexcept>

BuildTruthObjects::BuildTruthObjects(const char *name)
  : cafe::Processor(name),
    m_IsData(false),
    m_jetkey(),
    m_PhotonInOR(false),
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
  if ( ! store->record(truthjets.first,"TruthJets"+m_suffix).isSuccess() ) {
    throw std::runtime_error("Could not store TruthJets"+m_suffix);
  }
  //  if ( ! store->record<xAOD::ShallowAuxContainer>(truthjets.second,"TruthJets"+m_suffix+"Aux.").isSuccess() ) { //BEFORE
  if ( ! store->record(truthjets.second,"TruthJets"+m_suffix+"Aux.").isSuccess() ) {
    throw std::runtime_error("Could not store TruthJets"+m_suffix+"Aux.");
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
  if ( ! store->record(truthmuons2.first,"TruthMuons"+m_suffix).isSuccess() ) { 
    throw std::runtime_error("Could not store truthMuons"+m_suffix);
  }
  //if ( ! store->record<xAOD::ShallowAuxContainer>(truthmuons2.second,"TruthMuonsAux.").isSuccess()) { //BEFORE
  if ( ! store->record(truthmuons2.second,"TruthMuons"+m_suffix+"Aux.").isSuccess()) {
    throw std::runtime_error("Could not store truthMuons"+m_suffix+"Aux.");
  }
  
  // Electrons
  const xAOD::TruthParticleContainer* truthelectrons = 0 ; 
  if ( !event.retrieve(truthelectrons, "TruthElectrons").isSuccess() ){
    throw std::runtime_error("Could not retrieve truth particles with key TruthElectrons");
  }

  std::pair< xAOD::TruthParticleContainer*, xAOD::ShallowAuxContainer* > truthelectrons2 = xAOD::shallowCopyContainer( *truthelectrons );
  // if ( ! store->record<xAOD::TruthParticleContainer>(truthelectrons2.first,"TruthElectrons"+m_suffix).isSuccess() ) { //BEFORE
  if ( ! store->record(truthelectrons2.first,"TruthElectrons"+m_suffix).isSuccess() ) {
    throw std::runtime_error("Could not store TruthElectrons"+m_suffix);
  }
  //if ( ! store->record<xAOD::ShallowAuxContainer>(truthelectrons2.second,"TruthElectrons"+m_suffix+"Aux.").isSuccess() ) { //BEFORE
  if ( ! store->record(truthelectrons2.second,"TruthElectrons"+m_suffix+"Aux.").isSuccess() ) {
    throw std::runtime_error("Could not store TruthElectrons"+m_suffix+"Aux.");
  }

  
  // Photons
  const xAOD::TruthParticleContainer* truthphotons = 0 ;
  if ( !event.retrieve(truthphotons,"TruthPhotons").isSuccess() ){
    throw std::runtime_error("Could not retrieve truth particles with key TruthPhotons");
  }
  
  std::pair< xAOD::TruthParticleContainer*, xAOD::ShallowAuxContainer* > truthphotons2 = xAOD::shallowCopyContainer( *truthphotons );
  
  if ( ! store->record(truthphotons2.first,"TruthPhotons"+m_suffix).isSuccess() ) {
    throw std::runtime_error("Could not store TruthPhotons"+m_suffix);
  }
  if ( ! store->record(truthphotons2.second,"TruthPhotons"+m_suffix+"Aux.").isSuccess()) {
    throw std::runtime_error("Could not store TruthPhotons"+m_suffix+"Aux.");
  }


  
  // Taus
  const xAOD::TruthParticleContainer* truthtaus = 0 ;

  if (!event.retrieve(truthtaus,"TruthTaus").isSuccess()){
      throw std::runtime_error("Could not retrieve truth particles with key TruthTaus");
    }

  std::pair< xAOD::TruthParticleContainer*, xAOD::ShallowAuxContainer* > truthtaus2 = xAOD::shallowCopyContainer( *truthtaus );
  
  if ( ! store->record(truthtaus2.first,"TruthTaus"+m_suffix).isSuccess() ) {
    throw std::runtime_error("Could not store TruthTaus"+m_suffix);
  }
  if ( ! store->record(truthtaus2.second,"TruthTaus"+m_suffix+"Aux.").isSuccess()) {
    throw std::runtime_error("Could not store TruthTaus"+m_suffix+"Aux.");
  }
  
  

  // Overlap removal 

  if ( m_PhotonInOR ) {
    if (! OverlapRemoval(truthelectrons2.first, truthmuons2.first, outputjets, truthphotons2.first, false, 0.2, 0.4, 0.4, 0.01, 0.05, 0.2, 0.4, 0.000) ) throw std::runtime_error("Error in OverlapRemoval");
  }
  else{
  if (! OverlapRemoval(truthelectrons2.first, truthmuons2.first, outputjets, false, 0.2, 0.4, 0.4, 0.01, 0.05) ) throw std::runtime_error("Error in OverlapRemoval");
  }


  // Missing ET
  const xAOD::MissingETContainer* truthCont = 0 ; //new xAOD::MissingETContainer();
  if ( !event.retrieve(truthCont, "MET_Truth").isSuccess() ){
    throw std::runtime_error("Could not retrieve truth met with key MET_Truth");
  }
  
  //TVector2* MissingE1 = new TVector2(0.,0.);
  //TVector2* MissingE2 = new TVector2(0.,0.);
  //TVector2* MissingE3 = new TVector2(0.,0.);
  //TVector2* MissingE4 = new TVector2(0.,0.);
  TVector2* MissingET = new TVector2(0.,0.);
  //
  //
  xAOD::MissingETContainer::const_iterator met_it1 = truthCont->begin();
  //xAOD::MissingETContainer::const_iterator met_it2 = truthCont->begin()+1;
  //xAOD::MissingETContainer::const_iterator met_it3 = truthCont->begin()+2;
  //xAOD::MissingETContainer::const_iterator met_it4 = truthCont->begin()+3;

  MissingET->Set((*met_it1)->mpx(), (*met_it1)->mpy());

  //if ( met_it == truthCont->end() ) throw std::runtime_error("Could not find MET after running weird function");
  //MissingE1->Set((*met_it1)->mpx(), (*met_it1)->mpy());
  //MissingE2->Set((*met_it2)->mpx(), (*met_it2)->mpy());
  //MissingE3->Set((*met_it3)->mpx(), (*met_it3)->mpy());
  //MissingE4->Set((*met_it4)->mpx(), (*met_it4)->mpy());
  ////
  ////
  //std::cout << "TEST OF MET : " << std::endl;
  //std::cout << "COMP 1 : " << MissingE1->Mod() << "  " << MissingE1->Px() << "  " << MissingE1->Py() << std::endl;
  //std::cout << "COMP 2 : " << MissingE2->Mod() << "  " << MissingE2->Px() << "  " << MissingE2->Py() << std::endl;
  //std::cout << "COMP 3 : " << MissingE3->Mod() << "  " << MissingE3->Px() << "  " << MissingE3->Py() << std::endl;
  //std::cout << "COMP 4 : " << MissingE4->Mod() << "  " << MissingE4->Px() << "  " << MissingE4->Py() << std::endl;


  //if ( ! store->record<TVector2>(MissingET,"TruthMET"+m_suffix).isSuccess() ) {   //BEFORE
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

  xAOD::JetContainer::const_iterator jet_itr = jets->begin();
  xAOD::JetContainer::const_iterator jet_end = jets->end();
  for( ; jet_itr != jet_end; ++jet_itr ) {
    //bool jet_sel = (*jet_itr)->auxdecor<char>("baseline");//bool jet_sel = dec_baseline(**jet_itr);
    bool jet_sel; 
    if( (*jet_itr)->pt() > 20000 && std::abs((*jet_itr)->eta()) < 2.8 )
      jet_sel = 1 ; 
    if(jet_sel){
      (*jet_itr)->auxdecor<char>("passOR") = 1 ; //dec_passOR( **jet_itr ) = 1;
      Njetin++;
    }
    else
      (*jet_itr)->auxdecor<char>("passOR") = 0 ; //dec_passOR( **jet_itr ) = 0;
  }
  xAOD::TruthParticleContainer::const_iterator mu_itr = muons->begin();
  xAOD::TruthParticleContainer::const_iterator mu_end = muons->end();
  for( ; mu_itr != mu_end; ++mu_itr ) {
    bool mu_sel;
    //if(doHarmonization) mu_sel = (*mu_itr)->auxdecor<char>("signal"); //mu_sel = dec_signal(**mu_itr);
    //else mu_sel = (*mu_itr)->auxdecor<char>("baseline"); //mu_sel = dec_baseline(**mu_itr);
    if( (*mu_itr)->pt() > 10000 && std::abs((*mu_itr)->eta()) < 2.4 )
      mu_sel = 1 ; 
    if(mu_sel){
      (*mu_itr)->auxdecor<char>("passOR") = 1 ; //dec_passOR( **mu_itr ) = 1;
      Nmuin++;
    }
    else
      (*mu_itr)->auxdecor<char>("passOR") = 0 ; //dec_passOR( **mu_itr ) = 0;
  }
  // remove jets overlapping with (baseline/signal) electrons
  xAOD::TruthParticleContainer::const_iterator el_itr = electrons->begin();
  xAOD::TruthParticleContainer::const_iterator el_end = electrons->end();
  for( ; el_itr != el_end; ++el_itr ) {
    bool el_sel;
    //if(doHarmonization) el_sel = (*el_itr)->auxdecor<char>("signal"); // el_sel = dec_signal(**el_itr);
    //else el_sel = (*el_itr)->auxdecor<char>("baseline"); // el_sel = dec_baseline(**el_itr);
    if( (*el_itr)->pt() > 20000 && std::abs((*el_itr)->eta()) < 2.47 )
      el_sel = 1 ;
    if( !el_sel ){
      (*el_itr)->auxdecor<char>("passOR") = 0 ; //dec_passOR( **el_itr ) = 0;
      continue;
    }else{
      (*el_itr)->auxdecor<char>("passOR") = 1 ; //dec_passOR( **el_itr ) = 1;
      Nelin++;
    }

    //#######################
    
    xAOD::JetContainer::const_iterator jet_itr = jets->begin();
    xAOD::JetContainer::const_iterator jet_end = jets->end();
    
    for( ; jet_itr != jet_end; ++jet_itr ) {
      if((*jet_itr)->auxdecor<char>("passOR") != 1) continue ; // if( !dec_passOR(**jet_itr) ) continue;
      
      TLorentzVector el4vec = (*el_itr)->p4();
      TLorentzVector jet4vec = (*jet_itr)->p4();
      
      if (el4vec.DeltaR(jet4vec)<dRejet) {
	std::cout << " Rejecting jet at (eta,phi)=(" << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10) <<  (*jet_itr)->eta() <<"," << (*jet_itr)->phi() <<") "
		  << " due to electron at (eta,phi)=(" << (*el_itr)->eta() <<"," << (*el_itr)->phi() <<")" << std::endl ;
	(*jet_itr)->auxdecor<char>("passOR") = 0 ; //dec_passOR( **jet_itr ) = 0;
      }
    }
  } // END loop over electrons
  // Remove electrons and muons overlapping with jets
  el_itr = electrons->begin();
  el_end = electrons->end();
  for( ; el_itr != el_end; ++el_itr ) {
    
    if((*el_itr)->auxdecor<char>("passOR") !=1) continue ; // if( !dec_passOR(**el_itr) ) continue; 
    
    xAOD::JetContainer::const_iterator jet_itr = jets->begin();
    xAOD::JetContainer::const_iterator jet_end = jets->end();
    
    for( ; jet_itr != jet_end; ++jet_itr ) {
      
      if((*jet_itr)->auxdecor<char>("passOR") != 1) continue ;// if ( !dec_passOR( **jet_itr ) ) continue;
      TLorentzVector el4vec = (*el_itr)->p4();
      TLorentzVector jet4vec = (*jet_itr)->p4();
      
      if (el4vec.DeltaR(jet4vec)<dRjete) {
	std::cout << " Rejecting electron at (eta,phi)=(" <<  std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10) << (*el_itr)->eta() <<"," << (*el_itr)->phi() <<") "
		  << " due to jet at (eta,phi)=(" << (*jet_itr)->eta() <<"," << (*jet_itr)->phi() <<")"<< std::endl ;
	(*el_itr)->auxdecor<char>("passOR") = 0 ;//dec_passOR( **el_itr ) = 0;
      }
    }
  }
  
  mu_itr = muons->begin();
  mu_end = muons->end();
  for( ; mu_itr != mu_end; ++mu_itr ) {
    
    if((*mu_itr)->auxdecor<char>("passOR") != 1) continue;//if( !dec_passOR(**mu_itr) ) continue;
    xAOD::JetContainer::const_iterator jet_itr = jets->begin();
    xAOD::JetContainer::const_iterator jet_end = jets->end();
    for( ; jet_itr != jet_end; ++jet_itr ) {
      
      if((*jet_itr)->auxdecor<char>("passOR") != 1) continue;// if ( !dec_passOR( **jet_itr ) ) continue;
      TLorentzVector mu4vec = (*mu_itr)->p4();
      TLorentzVector jet4vec = (*jet_itr)->p4();
      //std::vector<int> nTrkVec;
      //(*jet_itr)->getAttribute(xAOD::JetAttribute::NumTrkPt500, nTrkVec);
      //int jet_nTrk = nTrkVec[0];
      if (mu4vec.DeltaR(jet4vec)<dRjetmu) {
	//if(doHarmonization && jet_nTrk<3){
	//std::cout << " Rejecting jet at (pT,eta,phi)=(" << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10)  << (*jet_itr)->pt() <<","<< (*jet_itr)->eta() <<"," << (*jet_itr)->phi() <<") with only nTrk=" << jet_nTrk
	//	    << " due to muon at (pT,eta,phi)=(" << (*mu_itr)->pt() <<","<< (*mu_itr)->eta() <<"," << (*mu_itr)->phi() <<")"<< std::endl ;
	//(*jet_itr)->auxdecor<char>("passOR") = 0 ; //dec_passOR( **jet_itr ) = 0;
	//}else{
	std::cout << " Rejecting muon at (eta,phi)=(" <<  std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10) << (*mu_itr)->eta() <<"," << (*mu_itr)->phi() <<") "
		  << " due to jet at (eta,phi)=(" << (*jet_itr)->eta() <<"," << (*jet_itr)->phi() <<")"<< std::endl ;
	(*mu_itr)->auxdecor<char>("passOR") = 0 ;//dec_passOR( **mu_itr ) = 0;
	//}
      }
    }
  }
  

  // Remove electrons and muons overlapping with each other 
  el_itr = electrons->begin();
  el_end = electrons->end();
  for( ; el_itr != el_end; ++el_itr ) {
    
    if((*el_itr)->auxdecor<char>("passOR") != 1) continue ;//if( !dec_passOR(**el_itr) ) continue;
    
    mu_itr = muons->begin();
    mu_end = muons->end();
    
    for( ; mu_itr != mu_end; ++mu_itr ) {
      
      if((*mu_itr)->auxdecor<char>("passOR") != 1) continue ;//if ( !dec_passOR( **mu_itr ) ) continue; 
      
      TLorentzVector el4vec = (*el_itr)->p4();
      TLorentzVector mu4vec = (*mu_itr)->p4();
      
      if (el4vec.DeltaR(mu4vec)<dRemu) {
	std::cout << " Rejecting both electron at (eta,phi)=(" << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10)  << (*el_itr)->eta() <<"," << (*el_itr)->phi() <<") "
		  << " and muon at (eta,phi)=(" << (*mu_itr)->eta() <<"," << (*mu_itr)->phi() <<")"<< std::endl ;
	(*el_itr)->auxdecor<char>("passOR") = 0 ;//dec_passOR( **el_itr ) = 0; 
	(*mu_itr)->auxdecor<char>("passOR") = 0 ;//dec_passOR( **mu_itr ) = 0;  
      }
    }
  }
  // Remove electrons overlapping with each other
  el_itr = electrons->begin();
  el_end = electrons->end();
  for( ; el_itr != el_end; ++el_itr ) {
    
    if((*el_itr)->auxdecor<char>("passOR") != 1) continue ;//if( !dec_passOR(**el_itr) ) continue;
    
    xAOD::TruthParticleContainer::const_iterator el2_itr = electrons->begin();
    xAOD::TruthParticleContainer::const_iterator el2_end = electrons->end();
    
    for( ; el2_itr != el2_end; ++el2_itr ) {
      
      if(el_itr == el2_itr) continue;
      if((*el2_itr)->auxdecor<char>("passOR") != 1) continue ;//if ( !dec_passOR( **el2_itr ) ) continue;
      
      TLorentzVector el4vec = (*el_itr)->p4();
      TLorentzVector el24vec = (*el2_itr)->p4();
      
      if (el4vec.DeltaR(el24vec)<dRee) {
	if((*el_itr)->pt() < (*el2_itr)->pt()){
	  std::cout << " Rejecting electron at (eta,phi)=(" << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10)  << (*el_itr)->eta() <<"," << (*el_itr)->phi() <<") "
		    << " and muon at (eta,phi)=(" << (*el2_itr)->eta() <<"," << (*el2_itr)->phi() <<")"<< std::endl ;
	  (*el_itr)->auxdecor<char>("passOR") = 0 ;//dec_passOR( **el_itr ) = 0;  
	}else{
	  std::cout << " Rejecting electron at (eta,phi)=(" << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10) << (*el2_itr)->eta() <<"," << (*el2_itr)->phi() <<") "
		    << " and muon at (eta,phi)=(" << (*el_itr)->eta() <<"," << (*el_itr)->phi() <<")"<< std::endl ;
	  (*el2_itr)->auxdecor<char>("passOR") = 0 ;//dec_passOR( **el2_itr ) = 0;   
	}
      }
    }
  }
  // Count number of objects after overlap removal
  int Nel=0;
  el_itr = electrons->begin();
  el_end = electrons->end();
  for( ; el_itr != el_end; ++el_itr ) {
    if((*el_itr)->auxdecor<char>("passOR") == 1) Nel++;//if(dec_passOR( **el_itr )) Nel++; 
  }
  
  int Nmu=0;
  mu_itr = muons->begin();
  mu_end = muons->end();
  for( ; mu_itr != mu_end; ++mu_itr ) {
    if((*mu_itr)->auxdecor<char>("passOR") == 1) Nmu++;//if(dec_passOR( **mu_itr )) Nmu++;  
  }
  
  int Njet=0;
  jet_itr = jets->begin();
  jet_end = jets->end();
  for( ; jet_itr != jet_end; ++jet_itr ) {
    if((*jet_itr)->auxdecor<char>("passOR") == 1) Njet++;//if(dec_passOR( **jet_itr )) Njet++;
  }
  //std::cout << " Before overlap removal: Nel=" << Nelin <<", Nmu="<< Nmuin <<", Njet=" << Njetin<< std::endl ;
  //std::cout << " After  overlap removal: Nel=" << Nel <<", Nmu="<< Nmu <<", Njet=" << Njet<< std::endl ;
  return true;
  
}

bool BuildTruthObjects::OverlapRemoval(const xAOD::TruthParticleContainer *electrons, const xAOD::TruthParticleContainer *muons, const xAOD::JetContainer *jets, const xAOD::TruthParticleContainer *photons, const bool doHarmonization, const double dRejet, const double dRjetmu, const double dRjete, double dRemu, double dRee, double dRphjet, double dReph, double dRmuph){

  xAOD::JetContainer::const_iterator jet_itr = jets->begin();
  xAOD::JetContainer::const_iterator jet_end = jets->end();
  for( ; jet_itr != jet_end; ++jet_itr ) {
    //bool jet_sel = (*jet_itr)->auxdecor<char>("baseline"); //bool jet_sel = dec_baseline(**jet_itr);
    bool jet_sel; 
    if( (*jet_itr)->pt() > 20000 && std::abs((*jet_itr)->eta()) < 2.8 )
      jet_sel = 1 ; 
    if(jet_sel)
      (*jet_itr)->auxdecor<char>("passOR") = 1 ;//dec_passOR( **jet_itr ) = 1;
    else
      (*jet_itr)->auxdecor<char>("passOR") = 0 ; //dec_passOR( **jet_itr ) = 0;
  }
  xAOD::TruthParticleContainer::const_iterator mu_itr = muons->begin();
  xAOD::TruthParticleContainer::const_iterator mu_end = muons->end();
  for( ; mu_itr != mu_end; ++mu_itr ) {
    bool mu_sel;
    //if(doHarmonization) mu_sel = (*mu_itr)->auxdecor<char>("signal"); //mu_sel = dec_signal(**mu_itr);      
    //else mu_sel = (*mu_itr)->auxdecor<char>("baseline"); //dec_baseline(**mu_itr); 
    if( (*mu_itr)->pt() > 10000 && std::abs((*mu_itr)->eta()) < 2.4 )
      mu_sel = 1 ; 
    if(mu_sel)
      (*mu_itr)->auxdecor<char>("passOR") = 1 ;  //dec_passOR( **mu_itr ) = 1;   
    else
      (*mu_itr)->auxdecor<char>("passOR") = 0 ;//dec_passOR( **mu_itr ) = 0;  
  }

  // remove jets overlapping with (baseline/signal) electrons
  xAOD::TruthParticleContainer::const_iterator el_itr = electrons->begin();
  xAOD::TruthParticleContainer::const_iterator el_end = electrons->end();
  for( ; el_itr != el_end; ++el_itr ) {
    bool el_sel;
    //if(doHarmonization) el_sel = (*el_itr)->auxdecor<char>("signal"); // el_sel = dec_signal(**el_itr);    
    //else el_sel = (*el_itr)->auxdecor<char>("baseline"); //dec_baseline(**el_itr); 
    if( (*el_itr)->pt() > 20000 && std::abs((*el_itr)->eta()) < 2.47 )
      el_sel = 1 ; 
    if( !el_sel ){
      (*el_itr)->auxdecor<char>("passOR") = 0 ;//dec_passOR( **el_itr ) = 0;    
      continue;
    }else{
      (*el_itr)->auxdecor<char>("passOR") = 1 ;//dec_passOR( **el_itr ) = 1; 
    }

    xAOD::JetContainer::const_iterator jet_itr = jets->begin();
    xAOD::JetContainer::const_iterator jet_end = jets->end();

    for( ; jet_itr != jet_end; ++jet_itr ) {
      if((*jet_itr)->auxdecor<char>("passOR") != 1) continue ;//if( !dec_passOR(**jet_itr) ) continue; 
      
      TLorentzVector el4vec = (*el_itr)->p4();
      TLorentzVector jet4vec = (*jet_itr)->p4();
      
      if (el4vec.DeltaR(jet4vec)<dRejet) {
	std::cout << " Rejecting jet at (eta,phi)=(" << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10)  << (*jet_itr)->eta() <<"," << (*jet_itr)->phi() <<") "
		  << " due to electron at (eta,phi)=(" << (*el_itr)->eta() <<"," << (*el_itr)->phi() <<")"<< std::endl ;
	(*jet_itr)->auxdecor<char>("passOR") = 0 ; //dec_passOR( **jet_itr ) = 0;
      }
    }
  } // END loop over electrons

  // remove jets overlapping with (baseline/signal) photons
  xAOD::TruthParticleContainer::const_iterator ph_itr = photons->begin();
  xAOD::TruthParticleContainer::const_iterator ph_end = photons->end();
  for( ; ph_itr != ph_end; ++ph_itr ) {
    bool ph_sel;
    //if(doHarmonization) ph_sel = (*ph_itr)->auxdecor<char>("signal"); //dec_signal(**ph_itr);   
    //else ph_sel = (*ph_itr)->auxdecor<char>("baseline"); //dec_baseline(**ph_itr);
    if( (*ph_itr)->pt() > 0 && std::abs((*ph_itr)->eta()) < 100 )
      ph_sel = 1 ; 
    if( !ph_sel ){
      (*ph_itr)->auxdecor<char>("passOR") = 0 ; //dec_passOR( **ph_itr ) = 0;    
      continue;
    }else{
      (*ph_itr)->auxdecor<char>("passOR") = 1 ; //dec_passOR( **ph_itr ) = 1;   
    }

    xAOD::JetContainer::const_iterator jet_itr = jets->begin();
    xAOD::JetContainer::const_iterator jet_end = jets->end();
    
    for( ; jet_itr != jet_end; ++jet_itr ) {
      if((*jet_itr)->auxdecor<char>("passOR") != 1) continue ;  //if( !dec_passOR(**jet_itr) ) continue;  
      
      TLorentzVector ph4vec = (*ph_itr)->p4();
      TLorentzVector jet4vec = (*jet_itr)->p4();
      
      if (ph4vec.DeltaR(jet4vec)<dRphjet) {
	std::cout << " Rejecting jet at (eta,phi)=(" << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10) << (*jet_itr)->eta() <<"," << (*jet_itr)->phi() <<") "
		  << " due to photon at (eta,phi)=(" << (*ph_itr)->eta() <<"," << (*ph_itr)->phi() <<")"<< std::endl ;
	(*jet_itr)->auxdecor<char>("passOR") = 0 ; //dec_passOR( **jet_itr ) = 0;       
      }
    }
  }// END loop over photons  
  
  // Remove electrons and muons overlapping with jets and photons   
  el_itr = electrons->begin();
  el_end = electrons->end();
  for( ; el_itr != el_end; ++el_itr ) {
    
    if((*el_itr)->auxdecor<char>("passOR") != 1) continue ;//if( !dec_passOR(**el_itr) ) continue; 
    
    xAOD::JetContainer::const_iterator jet_itr = jets->begin();
    xAOD::JetContainer::const_iterator jet_end = jets->end();
    
    for( ; jet_itr != jet_end; ++jet_itr ) {
      
      if((*jet_itr)->auxdecor<char>("passOR") != 1) continue ; //if ( !dec_passOR( **jet_itr ) ) continue;   
      
      TLorentzVector el4vec = (*el_itr)->p4();
      TLorentzVector jet4vec = (*jet_itr)->p4();
      
      if (el4vec.DeltaR(jet4vec)<dRjete) {
	std::cout <<  " Rejecting electron at (eta,phi)=(" << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10)  << (*el_itr)->eta() <<"," << (*el_itr)->phi() <<") "
		  << " due to jet at (eta,phi)=(" << (*jet_itr)->eta() <<"," << (*jet_itr)->phi() <<")"<< std::endl ;
	(*el_itr)->auxdecor<char>("passOR") = 0 ; //dec_passOR( **el_itr ) = 0;    
	
      }
    }
    
  }
  
  mu_itr = muons->begin();
  mu_end = muons->end();
  
  for( ; mu_itr != mu_end; ++mu_itr ) {
    
    if((*mu_itr)->auxdecor<char>("passOR") != 1) continue;//if( !dec_passOR(**mu_itr) ) continue;   
    
    xAOD::JetContainer::const_iterator jet_itr = jets->begin();
    xAOD::JetContainer::const_iterator jet_end = jets->end();
    
    for( ; jet_itr != jet_end; ++jet_itr ) {
      
      if((*jet_itr)->auxdecor<char>("passOR") != 1) continue; //if ( !dec_passOR( **jet_itr ) ) continue;                                                                                                                                        
      
      TLorentzVector mu4vec = (*mu_itr)->p4();
      TLorentzVector jet4vec = (*jet_itr)->p4();
      
      std::vector<int> nTrkVec;
      (*jet_itr)->getAttribute(xAOD::JetAttribute::NumTrkPt500, nTrkVec);
      int jet_nTrk = nTrkVec[0];
      
      if (mu4vec.DeltaR(jet4vec)<dRjetmu) {
	if(doHarmonization && jet_nTrk<3){
	  std::cout << " Rejecting jet at (pT,eta,phi)=(" << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10)  << (*jet_itr)->pt() <<","<< (*jet_itr)->eta() <<"," << (*jet_itr)->phi() <<") with only nTrk=" << jet_nTrk
		    << " due to muon at (pT,eta,phi)=(" << (*mu_itr)->pt() <<","<< (*mu_itr)->eta() <<"," << (*mu_itr)->phi() <<")"<< std::endl ;
	  (*jet_itr)->auxdecor<char>("passOR") = 0 ;//dec_passOR( **jet_itr ) = 0; 
	}else{
	  std::cout <<  " Rejecting muon at (eta,phi)=(" << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10)  << (*mu_itr)->eta() <<"," << (*mu_itr)->phi() <<") "
		    << " due to jet at (eta,phi)=(" << (*jet_itr)->eta() <<"," << (*jet_itr)->phi() <<")"<< std::endl ;
	  (*mu_itr)->auxdecor<char>("passOR") = 0 ; //dec_passOR( **mu_itr ) = 0; 
	}
      }
    }
  }
  
  // Remove electrons and muons overlapping with each other  
  el_itr = electrons->begin();
  el_end = electrons->end();
  for( ; el_itr != el_end; ++el_itr ) {
    
    if((*el_itr)->auxdecor<char>("passOR") != 1) continue ;//if( !dec_passOR(**el_itr) ) continue;  
    
    mu_itr = muons->begin();
    mu_end = muons->end();
    
    for( ; mu_itr != mu_end; ++mu_itr ) {
      if((*mu_itr)->auxdecor<char>("passOR") != 1) continue ; //if ( !dec_passOR( **mu_itr ) ) continue; 
      
      TLorentzVector el4vec = (*el_itr)->p4();
      TLorentzVector mu4vec = (*mu_itr)->p4();
      
      if (el4vec.DeltaR(mu4vec)<dRemu) {
	std::cout << " Rejecting both electron at (eta,phi)=(" << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10)  << (*el_itr)->eta() <<"," << (*el_itr)->phi() <<") "
		  << " and muon at (eta,phi)=(" << (*mu_itr)->eta() <<"," << (*mu_itr)->phi() <<")"<< std::endl ;
	(*el_itr)->auxdecor<char>("passOR") = 0 ;//dec_passOR( **el_itr ) = 0;
	(*mu_itr)->auxdecor<char>("passOR") = 0 ; //dec_passOR( **mu_itr ) = 0; 
      }
    }
  }
  
  // Remove electrons overlapping with each other
  el_itr = electrons->begin();
  el_end = electrons->end();
  for( ; el_itr != el_end; ++el_itr ) {
    
    if((*el_itr)->auxdecor<char>("passOR") != 1) continue ;  //if( !dec_passOR(**el_itr) ) continue; 
    
    xAOD::TruthParticleContainer::const_iterator el2_itr = electrons->begin();
    xAOD::TruthParticleContainer::const_iterator el2_end = electrons->end();
    
    for( ; el2_itr != el2_end; ++el2_itr ) {
      
      if(el_itr == el2_itr) continue;
      if((*el2_itr)->auxdecor<char>("passOR") != 1) continue ; //if ( !dec_passOR( **el2_itr ) ) continue;
      
      TLorentzVector el4vec = (*el_itr)->p4();
      TLorentzVector el24vec = (*el2_itr)->p4();
      
      if (el4vec.DeltaR(el24vec)<dRee) {
	if((*el_itr)->pt() < (*el2_itr)->pt()){
	  std::cout << " Rejecting electron at (eta,phi)=(" << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10) << (*el_itr)->eta() <<"," << (*el_itr)->phi() <<") "
		    << " and muon at (eta,phi)=(" << (*el2_itr)->eta() <<"," << (*el2_itr)->phi() <<")"<< std::endl ;
	  (*el_itr)->auxdecor<char>("passOR") = 0 ; //dec_passOR( **el_itr ) = 0;
	}else{
	  std::cout << " Rejecting electron at (eta,phi)=(" << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10)  << (*el2_itr)->eta() <<"," << (*el2_itr)->phi() <<") "
		    << " and muon at (eta,phi)=(" << (*el_itr)->eta() <<"," << (*el_itr)->phi() <<")"<< std::endl ;
	  (*el2_itr)->auxdecor<char>("passOR") = 0 ; //dec_passOR( **el2_itr ) = 0;
	}
      }
    }
  }
  
  // Remove photons if overlapping with electrons 
  
  ph_itr = photons->begin();
  ph_end = photons->end();
  for( ; ph_itr != ph_end; ++ph_itr ) {
    
    if((*ph_itr)->auxdecor<char>("passOR") != 1) continue ; //if( !dec_passOR( **ph_itr ) )
    continue;
    
    el_itr = electrons->begin();
    el_end = electrons->end();
    for( ; el_itr != el_end; ++el_itr ) {
      if((*el_itr)->auxdecor<char>("passOR") != 1) continue ; //if( !dec_passOR(**el_itr) ) continue;
      
      TLorentzVector el4vec = (*el_itr)->p4();
      TLorentzVector ph4vec = (*ph_itr)->p4();
      
      if (el4vec.DeltaR(ph4vec)<dReph) {
	std::cout << " Rejecting photon at (eta,phi)=(" << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10)  << (*ph_itr)->eta() <<"," << (*ph_itr)->phi() <<") "
		  << " due to electron at (eta,phi)=(" << (*el_itr)->eta() <<"," << (*el_itr)->phi() <<")"<< std::endl ;
	(*ph_itr)->auxdecor<char>("passOR") = 0 ; //dec_passOR( **ph_itr ) = 0;
      }
    }
  }
  
  // Remove photons if overlapping with muons
  ph_itr = photons->begin();
  ph_end = photons->end();
  for( ; ph_itr != ph_end; ++ph_itr ) {
    
    if((*ph_itr)->auxdecor<char>("passOR") != 1) continue ; //if( !dec_passOR( **ph_itr ) ) 
    continue;
    
    mu_itr = muons->begin();
    mu_end = muons->end();
    
    for( ; mu_itr != mu_end; ++mu_itr ) {
      if((*mu_itr)->auxdecor<char>("passOR") != 1) continue ; //if( !dec_passOR(**mu_itr) ) continue; 
      TLorentzVector mu4vec = (*mu_itr)->p4();
      TLorentzVector ph4vec = (*ph_itr)->p4();
      
      if (mu4vec.DeltaR(ph4vec)<dRmuph) {
	std::cout << " Rejecting photon at (eta,phi)=(" << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10) << (*ph_itr)->eta() <<"," << (*ph_itr)->phi() <<") "
		  << " due to muon at (eta,phi)=(" << (*mu_itr)->eta() <<"," << (*mu_itr)->phi() <<")"<< std::endl ;
	(*ph_itr)->auxdecor<char>("passOR") = 0 ; //dec_passOR( **ph_itr ) = 0;
      }
    }
  }
  
  // Count number of objects after overlap removal  
  int Nel=0;
  el_itr = electrons->begin();
  el_end = electrons->end();
  for( ; el_itr != el_end; ++el_itr ) {
    if((*el_itr)->auxdecor<char>("passOR") == 1) Nel++ ;//if(dec_passOR( **el_itr )) Nel++;
  }
  
  int Nmu=0;
  mu_itr = muons->begin();
  mu_end = muons->end();
  for( ; mu_itr != mu_end; ++mu_itr ) {
    if((*mu_itr)->auxdecor<char>("passOR") == 1) Nmu++ ;//if(dec_passOR( **mu_itr )) Nmu++;
  }
  
  int Njet=0;
  jet_itr = jets->begin();
  jet_end = jets->end();
  for( ; jet_itr != jet_end; ++jet_itr ) {
    if((*jet_itr)->auxdecor<char>("passOR") == 1) Njet++ ;//if(dec_passOR( **jet_itr )) Njet++; 
  }
  
  int Nph=0;
  ph_itr = photons->begin();
  ph_end = photons->end();
  for( ; ph_itr != ph_end; ++ph_itr ) {
    if((*ph_itr)->auxdecor<char>("passOR") == 1) Nph++ ;//if(dec_passOR( **ph_itr )) Nph++; 
  }
  
  std::cout << " After overlap removal: Nel=" << Nel <<", Nmu="<< Nmu <<", Njet=" << Njet << ", Nph=" << Nph << std::endl;
  return true;
  
}




ClassImp(BuildTruthObjects);

