

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
    m_suffix()
{
  cafe::Config config(name);
  m_IsData = config.get("IsData",false);
  m_jetkey = config.get("JetCollectionKey","xxxx");
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
  if ( ! store->record(truthmuons2.first,"TruthMuons").isSuccess() ) { 
    throw std::runtime_error("Could not store truthMuons");
  }
  //if ( ! store->record<xAOD::ShallowAuxContainer>(truthmuons2.second,"TruthMuonsAux.").isSuccess()) { //BEFORE
  if ( ! store->record(truthmuons2.second,"TruthMuonsAux.").isSuccess()) {
    throw std::runtime_error("Could not store truthMuonsAux.");
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
  
  if ( ! store->record(truthphotons2.first,"TruthPhotons").isSuccess() ) {
    throw std::runtime_error("Could not store TruthPhotons");
  }
  if ( ! store->record(truthphotons2.second,"TruthPhotonsAux.").isSuccess()) {
    throw std::runtime_error("Could not store TruthPhotonsAux.");
  }


  
  // Taus
  const xAOD::TruthParticleContainer* truthtaus = 0 ;

  if (!event.retrieve(truthtaus,"TruthTaus").isSuccess()){
      throw std::runtime_error("Could not retrieve truth particles with key TruthTaus");
    }

  std::pair< xAOD::TruthParticleContainer*, xAOD::ShallowAuxContainer* > truthtaus2 = xAOD::shallowCopyContainer( *truthtaus );
  
  if ( ! store->record(truthtaus2.first,"TruthTaus").isSuccess() ) {
    throw std::runtime_error("Could not store TruthTaus");
  }
  if ( ! store->record(truthtaus2.second,"TruthTausAux.").isSuccess()) {
    throw std::runtime_error("Could not store TruthTausAux.");
  }
  
  

  // Overlap removal 
  //  if ( ! m_SUSYObjTool->OverlapRemoval(truthelectrons2.first, truthmuons2.first, outputjets, false, 0.2, 0.4, 0.4, 0.01, 0.05).isSuccess() ) throw std::runtime_error("Error in OverlapRemoval");

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

    throw std::runtime_error("Could not store TruthMET");
  }  

return true;
}


ClassImp(BuildTruthObjects);

