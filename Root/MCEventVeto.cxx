
#include "ZeroLeptonRun2/MCEventVeto.h"
#include "ZeroLeptonRun2/MCEventVetoHelper.h"

#include "cafe/Config.h"

#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TActiveStore.h"
#include "xAODRootAccess/TStore.h"
#include "xAODRootAccess/tools/ReturnCheck.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/TruthEvent.h"
#include "xAODJet/JetContainer.h"
#include "xAODMissingET/MissingETContainer.h"

#include <stdexcept>

MCEventVeto::MCEventVeto(const char *name)
  : cafe::Processor(name), m_period(INVALID), m_TPCKey("")
{
  cafe::Config config(name);
  m_period = periodFromString(config.get("Period","p13tev"));
  if ( m_period == p7tev ) throw(std::domain_error("MCEventVeto does not support the 7tev run period"));
  if ( m_period == INVALID ) throw(std::domain_error("MCEventVeto: invalid run period specified"));

  m_TPCKey = config.get("TruthParticleContainerKey","TruthParticle");
}

bool MCEventVeto::processEvent(xAOD::TEvent& event)
{

  // Access all needed containers
  const xAOD::EventInfo* eventInfo = 0;
  RETURN_CHECK("MCEventVeto::processEvent", event.retrieve( eventInfo, "EventInfo") ); 

  bool isMC = false;
  //std::cout << " isMC " << eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) << std::endl;
  if(eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ){
    isMC = true; 
  }   
  if ( ! isMC ) return true;

  //const xAOD::TruthEventContainer* truthEventC = 0;
  //RETURN_CHECK("MCEventVeto::processEvent", event.retrieve( truthEventC, "TruthEvent" ) );

  const xAOD::TruthParticleContainer* truthPC = 0;
  RETURN_CHECK("MCEventVeto::processEvent", event.retrieve( truthPC, m_TPCKey ) );

  //const xAOD::JetContainer* jets = 0;
  //RETURN_CHECK("MCEventVeto::processEvent", event.retrieve( jets, "AntiKt4TruthJets" ) );

  const xAOD::MissingETContainer* metcontainer = 0;
  RETURN_CHECK("MCEventVeto::processEvent", event.retrieve( metcontainer, "MET_Truth" ));


  uint32_t mc_channel_number = eventInfo->mcChannelNumber();
  //uint32_t mc_channel_number = 147774;
  //bool isHighPtDijet = MCEventVetoHelper::isHighPtDijet(jets);
  //bool isHighPtJetMET = MCEventVetoHelper::isHighPtJetMET(mc_channel_number,jets,metcontainer);
  //unsigned int vetoQEDFSR = MCEventVetoHelper::vetoQEDFSR(mc_channel_number,truthPC);

  bool* mcaccept = new bool(true);
  unsigned int* veto = new unsigned int(0);
  if ( m_period == p8tev ) {
    *mcaccept = MCEventVetoHelper::mc12accept(*veto, mc_channel_number, truthPC, metcontainer);
  }
  else if ( m_period == p13tev ) {
    *mcaccept = MCEventVetoHelper::mc14accept(*veto, mc_channel_number, truthPC, metcontainer);
  }

  xAOD::TStore* store = xAOD::TActiveStore::store();
  RETURN_CHECK("MCEventVeto::processEvent",store->record<bool>(mcaccept,"mcAccept"));
  RETURN_CHECK("MCEventVeto::processEvent",store->record<unsigned int>(veto,"mcVetoCode"));

  return true;
}


ClassImp(MCEventVeto);

