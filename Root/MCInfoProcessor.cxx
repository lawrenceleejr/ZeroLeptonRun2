
#include "ZeroLeptonRun2/MCInfoProcessor.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TActiveStore.h"
#include "xAODRootAccess/TStore.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthVertex.h"
#include "cafe/Config.h"
#include "SUSYTools/SUSYCrossSection.h"

#include <vector>
#include <stdexcept>

MCInfoProcessor::MCInfoProcessor(const char *name): 
  cafe::Processor(name), m_truthPKey(), m_isSignal(false), m_mcDB(0)
{
  cafe::Config config(name);
  m_isSignal = config.get("IsSignal",false);
  m_truthPKey = config.get("TruthParticleContainerKey","xxxx");
  std::string mcDBFile = config.get("MCDBFile","SUSYTools/data/mc12_8TeV/");
  bool mcDBextended = config.get("MCDBExtended",false);
  m_mcDB = new SUSY::CrossSectionDB(mcDBFile,mcDBextended);
}



bool MCInfoProcessor::processEvent(xAOD::TEvent& event)
{
  static std::vector<float> normWeight;
  normWeight.resize(3,0.);

  const xAOD::EventInfo* eventInfo = 0;
  if ( ! event.retrieve( eventInfo, "EventInfo").isSuccess() ) throw std::runtime_error("Could not retrieve EventInfo");
  //uint32_t mc_channel_number = eventInfo->mcChannelNumber();

  unsigned int* finalstate = new unsigned int;
  *finalstate = 0;
  if ( m_isSignal ) *finalstate = hardProcess(event);
  //out() << " MC harprocess code " << *finalstate << std::endl; 

  //normWeights(event,normWeight,mc_channel_number,*finalstate);
  //out() << " MC event weight for id " << mc_channel_number << " is " << normWeight[0] << std::endl;

  xAOD::TStore* store = xAOD::TActiveStore::store();
  if ( !store->record<unsigned int>(finalstate,"HardProcess").isSuccess()){
    throw std::runtime_error("Could not store HarProcess");
  }

  return true;
}

unsigned int MCInfoProcessor::hardProcess(xAOD::TEvent& event) const
{
  const xAOD::TruthParticleContainer* mcparticles = 0;
  if ( ! event.retrieve(mcparticles, m_truthPKey).isSuccess() ) throw std::runtime_error("MCInfoProcessor::hardProcess :: Could not retrieve TruthParticleContainer with key "+m_truthPKey);
  const xAOD::TruthParticle*  firstSUSY = 0;
  const xAOD::TruthParticle*  secondSUSY = 0;
  for ( xAOD::TruthParticleContainer::const_iterator it = mcparticles->begin();
	it != mcparticles->end(); ++it ){
    int id = (*it)->absPdgId();
    if( (id>1000000 && id<1000007) || // squarkL
	(id>1000010 && id<1000017) || // sleptonL
	(id>2000000 && id<2000007) || // squarkR
	(id>2000010 && id<2000017) || // sleptonR
	(id>1000020 && id<1000040)    // gauginos
	) {
      //std::cout << " SUSY particle found " << id << " " <<  ( (*it)->hasProdVtx() ? (*it)->prodVtx()->incomingParticle(0)->absPdgId() : 0 ) << std::endl;
      if ( ((*it)->hasProdVtx() && (*it)->prodVtx()->incomingParticle(0)->absPdgId() < 1000000 )   // SUSY with SM parents
	   || !(*it)->hasProdVtx() // ugly but happens in Herwig++
	   ) {
	if ( !firstSUSY ) { firstSUSY = *it; }
	else { secondSUSY = *it; break;}
      }
    }
  }
  if ( !secondSUSY ) throw std::runtime_error("MCInfoProcessor::hardProcess :: Could not find two primary SUSY particles");

  // now check that the particles found are not propagators (meaning that there is only a child and has different pdgId)
  if ( firstSUSY->hasDecayVtx() && 
       firstSUSY->decayVtx()->nOutgoingParticles()==1 && 
       firstSUSY->pdgId() != firstSUSY->decayVtx()->outgoingParticle(0)->pdgId() ){
    firstSUSY = firstSUSY->decayVtx()->outgoingParticle(0);
  }
  if ( secondSUSY->hasDecayVtx() && 
       secondSUSY->decayVtx()->nOutgoingParticles()==1 && 
       secondSUSY->pdgId() != secondSUSY->decayVtx()->outgoingParticle(0)->pdgId() ){
    secondSUSY = secondSUSY->decayVtx()->outgoingParticle(0);
  }

  return SUSY::finalState(firstSUSY->pdgId(), secondSUSY->pdgId());

}

void MCInfoProcessor::normWeights(xAOD::TEvent& event, 
				  std::vector<float>& normWeight, 
				  uint32_t mc_channel_number,
				  unsigned int finalstate)
{
  static uint32_t lastChannel = 0;
  if ( m_isSignal ) {
    float rel_uncertainty = m_mcDB->rel_uncertainty(mc_channel_number,finalstate);
    float sumW = m_mcDB->sumweight(mc_channel_number,finalstate);
    float xsec = m_mcDB->xsectTimesEff(mc_channel_number,finalstate);
    out() << " process cross-sec info for final state " << finalstate << xsec << " " << sumW << " "<< rel_uncertainty << std::endl;
    if (sumW > 0.f ) {
      normWeight[0]= xsec/sumW;
      normWeight[1]= normWeight[0]*(1+ rel_uncertainty);
      normWeight[2]= normWeight[0]*(1- rel_uncertainty);
    }
  }
  else {
    if ( mc_channel_number != lastChannel ) { // avoid querying DB for every event
      lastChannel = mc_channel_number;
      float rel_uncertainty = m_mcDB->rel_uncertainty(mc_channel_number);
      float sumW = m_mcDB->sumweight(mc_channel_number);
      float xsec = m_mcDB->xsectTimesEff(mc_channel_number);
      out() << " process cross-sec info " << xsec << " " << sumW << " "<< rel_uncertainty << std::endl;
      if (sumW > 0.f ) {
	normWeight[0]= xsec/sumW;
	normWeight[1]= normWeight[0]*(1+ rel_uncertainty);
	normWeight[2]= normWeight[0]*(1- rel_uncertainty);
      }
    }
  }
}



ClassImp(MCInfoProcessor);

