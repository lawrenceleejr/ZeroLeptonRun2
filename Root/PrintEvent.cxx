
#include "ZeroLeptonRun2/PrintEvent.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TActiveStore.h"
#include "xAODRootAccess/TStore.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODJet/JetContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODMissingET/MissingETContainer.h"
#include "cafe/Config.h"
#include <stdexcept>
#include <iomanip>

PrintEvent::PrintEvent(const char *name)
  : cafe::Processor(name),
    m_jetKeys(),
    m_elKeys(),
    m_muKeys(),
    m_METKeys()
{
  cafe::Config config(name);
  m_jetKeys = config.getVString("JetKeys");
  m_elKeys = config.getVString("ElectronKeys");
  m_muKeys = config.getVString("MuonKeys");
  m_METKeys = config.getVString("METKeys");
}

bool PrintEvent::processEvent(xAOD::TEvent& event)
{
  //write_xAOD_event();
  xAOD::TStore* store = xAOD::TActiveStore::store();

  const xAOD::EventInfo* eventInfo = 0;
  if ( !event.retrieve( eventInfo, "EventInfo").isSuccess() ) return true; 

  bool isData = true;
  if(eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ){
    isData = false; 
  }
  uint32_t RunNumber = eventInfo->runNumber();
  unsigned long long EventNumber = eventInfo->eventNumber();
  uint32_t mc_channel_number = 0;
  if ( ! isData ) mc_channel_number = eventInfo->mcChannelNumber();
  out() << "---------------------------------------------------------------" << std::endl;
  out() << "---------------------------------------------------------------" << std::endl;
  out() << "Run " <<  RunNumber << " Event " << EventNumber << " isData " <<
    isData << " channel " << mc_channel_number << std::endl;
  if ( isData) out() << "StatusElement " <<  eventInfo->statusElement() << " extL1ID " << eventInfo->extendedLevel1ID()  << " L1Type " << eventInfo->level1TriggerType()  << std::endl;

  for ( std::size_t i = 0; i < m_jetKeys.size(); ++i ){
    const xAOD::JetContainer* jets = 0;
    std::string tag = m_jetKeys[i];
    bool inTS = inTStore(tag);
    if ( (inTS && store->retrieve(jets, tag).isSuccess()) ||
	 (!inTS && event.retrieve(jets, tag).isSuccess()) ) {
      out() << "------------- Jets : key = "<< tag  << m_jetKeys[i]<< std::endl;
      out() << "     pT       eta    phi      E         M   "  << std::endl;
      for ( xAOD::JetContainer::const_iterator it = jets->begin(); 
	    it != jets->end(); it++ ){
	out() << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10) << (*it)->pt() << " " << std::setw(6) << std::setprecision(3) << (*it)->eta() << " " << std::setw(6) << (*it)->phi() << " " << std::setw(10) << std::setprecision(1) << (*it)->e() << " " << std::setw(8) << (*it)->m() << std::endl;
      }
    }
  }

  for ( std::size_t i = 0; i < m_elKeys.size(); ++i ){
    const xAOD::ElectronContainer* electrons = 0;
    std::string tag = m_elKeys[i];
    bool inTS = inTStore(tag);
    if ( (inTS && store->retrieve(electrons, tag).isSuccess()) ||
	 (!inTS && event.retrieve(electrons, tag).isSuccess()) ){
      out() << "------------- Electrons : key = "<< tag << std::endl;
      out() << "     pT       eta    phi      E         M   "  << std::endl;
      for ( xAOD::ElectronContainer::const_iterator it = electrons->begin(); 
	    it != electrons->end(); it++ ){
	out() << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10) << (*it)->pt() << " " << std::setw(6) << std::setprecision(3) << (*it)->eta() << " " << std::setw(6) << (*it)->phi() << " " << std::setw(10) << std::setprecision(1) << (*it)->e() << " " << std::setw(6) << (*it)->m() << std::endl;
      }
    }
  }


  for ( std::size_t i = 0; i < m_muKeys.size(); ++i ){
    const xAOD::MuonContainer* muons = 0;
    std::string tag = m_muKeys[i];
    bool inTS = inTStore(tag);
    if ( (inTS && store->retrieve(muons, tag).isSuccess()) ||
	 (!inTS && event.retrieve(muons, tag).isSuccess()) ){
      out() << "------------- Muons : key = "<< tag << std::endl;
      out() << "     pT       eta    phi      E         M    type  "  << std::endl;
      for ( xAOD::MuonContainer::const_iterator it = muons->begin(); 
	    it != muons->end(); it++ ){
	out() << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10) << (*it)->pt() << " " << std::setw(6) << std::setprecision(3) << (*it)->eta() << " " << std::setw(6) << (*it)->phi() << " " << std::setw(10) << std::setprecision(1) << (*it)->e() << " " << std::setw(6) << (*it)->m() << " " << std::setw(6) << (int) (*it)->muonType() << std::endl;
      }
    }
  }

  for ( std::size_t i = 0; i < m_METKeys.size(); ++i ){
    const xAOD::MissingETContainer* met = 0;
    std::string tag = m_METKeys[i];
    bool inTS = inTStore(tag);
    if ( (inTS && store->retrieve(met, tag).isSuccess()) ||
	 (!inTS && event.retrieve(met, tag).isSuccess()) ){
      out() << "------------- MET : key = "<< tag << std::endl;
      out() << "      Term           px          py        MET     "  << std::endl;
      for ( xAOD::MissingETContainer::const_iterator it = met->begin();
	    it != met->end(); it++ ){
	out() << std::setw(16) << (*it)->name() << " " << std::setiosflags(std::ios::fixed) << std::setprecision(1) << std::setw(10) << (*it)->mpx() << " " << std::setw(10) << (*it)->mpy() << " " << std::setw(10) << (*it)->met() << std::endl;
      }
    }
  }

  return true;
}

bool PrintEvent::inTStore(std::string& tag)
{
  std::size_t pos = tag.find('<');
  if ( pos == std::string::npos || pos < 2 ) return false;
  if ( tag.substr(0,2) != "TS" ) return false;
  tag = tag.substr(pos+1);
  return true;
}

ClassImp(PrintEvent);

