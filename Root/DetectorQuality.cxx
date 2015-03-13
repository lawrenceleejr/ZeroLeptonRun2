#include "ZeroLeptonRun2/DetectorQuality.h"

#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TActiveStore.h"
#include "xAODRootAccess/TStore.h"
#include "xAODRootAccess/tools/ReturnCheck.h"
#include "xAODEventInfo/EventInfo.h"
#include "cafe/Config.h"

#include <iostream>
#include <string>

DetectorQuality::DetectorQuality(const char *name)
  : cafe::Processor(name), m_period(INVALID)
{
  cafe::Config config(name);
  m_period = periodFromString(config.get("Period","p8tev"));
}

bool DetectorQuality::processEvent(xAOD::TEvent& event)
{
  const xAOD::EventInfo* eventInfo = 0;
  if( ! event.retrieve( eventInfo, "EventInfo").isSuccess() ) return true;

  //std::cout << "Event flag : LAr " << eventInfo->errorState(xAOD::EventInfo::LAr) << " Tile " << eventInfo->errorState(xAOD::EventInfo::Tile)  << " Core " << eventInfo->eventFlags(xAOD::EventInfo::Core)  << std::endl;

  bool* badDetectorQuality = new bool;
  *badDetectorQuality = eventInfo->errorState(xAOD::EventInfo::LAr) != xAOD::EventInfo::NotSet;
  if ( m_period == p8tev )  *badDetectorQuality = *badDetectorQuality && eventInfo->errorState(xAOD::EventInfo::Tile)==xAOD::EventInfo::Error && ((eventInfo->eventFlags(xAOD::EventInfo::Core) & 0x40000) != 0);

  xAOD::TStore* store = xAOD::TActiveStore::store();
  RETURN_CHECK("DetectorQuality::processEvent",store->record<bool>(badDetectorQuality,"badDetectorQuality"));

  return true;
}

ClassImp(DetectorQuality);

