
#include "ZeroLeptonRun2/GRLProcessor.h"
#include "GoodRunsLists/GoodRunsListSelectionTool.h"
#include "GoodRunsLists/TGRLCollection.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TActiveStore.h"
#include "xAODRootAccess/TStore.h"
#include "xAODRootAccess/tools/ReturnCheck.h"
#include "xAODEventInfo/EventInfo.h"
#include "cafe/Config.h"

#include <vector>
#include <string>

GRLProcessor::GRLProcessor(const char *name)
  : cafe::Processor(name), m_GRLtool(0)
{
  cafe::Config config(name);
  std::string GRLXMLFile = config.get("GRLFile","");
  bool verbose = config.get("verbose",false);
  m_GRLtool = new GoodRunsListSelectionTool("cafeGRL");
  std::vector<std::string> GRLXMLs;
  GRLXMLs.push_back(GRLXMLFile);
  m_GRLtool->setProperty("GoodRunsListVec",GRLXMLs);
  m_GRLtool->setProperty("VerboseDetStatus",verbose);
  m_GRLtool->initialize();
  if ( verbose ) m_GRLtool->getGRLCollection()->Summary();
}
GRLProcessor::~GRLProcessor()
{
  if (  m_GRLtool ) delete  m_GRLtool;
}

bool GRLProcessor::processEvent(xAOD::TEvent& event)
{
  const xAOD::EventInfo* eventInfo = 0;
  if( ! event.retrieve( eventInfo, "EventInfo").isSuccess() ) return true;

  bool* passGRL = new bool(true);
  *passGRL = m_GRLtool->passRunLB(*eventInfo);

  xAOD::TStore* store = xAOD::TActiveStore::store();
  RETURN_CHECK("GRLProcessor::processEvent",store->record<bool>(passGRL,"passGRL"));


  return true;
}

ClassImp(GRLProcessor);

