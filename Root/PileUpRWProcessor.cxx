
#include "ZeroLeptonRun2/PileUpRWProcessor.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TActiveStore.h"
#include "xAODRootAccess/TStore.h"
#include "xAODRootAccess/tools/ReturnCheck.h"
#include "xAODEventInfo/EventInfo.h"
#include "cafe/Config.h"
#include "PileupReweighting/TPileupReweighting.h"
#include "PileupReweighting/PileupReweightingTool.h"

#include <string>
#include <vector>
#include <stdexcept>

PileUpRWProcessor::PileUpRWProcessor(const char *name)
    : cafe::Processor(name)
{
  cafe::Config config(name);
  m_forcedRunNumber = config.get("forcedRunNumber",-1);

  std::string PileUpMCFileName = config.get("PileUpMCFileName","PileupReweighting/mc14v1_defaults.prw.root");
  std::string PileUpDataFileName = config.get("PileUpDataFileName","SUSYTools/susy_data12_avgintperbx.root");

  std::vector<std::string> prwFiles;
  prwFiles.push_back(PileUpMCFileName);
  std::vector<std::string> lumicalcFiles;
  lumicalcFiles.push_back(PileUpDataFileName);

  m_PileupTool_CENTRAL = std::auto_ptr<CP::PileupReweightingTool>(new CP::PileupReweightingTool("PileUpReweightingTool_CENTRAL"));
  if ( !m_PileupTool_CENTRAL->setProperty("Prefix","CENTRAL_").isSuccess()) throw std::runtime_error("Could not initialise PileupReweightingTool");
  if ( !m_PileupTool_CENTRAL->setProperty("Input","EventInfo").isSuccess()) throw std::runtime_error("Could not initialise PileupReweightingTool");
  if ( !m_PileupTool_CENTRAL->setProperty("ConfigFiles",prwFiles).isSuccess()) throw std::runtime_error("Could not initialise PileupReweightingTool");
  if ( !m_PileupTool_CENTRAL->setProperty("LumiCalcFiles",lumicalcFiles).isSuccess()) throw std::runtime_error("Could not initialise PileupReweightingTool");
  if ( !m_PileupTool_CENTRAL->initialize().isSuccess()) throw std::runtime_error("Could not initialise PileupReweightingTool");

  m_PileupTool_UP = std::auto_ptr<CP::PileupReweightingTool>(new CP::PileupReweightingTool("PileUpReweightingTool_UP"));
  if ( !m_PileupTool_UP->setProperty("Prefix","UP_").isSuccess()) throw std::runtime_error("Could not initialise PileupReweightingTool");
  if ( !m_PileupTool_UP->setProperty("DataScaleFactor",1.1).isSuccess()) throw std::runtime_error("Could not initialise PileupReweightingTool");
  if ( !m_PileupTool_UP->setProperty("Input","EventInfo").isSuccess()) throw std::runtime_error("Could not initialise PileupReweightingTool");
  if ( !m_PileupTool_UP->setProperty("ConfigFiles",prwFiles).isSuccess()) throw std::runtime_error("Could not initialise PileupReweightingTool");
  if ( !m_PileupTool_UP->setProperty("LumiCalcFiles",lumicalcFiles).isSuccess()) throw std::runtime_error("Could not initialise PileupReweightingTool");
  if ( !m_PileupTool_UP->initialize().isSuccess()) throw std::runtime_error("Could not initialise PileupReweightingTool");

  m_PileupTool_DOWN = std::auto_ptr<CP::PileupReweightingTool>(new CP::PileupReweightingTool("PileUpReweightingTool_DOWN"));
  if ( !m_PileupTool_DOWN->setProperty("Prefix","DOWN_").isSuccess()) throw std::runtime_error("Could not initialise PileupReweightingTool");
  if ( !m_PileupTool_DOWN->setProperty("DataScaleFactor",0.9).isSuccess()) throw std::runtime_error("Could not initialise PileupReweightingTool");
  if ( !m_PileupTool_DOWN->setProperty("Input","EventInfo").isSuccess()) throw std::runtime_error("Could not initialise PileupReweightingTool");
  if ( !m_PileupTool_DOWN->setProperty("ConfigFiles",prwFiles).isSuccess()) throw std::runtime_error("Could not initialise PileupReweightingTool");
  if ( !m_PileupTool_DOWN->setProperty("LumiCalcFiles",lumicalcFiles).isSuccess()) throw std::runtime_error("Could not initialise PileupReweightingTool");
  if ( !m_PileupTool_DOWN->initialize().isSuccess()) throw std::runtime_error("Could not initialise PileupReweightingTool");
}

bool PileUpRWProcessor::processEvent(xAOD::TEvent& event)
{
  std::vector<float>* pileupWeight = new std::vector<float>(3,1.);

  const xAOD::EventInfo* eventInfo = 0;
  if( ! event.retrieve( eventInfo, "EventInfo").isSuccess() ) throw std::runtime_error("Could not retrieve EventInfo");

  int runnumber = -1;
  if ( m_forcedRunNumber > 0 ) runnumber = m_forcedRunNumber; 
  else runnumber = (int)eventInfo->runNumber();

  // no reweighting for mc14_13TeV nor in mc15_13TeV
  if ( runnumber != 222222  && runnumber != 222250  && runnumber != 222510) {
    if ( !m_PileupTool_CENTRAL->execute().isSuccess()) throw std::runtime_error("Could not execute PileupReweightingTool");
    if ( !m_PileupTool_UP->execute().isSuccess()) throw std::runtime_error("Could not execute PileupReweightingTool");
    if ( !m_PileupTool_DOWN->execute().isSuccess()) throw std::runtime_error("Could not execute PileupReweightingTool");
    (*pileupWeight)[0] = eventInfo->auxdata<double>("CENTRAL_PileupWeight");
    (*pileupWeight)[1] = eventInfo->auxdata<double>("UP_PileupWeight");
    (*pileupWeight)[2] = eventInfo->auxdata<double>("DOWN_PileupWeight");
  }

  //std::cout << "Pileup weights " << pileupWeight->at(0) << " " << pileupWeight->at(1) << " " << pileupWeight->at(2) << std::endl;

  xAOD::TStore* store = xAOD::TActiveStore::store();
  RETURN_CHECK("PileUpRWProcessor::processEvent",store->record(pileupWeight,"pileupWeights"));


  return true;
}

ClassImp(PileUpRWProcessor);

