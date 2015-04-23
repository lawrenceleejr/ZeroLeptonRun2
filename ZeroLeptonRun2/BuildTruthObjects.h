#ifndef ZeroLeptonRun2_BuildTruthObjects_H_
#define ZeroLeptonRun2_BuildTruthObjects_H_

#include "cafe/Processor.h"

#include <string>

//---------------------------------------------------------------------
// BuildTruthObjects use the collections of physics objects in the event
// The  output containers are stored in the transient store with keys "TruthJets", 
// "TruthMuons", "TruthElectrons", "TruthPhotons", "TruthTaus"
//
//---------------------------------------------------------------------
class BuildTruthObjects : public cafe::Processor {
public:
  BuildTruthObjects(const char *name);
  ~BuildTruthObjects();
  bool processEvent(xAOD::TEvent& event);

private:
  bool m_IsData;

  std::string m_jetkey; 
  std::string m_suffix;
public:
  ClassDef(BuildTruthObjects,0);
};

#endif // ZeroLeptonRun2_BuildTruthObjects_H_
