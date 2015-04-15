#ifndef ZeroLeptonRun2_BuildSUSYObjects_H_
#define ZeroLeptonRun2_BuildSUSYObjects_H_

#include "cafe/Processor.h"
#include "ZeroLeptonRun2/ZeroLeptonRunPeriod.h"
namespace ST{
  class SUSYObjDef_xAOD;
}

#include <string>

//---------------------------------------------------------------------
// BuildSUSYObjects use the collections of physics objects in the event
// and pass then through SUSYTools that applies calibration and add
// properties (as "decorations" to the objects). The output containers
// are stored in the transient store with keys "SUSYJets", "SUSYMuons",
// "SUSYElectrons", "SUSYPhotons", "SUSYTaus"
//
// If UseSmearedJets is true this is supposed to be called from within
// a loop on smeared iterations from an event.
//---------------------------------------------------------------------
class BuildSUSYObjects : public cafe::Processor {
public:
  BuildSUSYObjects(const char *name);
  ~BuildSUSYObjects();
  bool processEvent(xAOD::TEvent& event);

private:

  ST::SUSYObjDef_xAOD* m_SUSYObjTool;
  bool m_IsData;
  bool m_IsAtlfast;
  bool m_UseSmearedJets;
  bool m_UseSystematics;
  bool m_PhotonInOR;

  std::string m_jetkey; 
  std::string m_suffix;
  ZeroLeptonRunPeriod m_period;
  ZeroLeptonDerivationTag m_derivationTag;
  int m_JESNuisanceParameterSet;

public:
  ClassDef(BuildSUSYObjects,0);
};

#endif // ZeroLeptonRun2_BuildSUSYObjects_H_
