#ifndef ZeroLeptonRun2_BuildSUSYObjects_H_
#define ZeroLeptonRun2_BuildSUSYObjects_H_

#include "cafe/Processor.h"
#include "ZeroLeptonRun2/ZeroLeptonRunPeriod.h"
#include "AsgTools/ToolHandle.h"
#include "TauAnalysisTools/TauSelectionTool.h"
#include "TauAnalysisTools/TauEfficiencyCorrectionsTool.h"
#include "TauAnalysisTools/TauTruthMatchingTool.h"

namespace CP{
	class SystematicSet;
}
namespace ST{
  class SUSYObjDef_xAOD;
  class SystInfo;
}

#include <string>
#include <vector>

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
  void initSUSYTools();
  void fillTriggerInfo     (xAOD::TEvent& event) const;
  void fillTriggerJetAndMET(xAOD::TEvent& event) const;
  ST::SUSYObjDef_xAOD* m_SUSYObjTool;

  bool m_IsData;
  bool m_Is25ns;
  bool m_IsAtlfast;
  bool m_UseSmearedJets;
  bool m_DoSystematics;
  bool m_PhotonInOR;


  ToolHandle<TauAnalysisTools::ITauSelectionTool> m_tauSelTool;
  ToolHandle<TauAnalysisTools::ITauEfficiencyCorrectionsTool> m_tauEffTool;
  ToolHandle<TauAnalysisTools::ITauTruthMatchingTool> m_tauTruthMatchTool;
  std::vector<CP::SystematicSet> m_tauEffSystSetList;

  std::string m_jetkey;
  std::string m_taukey;
  std::string m_suffix;
  std::string m_suffixRecl;
  ZeroLeptonRunPeriod m_period;
  ZeroLeptonDerivationTag m_derivationTag;
  int m_JESNuisanceParameterSet;
  std::string m_ECKey;
  std::string m_PCKey;

  std::vector<ST::SystInfo> m_SystInfoList;
  std::vector<std::string> m_SystMatch;
  bool m_buildTriggerJetAndMET;
public:
  ClassDef(BuildSUSYObjects,0);
};

#endif // ZeroLeptonRun2_BuildSUSYObjects_H_
