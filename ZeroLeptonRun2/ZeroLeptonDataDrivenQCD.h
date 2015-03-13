#ifndef ZeroLeptonRun2_ZeroLeptonDataDrivenQCD_H_
#define ZeroLeptonRun2_ZeroLeptonDataDrivenQCD_H_

#include "cafe/Processor.h"
#include "ZeroLeptonRun2/Counter.h"
#include "ZeroLeptonRun2/QCDSeedSelection.h"
#include "ZeroLeptonRun2/PhysObjProxyFiller.h"
#include "ZeroLeptonRun2/PhysObjProxyUtils.h"
#include "ZeroLeptonRun2/ZeroLeptonUtils.h"
#include "ZeroLeptonRun2/ZeroLeptonCutVal.h"

#include <string>
#include <vector>
#include <list>
#include <memory>

namespace SUSY {
  class JetMCSmearingTool;
}
class TH2F;
class TFile;

//--------------------------------------------------------------------------
// ZeroLeptonDataDrivenQCD selectes QCD data events based on trigger
// and balance requirements, applies quality cut and save these "seed events"
// in xAOD format if saveSeedEvents is true, else it uses the JetSmearing
// package to smear jets, recalculates MET and save the event.
//
// Works a bit like a cafe::Controller except that the child
// Processors are called for each smearing draw.
//--------------------------------------------------------------------------
class ZeroLeptonDataDrivenQCD : public cafe::Processor {
public:
  ZeroLeptonDataDrivenQCD(const char *name);
  ~ZeroLeptonDataDrivenQCD();
  void begin();
  bool processEvent(xAOD::TEvent& event);
  void finish();
  virtual void inputFileOpened(TFile *file);
  virtual void inputFileClosing(TFile *file);

private:
  void InitialiseSmearing();
  void GetResponseMaps(TFile* smearFnFile,TH2F*& nominal, TH2F*& tailHigh, TH2F*& tailLow);

  Counter* m_counter;
  int m_seed;
  bool m_saveSeedEvents;
  std::string m_jetkey; 

  QCDSeedSelection m_QCDSeedSelector;
  PhysObjProxyFiller m_physobjsFiller;
  PhysObjProxyUtils m_proxyUtils;
  ZeroLeptonUtils m_ZLUtils;
  SUSY::JetMCSmearingTool* m_smearingTool;
  TFile *m_smearFnFile;
  TFile *m_smearFnFile2;
  std::string m_GaussianCoreSmearingType;
  std::vector<std::string> m_containers;

  ZeroLeptonCutVal m_cutVal;

  std::list<cafe::Processor*> m_processors;
  bool add(const std::list<cafe::Processor*>& procs);
  bool add(cafe::Processor *proc);

public:
  ClassDef(ZeroLeptonDataDrivenQCD,0);
};

#endif // ZeroLeptonRun2_ZeroLeptonDataDrivenQCD_H_
