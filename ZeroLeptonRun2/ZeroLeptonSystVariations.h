#ifndef ZeroLeptonRun2_ZeroLeptonSystVariations_H_
#define ZeroLeptonRun2_ZeroLeptonSystVariations_H_

#include "cafe/Processor.h"
#include "ZeroLeptonRun2/Counter.h"

#include <string>
#include <list>
#include <vector>
#include <memory>

class TFile;

//--------------------------------------------------------------------------
// ZeroLeptonSystVariations
//
// Works a bit like a cafe::Controller except that the child
// Processors are called for each systematics
//--------------------------------------------------------------------------
class ZeroLeptonSystVariations : public cafe::Processor {
public:
  ZeroLeptonSystVariations(const char *name);
  ~ZeroLeptonSystVariations();
  void begin();
  bool processEvent(xAOD::TEvent& event);
  void finish();
  virtual void inputFileOpened(TFile *file);
  virtual void inputFileClosing(TFile *file);

private:
  Counter* m_counter;

  std::list<cafe::Processor*> m_processors;
  bool add(const std::list<cafe::Processor*>& procs);
  bool add(cafe::Processor *proc);

public:
  ClassDef(ZeroLeptonSystVariations,0);
};

#endif // ZeroLeptonRun2_ZeroLeptonSystVariations_H_
