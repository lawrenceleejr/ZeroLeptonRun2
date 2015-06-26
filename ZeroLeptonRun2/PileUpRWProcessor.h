#ifndef ZeroLeptonRun2_PileUpRWProcessor_H_
#define ZeroLeptonRun2_PileUpRWProcessor_H_

namespace Root{
  class TPileupReweighting;
}
namespace CP{
  class PileupReweightingTool;
}

#include "cafe/Processor.h"
#include <memory>


//FIXME : at some point one may use the  PileupReweightingTool

class PileUpRWProcessor : public cafe::Processor {
public:
  PileUpRWProcessor(const char *name);
  bool processEvent(xAOD::TEvent& event);
private:
  std::auto_ptr<CP::PileupReweightingTool> m_PileupTool_CENTRAL;
  std::auto_ptr<CP::PileupReweightingTool> m_PileupTool_UP;
  std::auto_ptr<CP::PileupReweightingTool> m_PileupTool_DOWN;
  int m_forcedRunNumber;

public:
  ClassDef(PileUpRWProcessor,0);
};

#endif // ZeroLeptonRun2_PileUpRWProcessor_H_
