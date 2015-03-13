#ifndef ZeroLeptonRun2_MCInfoProcessor_H_
#define ZeroLeptonRun2_MCInfoProcessor_H_
#include "cafe/Processor.h"

namespace SUSY{
  class CrossSectionDB;
}
#include <vector>

class MCInfoProcessor : public cafe::Processor {
public:
  MCInfoProcessor(const char *name);
  bool processEvent(xAOD::TEvent& event);

private:
  void normWeights(xAOD::TEvent& event,std::vector<float>& normWeights,uint32_t mc_channel_number,unsigned int finalstate);
  unsigned int hardProcess(xAOD::TEvent& event) const;

  bool m_isSignal;
  SUSY::CrossSectionDB* m_mcDB;

public:
  ClassDef(MCInfoProcessor,0);
};

#endif // ZeroLeptonRun2_MCInfoProcessor_H_
