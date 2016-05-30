#ifndef ZeroLeptonRun2_MCInfoProcessor_H_
#define ZeroLeptonRun2_MCInfoProcessor_H_
#include "cafe/Processor.h"

namespace SUSY{
  class CrossSectionDB;
}
#include <vector>
#include <unordered_map>
#include <string>

class MCInfoProcessor : public cafe::Processor {
public:
  MCInfoProcessor(const char *name);
  bool processEvent(xAOD::TEvent& event);
  void finish();

private:
  void normWeights(xAOD::TEvent& event,std::vector<float>& normWeights,uint32_t mc_channel_number,unsigned int finalstate);
  unsigned int hardProcess(xAOD::TEvent& event) const;

  std::string m_truthPKey;
  bool m_isSignal;
  SUSY::CrossSectionDB* m_mcDB;

  bool m_buildDB;
  uint32_t m_runNumber;
  uint32_t m_channelNumber;
  unsigned int m_eventCounter;
  std::unordered_map<unsigned int,unsigned int> m_hardproc;

public:
  ClassDef(MCInfoProcessor,0);
};

#endif // ZeroLeptonRun2_MCInfoProcessor_H_
