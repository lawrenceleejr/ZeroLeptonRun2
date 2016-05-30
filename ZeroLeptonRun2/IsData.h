#ifndef ZeroLeptonRun2_IsData_H_
#define ZeroLeptonRun2_IsData_H_

#include "cafe/Processor.h"

// -------------------------------------------------------------
// IsData: processEvent() return true if event from data
// -------------------------------------------------------------
class IsData : public cafe::Processor {
public:
  IsData(const char *name);
  bool processEvent(xAOD::TEvent& event);
private:
  bool m_enforce;
  bool m_expected;

public:
  ClassDef(IsData,0);
};

#endif // ZeroLeptonRun2_IsData_H_
