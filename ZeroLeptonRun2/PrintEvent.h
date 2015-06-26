#ifndef ZeroLeptonRun2_PrintEvent_H_
#define ZeroLeptonRun2_PrintEvent_H_

#include "cafe/Processor.h"
#include <vector>
#include <string>

//--------------------------------------------------------------------------
// A simple processor to print event properties
//--------------------------------------------------------------------------
class PrintEvent : public cafe::Processor {
public:
  PrintEvent(const char *name);
  bool processEvent(xAOD::TEvent& event);

private:
  std::vector<std::string> m_jetKeys;
  std::vector<std::string> m_elKeys;
  std::vector<std::string> m_muKeys;
  std::vector<std::string> m_METKeys;

  // look for character '<' in the string and if the first two letters are "TS"
  // then change tag to the string following the '<' character and return true
  // return false otherwise
  bool inTStore(std::string& tag);

public:
  ClassDef(PrintEvent,0);
};

#endif // ZeroLeptonRun2_PrintEvent_H_
