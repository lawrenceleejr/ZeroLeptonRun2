#ifndef ZeroLeptonRun2_GRLProcessor_H_
#define ZeroLeptonRun2_GRLProcessor_H_

#include "cafe/Processor.h"

class GoodRunsListSelectionTool;


//------------------------------------------------
// GRLProcessor:
//   Select event based on Good Run lists
//------------------------------------------------

class GRLProcessor : public cafe::Processor {
public:
  GRLProcessor(const char *name);
  ~GRLProcessor();
  bool processEvent(xAOD::TEvent& event);
private:
  GoodRunsListSelectionTool* m_GRLtool;
  bool m_passAll;

public:
  ClassDef(GRLProcessor,0);
};

#endif // ZeroLeptonRun2_GRLProcessor_H_
