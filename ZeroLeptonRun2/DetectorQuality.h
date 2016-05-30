#ifndef ZeroLeptonRun2_DetectorQuality_H_
#define ZeroLeptonRun2_DetectorQuality_H_

#include "cafe/Processor.h"
#include "ZeroLeptonRun2/ZeroLeptonRunPeriod.h"
#include "ZeroLeptonRun2/ZeroLeptonNTVars.h"
class TTree;

//------------------------------------------------
// DetectorQuality:
//   select detector and processing quality flags
//   and write decision to keep event or not in
//   transient store
// -----------------------------------------------

class DetectorQuality : public cafe::Processor {
public:
    DetectorQuality(const char *name);
    bool processEvent(xAOD::TEvent& event);
    void begin();

private:

    ZeroLeptonRunPeriod m_period;
    bool m_doNtuple;
    DQVars m_dqv;
    TTree* m_tree;

public:
    ClassDef(DetectorQuality,1);
};

#endif // ZeroLeptonRun2_DetectorQuality_H_
