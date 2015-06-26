#ifndef ZeroLeptonRun2_DetectorQuality_H_
#define ZeroLeptonRun2_DetectorQuality_H_

#include "cafe/Processor.h"
#include "ZeroLeptonRun2/ZeroLeptonRunPeriod.h"

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

private:
    ZeroLeptonRunPeriod m_period;

public:
    ClassDef(DetectorQuality,0);
};

#endif // ZeroLeptonRun2_DetectorQuality_H_
