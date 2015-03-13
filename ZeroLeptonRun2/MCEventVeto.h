#ifndef ZeroLeptonRun2_MCEventVeto_H_
#define ZeroLeptonRun2_MCEventVeto_H_

#include "cafe/Processor.h"
#include "ZeroLeptonRun2/ZeroLeptonRunPeriod.h"


// -------------------------------------------------------------------
//  MCEventVeto:
//     reject event based on MC sample and truth event properties,
//     e.g. to reject known buggy events or events that would cause
//     double counting
// -------------------------------------------------------------------


class MCEventVeto : public cafe::Processor {
public:
    MCEventVeto(const char *name);
    bool processEvent(xAOD::TEvent& event);

private:
    ZeroLeptonRunPeriod m_period;

public:
    ClassDef(MCEventVeto,0);
};

#endif // ZeroLeptonRun2_MCEventVeto_H_
