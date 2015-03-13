#ifndef ZeroLeptonUtils_h
#define ZeroLeptonUtils_h


namespace xAOD {
  class TEvent;
}

#include "TVector2.h"
#include "xAODTracking/VertexFwd.h"

//--------------------------------------------------------------------------
// Utility functions that do not depend on physics object proxy and my
// require access to the event store
//
//--------------------------------------------------------------------------
class ZeroLeptonUtils
{
 public:
  ZeroLeptonUtils(bool IsData): m_IsData(IsData) {}

  bool NegCellCleaning(xAOD::TEvent& event, const TVector2& missingET) const;

  void trackMET(xAOD::TEvent& event, double& met, double& phi) const;

  // copied from SUSYObjDef_xAOD since it is not static there
  static const xAOD::Vertex* GetPrimVtx(xAOD::TEvent& event);

 private:
  bool m_IsData;
};
#endif
