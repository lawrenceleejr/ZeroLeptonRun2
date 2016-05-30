#ifndef ZeroLeptonUtils_h
#define ZeroLeptonUtils_h


namespace xAOD {
  class TEvent;
}

#include "TVector2.h"
#include "xAODTracking/VertexFwd.h"
#include "xAODJet/JetContainerInfo.h"

//--------------------------------------------------------------------------
// Utility functions that do not depend on physics object proxy and may
// require access to the event store
//
//--------------------------------------------------------------------------
class ZeroLeptonUtils
{
 public:
  ZeroLeptonUtils(bool IsData, const std::string& metKey = "MET_ZL");

  void trackMET(xAOD::TEvent& event, double& met, double& phi) const;

  // copied from SUSYObjDef_xAOD since it is not static there
  static const xAOD::Vertex* GetPrimVtx(xAOD::TEvent& event);

  static  xAOD::JetInput::Type  JetTypeFromString(const std::string& algname);

 private:
  bool m_IsData;

  std::string m_MET_Track_key;
  std::string m_MET_key;
};
#endif
