#include "ZeroLeptonRun2/ZeroLeptonUtils.h"

#include "xAODRootAccess/TEvent.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/Vertex.h"

#include <stdexcept>
#include <cmath>

ZeroLeptonUtils::ZeroLeptonUtils(bool IsData, ZeroLeptonDerivationTag tag, 
				 const std::string& metKey): 
  m_IsData(IsData), m_derivationTag(tag), m_MET_Track_key("MET_Track"), 
  m_MET_key(metKey)
{
  if ( m_derivationTag == p1872 ) {
    m_MET_Track_key = "MET_TrackFix";
  }
}

bool ZeroLeptonUtils::NegCellCleaning(xAOD::TEvent& event, const TVector2& missingET) const
{
  bool isNegCell=false; 

  double MET_phi = TMath::ATan2(missingET.Y(),missingET.X()); 

  const xAOD::MissingETContainer* metcontainer = 0;
  if ( ! event.retrieve( metcontainer, m_MET_key ).isSuccess() ){
    throw std::runtime_error("Could not retrieve MissingETContainer with key "+m_MET_key);
  }
  xAOD::MissingETContainer::const_iterator met = metcontainer->find("SoftClus");
  if (met == metcontainer->end()) {
    throw std::runtime_error("MissingETContainer with key "+m_MET_key+" has no SoftClus component");
  }

  if ( ( (*met)->met()/missingET.Mod()) * std::cos((*met)->phi()-MET_phi) >= 0.5) isNegCell=true;
  
  return isNegCell ;

}

void  ZeroLeptonUtils::trackMET(xAOD::TEvent& event, double& met, double& phi) const
{
  const xAOD::MissingETContainer* metcontainer = 0;
  if ( !event.retrieve( metcontainer, m_MET_Track_key ).isSuccess() ){
    throw std::runtime_error("Could not retrieve MissingETContainer with key "+m_MET_Track_key);
  }
  // FIXME should be "Track" or "PVTrack" ?
  xAOD::MissingETContainer::const_iterator trackmet = metcontainer->find("Track");
  if ( trackmet ==  metcontainer->end() ) {
    throw std::runtime_error("Could not extract MET with name Track from container");
  }
  met = (*trackmet)->met();
  phi = (*trackmet)->phi();
}


const xAOD::Vertex* ZeroLeptonUtils::GetPrimVtx(xAOD::TEvent& event)
{
  const xAOD::VertexContainer* vertices(0);
  if( event.retrieve( vertices, "PrimaryVertices" ).isSuccess() ) {
    for( const auto& vx : *vertices ) {
      if(vx->vertexType() == xAOD::VxType::PriVtx) return vx;
    }
  }
  return NULL;
}

xAOD::JetInput::Type  ZeroLeptonUtils::JetTypeFromString(const std::string& algname)
{
  if ( algname.find(xAOD::JetInput::typeName(xAOD::JetInput::LCTopo)) != std::string::npos ) return xAOD::JetInput::LCTopo;
  if ( algname.find(xAOD::JetInput::typeName(xAOD::JetInput::EMTopo)) != std::string::npos ) return xAOD::JetInput::EMTopo;

  return xAOD::JetInput::Uncategorized;
}

