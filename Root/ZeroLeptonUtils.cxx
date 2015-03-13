#include "ZeroLeptonRun2/ZeroLeptonUtils.h"

#include "xAODRootAccess/TEvent.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODTracking/Vertex.h"

#include <stdexcept>
#include <cmath>

bool ZeroLeptonUtils::NegCellCleaning(xAOD::TEvent& event, const TVector2& missingET) const
{
  bool isNegCell=false; 

  double MET_phi = TMath::ATan2(missingET.Y(),missingET.X()); 

  const xAOD::MissingETContainer* metcontainer = 0;
  if ( ! event.retrieve( metcontainer, "MET_RefFinal" ).isSuccess() ){
    throw std::runtime_error("Could not retrieve MissingETContainer with key MET_RefFinal");
  }
  xAOD::MissingETContainer::const_iterator met = metcontainer->find("SoftClus");
  if (met == metcontainer->end()) {
    throw std::runtime_error("MissingETContainer with key MET_RefFinal has no SoftClus component");
  }

  if ( ( (*met)->met()/missingET.Mod()) * std::cos((*met)->phi()-MET_phi) >= 0.5) isNegCell=true;
  
  return isNegCell ;

}

void  ZeroLeptonUtils::trackMET(xAOD::TEvent& event, double& met, double& phi) const
{
  const xAOD::MissingETContainer* metcontainer = 0;
  if ( !event.retrieve( metcontainer, "MET_Track" ).isSuccess() ){
    throw std::runtime_error("Coulnot retrieve MissingETContainer with key MET_Track");
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
