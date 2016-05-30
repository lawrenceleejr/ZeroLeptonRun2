
#include "ZeroLeptonRun2/TileTriggerCheck.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODJet/JetContainer.h"
#include "TH1F.h"

#include <cmath>

TileTriggerCheck::TileTriggerCheck(const char *name)
  : cafe::Processor(name), m_TileEnergy(0)
{
}

void TileTriggerCheck::begin()
{
  m_TileEnergy = new TH1F("JetTileEnergy","Jet energy in Tile",120,-200.,1000.);
  m_TileEnergy->SetDirectory(getDirectory());
}

bool TileTriggerCheck::processEvent(xAOD::TEvent& event)
{
  const xAOD::JetContainer* jets(0);
  if ( ! event.retrieve(jets,"AntiKt4EMTopoJets").isSuccess() ) throw std::runtime_error("TileTriggerCheck: Could not retrieve AntiKt4EMTopoJets");
  for ( const auto& jet : *jets ){
    float emf = 0.;
    jet->getAttribute(xAOD::JetAttribute::EMFrac,emf);
    float hecf = 0.;
    jet->getAttribute(xAOD::JetAttribute::HECFrac,hecf);
    //not in SUSY1 dervations 
    float EM_E = jet->jetP4(xAOD::JetEMScaleMomentum).E();
    //float EM_E = jet->e();
    
    //std::cout << "jet pt " << jet->pt() << " " << EM_E << " " << emf << " " << hecf << std::endl;

    m_TileEnergy->Fill(std::min(EM_E*(1.-emf-hecf)*0.001,999.99));
  }

  return true;
}

ClassImp(TileTriggerCheck);

