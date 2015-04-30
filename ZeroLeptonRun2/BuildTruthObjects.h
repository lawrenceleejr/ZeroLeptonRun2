#ifndef ZeroLeptonRun2_BuildTruthObjects_H_
#define ZeroLeptonRun2_BuildTruthObjects_H_

#include "cafe/Processor.h"


#include "xAODJet/JetContainer.h"
#include "xAODJet/JetAuxContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthParticle.h"
#include "xAODTau/TauJetContainer.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODMissingET/MissingETAuxContainer.h"


#include <string>

//---------------------------------------------------------------------
// BuildTruthObjects use the collections of physics objects in the event
// The  output containers are stored in the transient store with keys "TruthJets", 
// "TruthMuons", "TruthElectrons", "TruthPhotons", "TruthTaus"
//
//---------------------------------------------------------------------
class BuildTruthObjects : public cafe::Processor {
public:
  BuildTruthObjects(const char *name);
  ~BuildTruthObjects();
  bool processEvent(xAOD::TEvent& event);
  bool OverlapRemoval(const xAOD::TruthParticleContainer *electrons, const xAOD::TruthParticleContainer *muons, const xAOD::JetContainer *jets, bool doHarmonization, double dRejet, double dRjetmu, double dRjete, double dRemu, double dRee);
  bool OverlapRemoval(const xAOD::TruthParticleContainer *electrons, const xAOD::TruthParticleContainer *muons, const xAOD::JetContainer *jets, const xAOD::TruthParticleContainer *photons, const bool doHarmonization, const double dRejet, const double dRjetmu, const double dRjete, double dRemu, double dRee, double dRphjet, double dReph, double dRmuph);

private:
  bool m_IsData;

  std::string m_jetkey; 
  std::string m_suffix;
public:
  ClassDef(BuildTruthObjects,0);
};

#endif // ZeroLeptonRun2_BuildTruthObjects_H_
