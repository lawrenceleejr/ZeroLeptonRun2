#ifndef ZeroLeptonRun2_MCEventVetoHelper_H_
#define ZeroLeptonRun2_MCEventVetoHelper_H_

#include "xAODJet/JetContainer.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODTruth/TruthParticleContainer.h"

// -------------------------------------------------------------------
//  MCEventVetoHelper:
//    helper class for MCEventVeto to work around the limitations of
//    xAOD EDM vs CINT
// -------------------------------------------------------------------


class MCEventVetoHelper
{
 public:
  MCEventVetoHelper() {};

  static bool isHighPtDijet(const xAOD::JetContainer* jets);

  static bool isHighPtJetMET(uint32_t mc_channel_number, 
			     const xAOD::JetContainer* jets,
			     const xAOD::MissingETContainer* metc);

  // general function for truth veto on mc12 samples
  static bool mc12accept(unsigned int& veto,
			 uint32_t mc_channel_number, 
			 const xAOD::TruthParticleContainer* mcparticles, 
			 const xAOD::MissingETContainer* metc);

  // helper function for mc12accept
  static unsigned int vetoQEDFSR(uint32_t mc_channel_number, 
				 const xAOD::TruthParticleContainer* mcparticles);
  static bool mc12AlpgenZjets_accept(uint32_t mc_channel_number, 
				     const xAOD::TruthParticleContainer* mcparticles);
  static bool mc12AlpgenYjets_accept(uint32_t mc_channel_number, 
				     const xAOD::TruthParticleContainer* mcparticles);
  static void mc12AlpgenJimmyW_accept(unsigned int& veto,
				      uint32_t mc_channel_number, 
				      const xAOD::TruthParticleContainer* mcparticles, 
				      const xAOD::MissingETContainer* metc);
  static void mc12SherpaZnunu_accept(unsigned int& veto,
				     uint32_t mc_channel_number, 
				     const xAOD::TruthParticleContainer* mcparticles);
  static void mc12SherpaYjets_accept(unsigned int& veto,
				     uint32_t mc_channel_number, 
				     const xAOD::TruthParticleContainer* mcparticles);
  static bool trueBosonFromWorZplusJetsMCSample(TLorentzVector&trueBoson, 
						uint32_t mc_channel_number,
						const xAOD::TruthParticleContainer* mcparticles);
  static bool mc12SherpaWZjets_accept(unsigned int& veto,
				      uint32_t mc_channel_number, 
				      const xAOD::TruthParticleContainer* mcparticles);

  static void mc12HerwigVVjets_accept(unsigned int& veto,
				      uint32_t mc_channel_number, 
				      const xAOD::TruthParticleContainer* mcparticles);

 private:

};

#endif
