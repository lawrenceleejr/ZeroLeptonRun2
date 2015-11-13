
#include "ZeroLeptonRun2/PhysObjProxyFillerTruth.h"
#include "ZeroLeptonRun2/PhysObjProxies.h"
#include "ZeroLeptonRun2/PtOrder.h"
//#include "xAODTruth/xAODTruthHelpers.h"
#include "MCTruthClassifier/MCTruthClassifierDefs.h"

#include "xAODRootAccess/TActiveStore.h"
#include "xAODRootAccess/TStore.h"

#include "xAODMuon/MuonContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODBTagging/BTagging.h"
//#include "xAODEgamma/PhotonContainer.h"
//#include "xAODTau/TauJetContainer.h"
#include "AthContainers/AuxElement.h"

#include <iostream>

static const SG::AuxElement::ConstAccessor<unsigned int> acc_truthType("classifierParticleType");

PhysObjProxyFillerTruth::PhysObjProxyFillerTruth(float jetPtCut, float elPtCut, float muonPtCut, float phPtCut, const std::string suffix):
  m_jetPtCut(jetPtCut), m_elPtCut(elPtCut), m_muonPtCut(muonPtCut), m_phPtCut(phPtCut), m_suffix(suffix)
{
}


void PhysObjProxyFillerTruth::FillJetProxies(std::vector<JetProxy>& good_jets,
					std::vector<JetProxy>& b_jets) const
{
  good_jets.clear();
  xAOD::TStore* store = xAOD::TActiveStore::store();
  const xAOD::JetContainer* jets = 0;
  if ( !store->retrieve(jets, "myTruthJets"+m_suffix).isSuccess() ) {
    throw std::runtime_error("Could not retrieve JetContainer with key myTruthJets"+m_suffix);
  }
  for ( xAOD::JetContainer::const_iterator it = jets->begin();
	it != jets->end(); ++it ){
    if ( (*it)->pt() <= m_jetPtCut ) continue;
    if ( (*it)->auxdecor<char>("passOR") == 0) continue;
    if ( std::abs((*it)->eta()) < 2.8 ) {
      // std::cout << " Truth jet " <<  (*it)->pt() << " " <<  (*it)->eta() << " " << (*it)->phi()  << " Cut:  " << m_jetPtCut << std::endl;
      good_jets.push_back(JetProxy(*it));
      // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/Run2JetMoments
      int tagInfo = (*it)->getAttribute<int>("ConeTruthLabelID");
      if ( tagInfo==5 ) b_jets.push_back(JetProxy(*it));
    }
    //}
  }

  // sort collections (calibration might have changed the order)
  std::sort(good_jets.begin(),good_jets.end(),PtOrder<JetProxy>);
}



void PhysObjProxyFillerTruth::FillElectronProxies(std::vector<ElectronTruthProxy>& baseline_electrons,
					     std::vector<ElectronTruthProxy>& isolated_baseline_electrons,
					     std::vector<ElectronTruthProxy>& isolated_signal_electrons)
{
  baseline_electrons.clear();
  isolated_baseline_electrons.clear();
  isolated_signal_electrons.clear();
  xAOD::TStore* store = xAOD::TActiveStore::store();
  const xAOD::TruthParticleContainer* electrons = 0;
  if ( !store->retrieve(electrons, "myTruthElectrons"+m_suffix).isSuccess() ){
    throw std::runtime_error("Could not retrieve ElectronContainer with key myTruthElectrons");
  }

  for ( xAOD::TruthParticleContainer::const_iterator it = electrons->begin();
	it != electrons->end(); ++it ){
    // std::cout << " Truth el " <<  (*it)->pt() << " " << (*it)->eta() << " Cut : " << m_elPtCut << std::endl;
    if ( (*it)->pt() < m_elPtCut ) continue;
    if ( (*it)->auxdecor<char>("passOR") == 0) continue;
    if ( std::abs((*it)->eta()) < 2.47 ) {
     // https://svnweb.cern.ch/trac/atlasoff/browser/PhysicsAnalysis/MCTruthClassifier/trunk/MCTruthClassifier/MCTruthClassifierDefs.h
      if(acc_truthType(**it)==MCTruthPartClassifier::IsoElectron){
	baseline_electrons.push_back(ElectronTruthProxy(*it));
	isolated_baseline_electrons.push_back(ElectronTruthProxy(*it));
	isolated_signal_electrons.push_back(ElectronTruthProxy(*it));
      }
    }
  }

  // sort collections (calibration might have changed the order)
  std::sort(baseline_electrons.begin(),baseline_electrons.end(),PtOrder<ElectronTruthProxy>);
  std::sort(isolated_baseline_electrons.begin(),isolated_baseline_electrons.end(),PtOrder<ElectronTruthProxy>);
  std::sort(isolated_signal_electrons.begin(),isolated_signal_electrons.end(),PtOrder<ElectronTruthProxy>);
}



void PhysObjProxyFillerTruth::FillMuonProxies(std::vector<MuonTruthProxy>& baseline_muons,
					 std::vector<MuonTruthProxy>& isolated_baseline_muons,
					 std::vector<MuonTruthProxy>& isolated_signal_muons)
{
  baseline_muons.clear();
  isolated_baseline_muons.clear();
  isolated_signal_muons.clear();
  xAOD::TStore* store = xAOD::TActiveStore::store();
  const xAOD::TruthParticleContainer* muons = 0;
  if ( !store->retrieve(muons, "myTruthMuons"+m_suffix).isSuccess() ){
    throw std::runtime_error("Could not retrieve MuonContainer with key myTruthMuons");
  }

  for ( xAOD::TruthParticleContainer::const_iterator it = muons->begin();
	it != muons->end(); ++it ){
    //std::cout << " Truth muon " <<  (*it)->pt() << " " << (*it)->eta() << " Cut : " << m_muonPtCut << std::endl;
    if ( (*it)->pt() < m_muonPtCut ) continue;
    if ( (*it)->auxdecor<char>("passOR") == 0) continue;
    if ( std::abs((*it)->eta()) < 2.4 ) {
      // https://svnweb.cern.ch/trac/atlasoff/browser/PhysicsAnalysis/MCTruthClassifier/trunk/MCTruthClassifier/MCTruthClassifierDefs.h
      if(acc_truthType(**it)==MCTruthPartClassifier::IsoMuon){
	baseline_muons.push_back(MuonTruthProxy(*it));
	isolated_baseline_muons.push_back(MuonTruthProxy(*it));
	isolated_signal_muons.push_back(MuonTruthProxy(*it));
      }
    }
  }

  // sort collections (calibration might have changed the order)
  std::sort(baseline_muons.begin(),baseline_muons.end(),PtOrder<MuonTruthProxy>);
  std::sort(isolated_baseline_muons.begin(),isolated_baseline_muons.end(),PtOrder<MuonTruthProxy>);
  std::sort(isolated_signal_muons.begin(),isolated_signal_muons.end(),PtOrder<MuonTruthProxy>);
}

void PhysObjProxyFillerTruth::FillPhotonProxies(std::vector<PhotonTruthProxy>& baseline_photons,
						std::vector<PhotonTruthProxy>& isolated_baseline_photons,
						std::vector<PhotonTruthProxy>& isolated_signal_photons)
{
  baseline_photons.clear();
  isolated_baseline_photons.clear();
  isolated_signal_photons.clear();
  xAOD::TStore* store = xAOD::TActiveStore::store();
  const xAOD::TruthParticleContainer* photons = 0;
  if ( !store->retrieve(photons, "myTruthPhotons"+m_suffix).isSuccess() ){
    throw std::runtime_error("Could not retrieve PhotonContainer with key myTruthPhotons");
  }
  for ( xAOD::TruthParticleContainer::const_iterator it = photons->begin();
        it != photons->end(); ++it ){
    if ( (*it)->pt() < m_phPtCut ) continue;
    if ( (*it)->auxdecor<char>("passOR") == 0) continue;
    if ( std::abs((*it)->eta()) < 2.37 ) {
      if(acc_truthType(**it)==MCTruthPartClassifier::IsoPhoton){
	// std::cout << "photon truth type          : " << acc_truthType(**it) << std::endl;
	// std::cout << " Truth photon " <<  (*it)->pt() << " " <<  (*it)->eta() << " " << (*it)->phi()  << " Cut:  " << m_phPtCut << std::endl;
        baseline_photons.push_back(PhotonTruthProxy(*it));
        isolated_baseline_photons.push_back(PhotonTruthProxy(*it));
        isolated_signal_photons.push_back(PhotonTruthProxy(*it));
      }
    }
  }

  std::cout << "Event has " << isolated_signal_photons.size() << " good photons." << std::endl;

  std::sort(baseline_photons.begin(),baseline_photons.end(),PtOrder<PhotonTruthProxy>);
  std::sort(isolated_baseline_photons.begin(),isolated_baseline_photons.end(),PtOrder<PhotonTruthProxy>);
  std::sort(isolated_signal_photons.begin(),isolated_signal_photons.end(),PtOrder<PhotonTruthProxy>);
}

//  LocalWords:  classifierParticleType
