
#include "ZeroLeptonRun2/PhysObjProxyFiller.h"
#include "ZeroLeptonRun2/PhysObjProxies.h"
#include "ZeroLeptonRun2/PtOrder.h"

#include "xAODRootAccess/TActiveStore.h"
#include "xAODRootAccess/TStore.h"

#include "xAODMuon/MuonContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODEgamma/ElectronContainer.h"
//#include "xAODEgamma/PhotonContainer.h"
#include "xAODTau/TauJetContainer.h"
#include "xAODTau/TauJet.h"

#include <iostream>


PhysObjProxyFiller::PhysObjProxyFiller(float jetPtCut, float elPtCut, float muonPtCut, const std::string suffix):
  m_jetPtCut(jetPtCut), m_elPtCut(elPtCut), m_muonPtCut(muonPtCut), m_suffix(suffix)
{
}


void PhysObjProxyFiller::FillJetProxies(std::vector<JetProxy>& good_jets,
					std::vector<JetProxy>& bad_jets,
					std::vector<JetProxy>& b_jets) const
{
  good_jets.clear();
  bad_jets.clear();
  xAOD::TStore* store = xAOD::TActiveStore::store();
  const xAOD::JetContainer* jets = 0;
  if ( !store->retrieve(jets, "SUSYJets"+m_suffix).isSuccess() ) {
    throw std::runtime_error("Could not retrieve JetContainer with key SUSYJets"+m_suffix);
  }
  for ( xAOD::JetContainer::const_iterator it = jets->begin();
	it != jets->end(); ++it ){
    //std::cout << " SUSY jet " <<  (*it)->pt() << " " <<  (*it)->eta() << " " << (*it)->phi()  << std::endl;
    if ( (*it)->pt() <= m_jetPtCut ) continue;
    if ( (*it)->auxdecor<char>("passOR") == 0) continue;
    if ( (*it)->auxdecor<char>("bad") == 0  ) {
      if ( std::abs((*it)->eta()) < 2.8 ) {
	good_jets.push_back(JetProxy(*it));
	if ( (*it)->auxdecor<char>("bjet") == 1) b_jets.push_back(JetProxy(*it));
      }
    }
    else {
      // pT cut no longer applied in SUSYTools for bad jets
      if ( (*it)->pt() > 20000.) bad_jets.push_back(JetProxy(*it));
    }
  }

  // sort collections (calibration might have changed the order)
  std::sort(good_jets.begin(),good_jets.end(),PtOrder<JetProxy>);
  std::sort(bad_jets.begin(),bad_jets.end(),PtOrder<JetProxy>);
}



void PhysObjProxyFiller::FillElectronProxies(std::vector<ElectronProxy>& baseline_electrons,
					     std::vector<ElectronProxy>& isolated_baseline_electrons,
					     std::vector<ElectronProxy>& isolated_signal_electrons)
{
  baseline_electrons.clear();
  isolated_baseline_electrons.clear();
  isolated_signal_electrons.clear();
  xAOD::TStore* store = xAOD::TActiveStore::store();
  const xAOD::ElectronContainer* electrons = 0;
  if ( !store->retrieve(electrons, "SUSYElectrons"+m_suffix).isSuccess() ){
    throw std::runtime_error("Could not retrieve ElectronContainer with key SUSYElectrons"+m_suffix);
  }

  for ( xAOD::ElectronContainer::const_iterator it = electrons->begin();
	it != electrons->end(); ++it ){
    //std::cout << " SUSY el " <<  (*it) << " " << (*it)->pt() << " " << (int)((*it)->auxdecor<char>("passOR"))  << " " << (int)((*it)->auxdecor<char>("baseline"))  << " " << (int)((*it)->auxdecor<char>("signal")) << " " << (*it)->charge() << std::endl;
    if ( (*it)->pt() < m_elPtCut ) continue;
    if ( (*it)->auxdecor<char>("baseline") == 0 ) continue;
    baseline_electrons.push_back(ElectronProxy(*it));
    if ( (*it)->auxdecor<char>("passOR") == 1) {
      isolated_baseline_electrons.push_back(ElectronProxy(*it));
      if ( (*it)->auxdecor<char>("signal") == 1) {
	isolated_signal_electrons.push_back(ElectronProxy(*it));
      }
    }
  }

  // sort collections (calibration might have changed the order)
  std::sort(baseline_electrons.begin(),baseline_electrons.end(),PtOrder<ElectronProxy>);
  std::sort(isolated_baseline_electrons.begin(),isolated_baseline_electrons.end(),PtOrder<ElectronProxy>);
  std::sort(isolated_signal_electrons.begin(),isolated_signal_electrons.end(),PtOrder<ElectronProxy>);
}



void PhysObjProxyFiller::FillMuonProxies(std::vector<MuonProxy>& baseline_muons,
					 std::vector<MuonProxy>& isolated_baseline_muons,
					 std::vector<MuonProxy>& isolated_signal_muons)
{
  baseline_muons.clear();
  isolated_baseline_muons.clear();
  isolated_signal_muons.clear();
  xAOD::TStore* store = xAOD::TActiveStore::store();
  const xAOD::MuonContainer* muons = 0;
  if ( !store->retrieve(muons, "SUSYMuons"+m_suffix).isSuccess() ){
    throw std::runtime_error("Could not retrieve MuonContainer with key SUSYMuons"+m_suffix);
  }

  for ( xAOD::MuonContainer::const_iterator it = muons->begin();
	it != muons->end(); ++it ){
    //std::cout << " SUSY muon " <<  (*it) << " " << (*it)->pt() << " " << (int)((*it)->auxdecor<char>("passOR"))  << " " << (int)((*it)->auxdecor<char>("baseline"))  << " " << (int)((*it)->auxdecor<char>("signal")) << " " << (*it)->charge() << std::endl;
    if ( (*it)->pt() < m_muonPtCut ) continue;
    if ( (*it)->auxdecor<char>("baseline") == 0 ) continue;
    baseline_muons.push_back(MuonProxy(*it));
    if (  (*it)->auxdecor<char>("passOR") == 1) {
      isolated_baseline_muons.push_back(MuonProxy(*it));
      if ( (*it)->auxdecor<char>("signal") == 1) {
	isolated_signal_muons.push_back(MuonProxy(*it));
      }
    }
  }

  // sort collections (calibration might have changed the order)
  std::sort(baseline_muons.begin(),baseline_muons.end(),PtOrder<MuonProxy>);
  std::sort(isolated_baseline_muons.begin(),isolated_baseline_muons.end(),PtOrder<MuonProxy>);
  std::sort(isolated_signal_muons.begin(),isolated_signal_muons.end(),PtOrder<MuonProxy>);
}



void PhysObjProxyFiller::FillTauProxies(std::vector<TauProxy>& baseline_taus,
					std::vector<TauProxy>& signal_taus)
{
  baseline_taus.clear();
  signal_taus.clear();
  xAOD::TStore* store = xAOD::TActiveStore::store();
  const xAOD::TauJetContainer* taus = 0;
  if ( !store->retrieve(taus, "SUSYTaus").isSuccess() ){
    throw std::runtime_error("Could not retrieve TauJetContainer with key SUSYTaus");
  }

  for ( xAOD::TauJetContainer::const_iterator it = taus->begin();
	it != taus->end(); ++it ){
    if ( (*it)->auxdecor<char>("baseline") == 0 ) continue;
    baseline_taus.push_back(TauProxy(*it));
    if ( (*it)->auxdecor<char>("signal") == 1) {
      signal_taus.push_back(TauProxy(*it));
    }
  }

  // sort collections (calibration might have changed the order)
  std::sort(baseline_taus.begin(),baseline_taus.end(),PtOrder<TauProxy>);
  std::sort(signal_taus.begin(),signal_taus.end(),PtOrder<TauProxy>);
}
