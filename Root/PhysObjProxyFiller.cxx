
#include "ZeroLeptonRun2/PhysObjProxyFiller.h"
#include "ZeroLeptonRun2/PhysObjProxies.h"
#include "ZeroLeptonRun2/PtOrder.h"

#include "xAODRootAccess/TActiveStore.h"
#include "xAODRootAccess/TStore.h"

#include "xAODMuon/MuonContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODTau/TauJetContainer.h"
#include "xAODTau/TauJet.h"


#include "xAODCore/ShallowCopy.h"
#include "xAODJetReclustering/JetReclusteringTool.h"
#include "xAODJet/JetAuxContainer.h"
#include "xAODCore/AuxContainerBase.h"

#include <iostream>


PhysObjProxyFiller::PhysObjProxyFiller(float jetPtCut, float elPtCut, float muonPtCut, float phPtCut, const std::string suffix, bool doRecl, const std::string suffixRecl, const std::string suffixSyst, const std::string suffixSmear):
  m_jetPtCut(jetPtCut), m_elPtCut(elPtCut), m_muonPtCut(muonPtCut), m_phPtCut(phPtCut),  m_doRecl(doRecl), m_suffix(suffix),m_suffixRecl(suffixRecl),m_suffixSyst(suffixSyst),m_suffixSmear(suffixSmear)
{
}

void PhysObjProxyFiller::setJRT()
{
  std::map<std::string, JetReclusteringTool*>::const_iterator pos = m_jrtMap.find(m_suffixSyst);
  if ( pos == m_jrtMap.end() ) {
    m_jetReclusteringTool = new JetReclusteringTool("ZLJetReclusteringTool"+m_suffixRecl+m_suffixSyst);
    m_jetReclusteringTool->setProperty("InputJetContainer",  "SUSYJetsNEW"+m_suffixRecl+m_suffixSyst);
    m_jetReclusteringTool->setProperty("OutputJetContainer", "SUSYJetsRecl"+m_suffixRecl+m_suffixSyst);
    m_jetReclusteringTool->setProperty("InputJetPtMin",      25.0);
    m_jetReclusteringTool->setProperty("RCJetPtMin",         50.0);
    m_jetReclusteringTool->setProperty("RCJetPtFrac",        0.05);
    m_jetReclusteringTool->initialize();
    m_jrtMap[m_suffixSyst] = m_jetReclusteringTool;
  }
  else {
    m_jetReclusteringTool =  m_jrtMap[m_suffixSyst];
  }
}

void PhysObjProxyFiller::FillJetProxies(std::vector<JetProxy>& good_jets,
					std::vector<JetProxy>& bad_jets,
					std::vector<JetProxy>& b_jets) const
{
  good_jets.clear();
  bad_jets.clear();
  xAOD::TStore* store = xAOD::TActiveStore::store();
  const xAOD::JetContainer* jets = 0;
  if ( !store->retrieve(jets, "SUSYJets"+m_suffix+m_suffixSmear).isSuccess() ) {
    throw std::runtime_error("Could not retrieve JetContainer with key SUSYJets"+m_suffix+m_suffixSmear);
  }
  for ( xAOD::JetContainer::const_iterator it = jets->begin();
	it != jets->end(); ++it ){
    //std::cout << " SUSY jet " <<  (*it)->pt() << " " <<  (*it)->eta() << " " << (*it)->phi()  << std::endl;
    if ( (*it)->pt() <= m_jetPtCut ) continue;
    if ( (*it)->auxdecor<char>("passOR") == 0) continue;
    if ( (*it)->auxdecor<char>("bad") == 0  ) {
      if ( std::abs((*it)->eta()) < 2.8 ) {
	good_jets.push_back(JetProxy(*it));
        if ( (*it)->auxdecor<char>("bjet") == 1
             && std::abs((*it)->eta()) < 2.5
             && (*it)->pt() >= 50000.) b_jets.push_back(JetProxy(*it));
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

//----------------------Begin Fat Jets---------------------------------//
void PhysObjProxyFiller::FillFatJetProxies(std::vector<JetProxy>& good_fat_jets,
					   std::vector<JetProxy>& bad_fat_jets,
					   std::vector<float>& vD2_fat,
					   std::vector<bool>& visWmedium_fat) const
{
  good_fat_jets.clear();
  bad_fat_jets.clear();
  vD2_fat.clear();
  xAOD::TStore* store = xAOD::TActiveStore::store();
  const xAOD::JetContainer* FatJets = 0;
  if ( !store->retrieve(FatJets, "SUSYFatJets"+m_suffix).isSuccess() ) {
    throw std::runtime_error("Could not retrieve JetContainer with key SUSYFatJets"+m_suffix);
  }
  
  std::vector<JetProxy> good_fat_jets_temp;
  std::vector<float> vD2_fat_temp;
  
  for ( xAOD::JetContainer::const_iterator it = FatJets->begin();
	it != FatJets->end(); ++it ){
    //std::cout << " SUSY jet " <<  (*it)->pt() << " " <<  (*it)->eta() << " " << (*it)->phi()  << std::endl;
    if ( (*it)->pt() <= m_jetPtCut ) continue;
    if ( (*it)->auxdecor<char>("passOR") == 0) continue;
    if ( (*it)->auxdecor<char>("bad") == 0  ) {
      if ( std::abs((*it)->eta()) < 2.8 ) {
	good_fat_jets.push_back(JetProxy(*it));
        good_fat_jets_temp.push_back(JetProxy(*it));

	float ECF1=-888;
	(*it)->getAttribute("ECF1",ECF1);
	
	float ECF2=-888;
	(*it)->getAttribute("ECF2",ECF2);
	
	float ECF3=-888;
	(*it)->getAttribute("ECF3",ECF3);	

	float D2_fat=-888;
	if(fabs(ECF2) > 0)
	  D2_fat = ECF3 * pow(ECF1,3.0)/pow(ECF2,3.0);

	vD2_fat_temp.push_back(D2_fat);

	if((*it)->auxdecor<int>("Wtag") == 1)
	  visWmedium_fat.push_back(1);
	else
	  visWmedium_fat.push_back(0);	
      }
    }
    
    else {
      // pT cut no longer applied in SUSYTools for bad jets
      if ( (*it)->pt() > 20000.) bad_fat_jets.push_back(JetProxy(*it));
    }
  }    
  // sort collections (calibration might have changed the order)
  std::sort(good_fat_jets.begin(),good_fat_jets.end(),PtOrder<JetProxy>);
  std::sort(bad_fat_jets.begin(),bad_fat_jets.end(),PtOrder<JetProxy>);
  
  // sort D2 according to fat jets 
  for(size_t i=0;i<good_fat_jets.size();i++){
    for(size_t j=0;j<good_fat_jets_temp.size();j++){
      if(abs(good_fat_jets[i].Pt()-good_fat_jets_temp[j].Pt())<0.00001){
  	vD2_fat.push_back(vD2_fat_temp[j]);
      }
    }
  }
}
//----------------------End Fat Jets---------------------------------//


void PhysObjProxyFiller::FillJetReclProxies(std::vector<JetProxy>& good_jets_recl,
                                            std::vector<float>& vD2)
{
  
  good_jets_recl.clear();
  vD2.clear();
  xAOD::TStore* store = xAOD::TActiveStore::store();
  const xAOD::JetContainer* jets = 0;
  if ( !store->retrieve(jets, "SUSYJets"+m_suffix).isSuccess() ) {
    throw std::runtime_error("Could not retrieve JetContainer with key SUSYJets"+m_suffix);
  }

  // NEEDED FOR RECLUSTERING
  xAOD::JetContainer* jetsNEW = new xAOD::JetContainer();
  xAOD::AuxContainerBase* jetsNEWAux = new xAOD::AuxContainerBase();
  jetsNEW->setStore(jetsNEWAux);
  for ( xAOD::JetContainer::const_iterator it = jets->begin();
        it != jets->end(); ++it ){
    if ( (*it)->pt() <= m_jetPtCut ) continue;
    if ( (*it)->auxdecor<char>("passOR") == 0) continue;
    if ( (*it)->auxdecor<char>("bad") == 0  ) {
      if ( std::abs((*it)->eta()) < 2.8 ) {
	xAOD::Jet* myjet = new xAOD::Jet();
        myjet->makePrivateStore(*it);
        jetsNEW->push_back(myjet) ;
      } 
    }
  }
  if ( ! store->record(jetsNEW,"SUSYJetsNEW"+m_suffixRecl+m_suffixSyst).isSuccess() ) throw std::runtime_error("Could not register SUSYJetsNEW"+m_suffixRecl+m_suffixSyst) ;
  if ( ! store->record(jetsNEWAux,"SUSYJetsNEWAux"+m_suffixRecl+m_suffixSyst).isSuccess() ) throw std::runtime_error("Could not register SUSYJetsNEWAux"+m_suffixRecl+m_suffixSyst) ;

  const xAOD::JetContainer* jetsrecl = 0;
  if(m_doRecl){
    setJRT();
    m_jetReclusteringTool->execute();
    if ( !store->retrieve(jetsrecl, "SUSYJetsRecl"+m_suffixRecl+m_suffixSyst).isSuccess() ) {
      throw std::runtime_error("Could not retrieve JetContainer with key SUSYJetsRecl"+m_suffixRecl+m_suffixSyst);
    }
    
    // Boson Tagging 
    //static JetSubStructureUtils::BosonTag WTagger("medium", "smooth", "$ROOTCOREBIN/data/JetSubStructureUtils/config_13TeV_20150528_Wtagging.dat", true, true);
    //static JetSubStructureUtils::BosonTag ZTagger("medium", "smooth", "$ROOTCOREBIN/data/JetSubStructureUtils/config_13TeV_20150528_Ztagging.dat", true, true);
    
    // need temporary objects for re-ordering
    int nrecl = 0;
    std::vector<JetProxy> good_jets_recl_temp;
    std::vector<float> vD2_temp;

    for ( xAOD::JetContainer::const_iterator itrecl = jetsrecl->begin();
          itrecl != jetsrecl->end(); ++itrecl ){

      good_jets_recl.push_back(JetProxy(*itrecl));
      good_jets_recl_temp.push_back(JetProxy(*itrecl));
      
      // components of the reclustered jet
      const xAOD::Jet* subjet(nullptr);
      int nsubjet = 0 ;
      std::vector<TLorentzVector> vtlvsubjet;

      for(auto constit: (*itrecl)->getConstituents()){
        subjet = static_cast<const xAOD::Jet*>(constit->rawConstituent());
        TLorentzVector tlvsubjet;
        tlvsubjet.SetPtEtaPhiE(constit->pt(),constit->eta(),constit->phi(),constit->E());
        vtlvsubjet.push_back(tlvsubjet);
        nsubjet++;
      }
      
      float ecf1=0;
      float ecf2=0;
      float ecf3=0;
      float D2=0;
      
      // https://cds.cern.ch/record/2014762/files/ATL-COM-PHYS-2015-378.pdf, L.184, Equation 7
      // i = subjets of JET
      // ECF1 = Sum_(i)    [ pT_i ]
      // ECF2 = Sum_(i<j)  [ pT_i*pT_j*dRij ]
      // ECF3 = SUM_(i<j<k)[ pT_i*pT_j*pT_k*dRij*dRik*dRjk ]
      
      for(size_t i=0;i<vtlvsubjet.size();i++){
        ecf1 += (vtlvsubjet.at(i)).Pt();
        for(size_t j=i+1;j<vtlvsubjet.size();j++){
          if(i<vtlvsubjet.size()-1){
            float dRij = (vtlvsubjet.at(i)).DeltaR(vtlvsubjet.at(j));
            ecf2 += (vtlvsubjet.at(i)).Pt()*(vtlvsubjet.at(j)).Pt()*dRij;
          }
          for(size_t k=j+1;k<vtlvsubjet.size();k++){
            if(j<vtlvsubjet.size()-1){
              float dRij = (vtlvsubjet.at(i)).DeltaR(vtlvsubjet.at(j));
              float dRik = (vtlvsubjet.at(i)).DeltaR(vtlvsubjet.at(k));
              float dRjk = (vtlvsubjet.at(j)).DeltaR(vtlvsubjet.at(k));
              ecf3 += (vtlvsubjet.at(i)).Pt()*(vtlvsubjet.at(j)).Pt()*(vtlvsubjet.at(k)).Pt()*dRij*dRik*dRjk;
            }
          } // 3rd loop on subjets
	} // 2nd loop on subjets
      } // 1st loop on subjets
      
      // D2 calculation
      // ecf3_beta2 * pow(ecf1_beta2, 3.0) / pow(ecf2_beta2, 3.0)
      if(fabs(ecf2) > 1e-8)
        D2 = ecf3 * pow(ecf1,3.0)/pow(ecf2,3.0);
      else
        D2 = -999.0;
      
      vD2_temp.push_back(D2);
      nrecl++;
    } // loop on recl jets
    
    // sort jet reclustered collections
    std::sort(good_jets_recl.begin(),good_jets_recl.end(),PtOrder<JetProxy>);
    
    // sort D2 according to jet reclustered 
    for(size_t i=0;i<good_jets_recl.size();i++){
      for(size_t j=0;j<good_jets_recl_temp.size();j++){
        if(abs(good_jets_recl[i].Pt()-good_jets_recl_temp[j].Pt())<0.00001){
          vD2.push_back(vD2_temp[j]);
        }
      }
    }
  }// do Reclustering
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

void PhysObjProxyFiller::FillPhotonProxies(std::vector<PhotonProxy>& baseline_photons,
					   std::vector<PhotonProxy>& isolated_baseline_photons,
					   std::vector<PhotonProxy>& isolated_signal_photons)
{
  baseline_photons.clear();
  isolated_baseline_photons.clear();
  isolated_signal_photons.clear();
  xAOD::TStore* store = xAOD::TActiveStore::store();
  const xAOD::PhotonContainer* photons = 0;
  if ( !store->retrieve(photons, "SUSYPhotons"+m_suffix).isSuccess() ){
    throw std::runtime_error("Could not retrieve PhotonContainer with key SUSYPhotons"+m_suffix);
  }

  for ( xAOD::PhotonContainer::const_iterator it = photons->begin();
        it != photons->end(); ++it ){

    if ( (*it)->pt() < m_phPtCut ) continue;
    if ( std::abs((*it)->eta())> 2.37 ) continue; 
    if ( (*it)->auxdecor<char>("baseline") == 0 ) continue;
    baseline_photons.push_back(PhotonProxy(*it));
    if ( (*it)->auxdecor<char>("passOR") == 1) {
      isolated_baseline_photons.push_back(PhotonProxy(*it));
      if ( (*it)->auxdecor<char>("signal") == 1) {
	isolated_signal_photons.push_back(PhotonProxy(*it));

	// isolation information 
	// Photon isolation /cvmfs/atlas.cern.ch/repo/sw/ASG/AnalysisBase/2.3.14/ElectronIsolationSelection/Root/IsolationSelectionTool.cxx 

	//if(!isData){
	//  vtruthType.push_back((*it)->auxdata<int>("truthType"));
	//  vtruthOrigin.push_back((*it)->auxdata<int>("truthOrigin"));
	//}
	//else{
	//  vtruthType.push_back(-1000);
	//  vtruthOrigin.push_back(-1000);
	//}
      }
    }
  }
  std::sort(baseline_photons.begin(),baseline_photons.end(),PtOrder<PhotonProxy>);
  std::sort(isolated_baseline_photons.begin(),isolated_baseline_photons.end(),PtOrder<PhotonProxy>);
  std::sort(isolated_signal_photons.begin(),isolated_signal_photons.end(),PtOrder<PhotonProxy>);
}
