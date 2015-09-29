#include "ZeroLeptonRun2/PhysObjProxies.h"

#include "xAODJet/Jet.h"
#include "xAODEgamma/Electron.h"
#include "xAODEgamma/Photon.h"
#include "xAODMuon/Muon.h"
#include "xAODTau/TauJet.h"
#include "xAODTruth/TruthParticle.h"

JetProxy::JetProxy():
  TLorentzVector(),
  m_isBaseline(false),
  m_isSignal(false),
  m_isBad(true),
  m_passOR(true),
  m_isBJet(false),
  m_jet(0)
{
}

JetProxy::JetProxy(const TLorentzVector& in, bool isBaseline, bool isSignal, bool isBad, bool passOR, bool isBJet):
  TLorentzVector(in),
  m_isBaseline(isBaseline),
  m_isSignal(isSignal),
  m_isBad(isBad),
  m_passOR(passOR),
  m_isBJet(isBJet),
  m_jet(0)
{
}

JetProxy::JetProxy(const xAOD::Jet* jet):
  TLorentzVector(jet->p4())
{
  m_isBaseline = jet->auxdecor<char>("baseline")==1;
  m_isSignal   = jet->auxdecor<char>("signal")==1;
  m_isBad      = jet->auxdecor<char>("bad")==1;
  m_passOR     = jet->auxdecor<char>("passOR")==1;
  m_isBJet     = jet->auxdecor<char>("bjet")==1;
  m_jet        = jet;
}


ClassImp(JetProxy);



ElectronProxy::ElectronProxy():
  TLorentzVector(),
  m_isBaseline(false),
  m_isSignal(false),
  m_passOR(true),
  m_sf(0.f),
  m_el(0)
{
}

ElectronProxy::ElectronProxy(const xAOD::Electron* el):
  TLorentzVector(el->p4())
{
  m_isBaseline = el->auxdecor<char>("baseline")==1;
  m_isSignal   = el->auxdecor<char>("signal")==1;
  m_passOR     = el->auxdecor<char>("passOR")==1;
  m_sf         = el->auxdecor<float>("sf");
  m_el         = el;
}
ClassImp(ElectronProxy);

PhotonProxy::PhotonProxy():
  TLorentzVector(),
  m_isBaseline(false),
  m_isSignal(false),
  m_passOR(true),
  m_sf(0.f),
  m_ph(0)
{
}

PhotonProxy::PhotonProxy(const TLorentzVector& input):
  TLorentzVector(input),
  m_isBaseline(true),
  m_isSignal(true),
  m_passOR(true),
  m_sf(0.f),
  m_ph(0)
{
}

PhotonProxy::PhotonProxy(const xAOD::Photon* ph):
  TLorentzVector(ph->p4())
{
  m_isBaseline = ph->auxdecor<char>("baseline")==1;
  m_isSignal   = ph->auxdecor<char>("signal")==1;
  m_passOR     = ph->auxdecor<char>("passOR")==1;
  m_sf         = ph->auxdecor<float>("sf");
  m_ph         = ph;
}
ClassImp(PhotonProxy);


MuonProxy::MuonProxy():
  TLorentzVector(),
  m_isBaseline(false),
  m_isSignal(false),
  m_passOR(true),
  m_isCosmic(true),
  m_isBad(true),
  m_sf(0.f),
  m_muon(0)
{
}

MuonProxy::MuonProxy(const xAOD::Muon* muon):
  TLorentzVector(muon->p4())
{
  m_isBaseline = muon->auxdecor<char>("baseline")==1;
  m_isSignal   = muon->auxdecor<char>("signal")==1;
  m_passOR     = muon->auxdecor<char>("passOR")==1;
  m_isCosmic   = muon->auxdecor<char>("cosmic")==1;
  m_isBad      = muon->auxdecor<char>("bad")==1;
  m_sf         = muon->auxdecor<float>("sf");
  m_muon       = muon;
}
ClassImp(MuonProxy);

ElectronTruthProxy::ElectronTruthProxy():
  TLorentzVector(),
  m_eltruth(0)
{
}

ElectronTruthProxy::ElectronTruthProxy(const xAOD::TruthParticle* eltruth):
  TLorentzVector(eltruth->p4())
{
  m_eltruth         = eltruth;
}
ClassImp(ElectronTruthProxy);


MuonTruthProxy::MuonTruthProxy():
  TLorentzVector(),
  m_muontruth(0)
{
}

MuonTruthProxy::MuonTruthProxy(const xAOD::TruthParticle* muontruth):
  TLorentzVector(muontruth->p4())
{
  m_muontruth       = muontruth;
}
ClassImp(MuonTruthProxy);


TauProxy::TauProxy():
  TLorentzVector(),
  m_isBaseline(false),
  m_isSignal(false),
  m_sf(0.f),
  m_sfStatUp(0.f),
  m_sfStatDown(0.f),
  m_sfSystUp(0.f),
  m_sfSystDown(0.f),
  m_tau(0) 
{
}

TauProxy::TauProxy(const xAOD::TauJet* tau):
  TLorentzVector(tau->p4())
{
  m_isBaseline = tau->auxdecor<char>("baseline")==1;
  m_isSignal   = tau->auxdecor<char>("signal")==1;
  m_sf         = tau->auxdecor<float>("SFJetID");
  m_sfStatUp   = tau->auxdecor<float>("SFJetIDStatUp");
  m_sfStatDown = tau->auxdecor<float>("SFJetIDStatDown");
  m_sfSystUp   = tau->auxdecor<float>("SFJetIDSystUp");
  m_sfSystDown = tau->auxdecor<float>("SFJetIDSystDown");

  m_tau       = tau;
}
ClassImp(TauProxy);

