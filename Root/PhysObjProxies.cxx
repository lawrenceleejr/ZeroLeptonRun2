#include "ZeroLeptonRun2/PhysObjProxies.h"

#include "xAODJet/Jet.h"
#include "xAODEgamma/Electron.h"
#include "xAODMuon/Muon.h"

JetProxy::JetProxy():
  TLorentzVector(),
  m_isBaseline(false),
  m_isBad(true),
  m_passOR(true),
  m_isBJet(false),
  m_jet(0)
{
}

JetProxy::JetProxy(const TLorentzVector& in, bool isBaseline, bool isBad, bool passOR, bool isBJet):
  TLorentzVector(in),
  m_isBaseline(isBaseline),
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
  m_el(0)
{
}

ElectronProxy::ElectronProxy(const xAOD::Electron* el):
  TLorentzVector(el->p4())
{
  m_isBaseline = el->auxdecor<char>("baseline")==1;
  m_isSignal   = el->auxdecor<char>("signal")==1;
  m_passOR     = el->auxdecor<char>("passOR")==1;
  m_el         = el;
}
ClassImp(ElectronProxy);



MuonProxy::MuonProxy():
  TLorentzVector(),
  m_isBaseline(false),
  m_isSignal(false),
  m_passOR(true),
  m_isCosmic(true),
  m_muon(0)
{
}

MuonProxy::MuonProxy(const xAOD::Muon_v1* muon):
  TLorentzVector(muon->p4())
{
  m_isBaseline = muon->auxdecor<char>("baseline")==1;
  m_isSignal   = muon->auxdecor<char>("signal")==1;
  m_passOR     = muon->auxdecor<char>("passOR")==1;
  m_isCosmic   = muon->auxdecor<char>("cosmic")==1;
  m_muon       = muon;
}
ClassImp(MuonProxy);
