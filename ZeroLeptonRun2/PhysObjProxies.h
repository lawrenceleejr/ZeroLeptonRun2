#ifndef PhysObjProxies_h
#define PhysObjProxies_h

#include "TLorentzVector.h"

//FIXME once forward declaration works
namespace xAOD {
  class Jet_v1;
  typedef Jet_v1 Jet;
}
namespace xAOD {
  class Muon_v1;
}
namespace xAOD {
  class TruthParticle_v1;
  typedef TruthParticle_v1 TruthParticle;
}
#include "xAODEgamma/ElectronFwd.h"
#include "xAODEgamma/PhotonFwd.h"


//----------------------------------------------------------------------
// These classes act as proxies for the xAOD physics objects in order
// to be able to add a lepton to jets in CR and present properties
// defined in SUSYOBjDef_xAOD as class members.
//
// The alternative would be to use xAOD::IParticle* and dynamic_cast<>
// with a code that is less readable and possibly slow.
//----------------------------------------------------------------------
class JetProxy: public TLorentzVector
{
 public:
  JetProxy();
  JetProxy(const TLorentzVector& in, bool isBaseline, bool isBad, bool passOR, bool isBJet);
  JetProxy(const xAOD::Jet* jet);

  inline bool isBaseline() const {return m_isBaseline;}
  inline bool isBad() const {return m_isBad;}
  inline bool isBJet() const {return m_isBJet;}
  inline bool passOVerlapRemoval() const {return m_passOR;}
  inline const xAOD::Jet* jet() const {return m_jet;}

 private:
  bool m_isBaseline;
  bool m_isBad;
  bool m_passOR;
  bool m_isBJet;
  const xAOD::Jet* m_jet;

public:
    ClassDef(JetProxy,0);
};


class ElectronProxy: public TLorentzVector
{
 public:
  ElectronProxy();
  ElectronProxy(const xAOD::Electron* el);

  inline bool isBaseline() const {return m_isBaseline;}
  inline bool isSignal() const {return m_isSignal;}
  inline bool passOVerlapRemoval() const {return m_passOR;}
  inline const xAOD::Electron* electron() const {return m_el;}

 private:
  bool m_isBaseline;
  bool m_isSignal;
  bool m_passOR;
  const xAOD::Electron* m_el;

public:
    ClassDef(ElectronProxy,0);
};

class PhotonProxy: public TLorentzVector
{
 public:
  PhotonProxy();
  PhotonProxy(const xAOD::Photon* ph);
  PhotonProxy(const TLorentzVector& input);

  inline bool isBaseline() const {return m_isBaseline;}
  inline bool isSignal() const {return m_isSignal;}
  inline bool passOVerlapRemoval() const {return m_passOR;}
  inline const xAOD::Photon* photon() const {return m_ph;}

 private:
  bool m_isBaseline;
  bool m_isSignal;
  bool m_passOR;
  const xAOD::Photon* m_ph;

public:
    ClassDef(PhotonProxy,1);
};

class MuonProxy: public TLorentzVector
{
 public:
  MuonProxy();
  MuonProxy(const xAOD::Muon_v1* muon);

  inline bool isBaseline() const {return m_isBaseline;}
  inline bool isSignal() const {return m_isSignal;}
  inline bool passOVerlapRemoval() const {return m_passOR;}
  inline bool isCosmic() const {return m_isCosmic;}
  inline const xAOD::Muon_v1* muon() const {return m_muon;}

 private:
  bool m_isBaseline;
  bool m_isSignal;
  bool m_passOR;
  bool m_isCosmic;
  const xAOD::Muon_v1* m_muon;

public:
    ClassDef(MuonProxy,0);
};

class ElectronTruthProxy: public TLorentzVector
{
 public:
  ElectronTruthProxy();
  ElectronTruthProxy(const xAOD::TruthParticle* eltruth);
  inline const xAOD::TruthParticle* eltruth() const {return m_eltruth;}

 private:
  const xAOD::TruthParticle* m_eltruth;

 public:
  ClassDef(ElectronTruthProxy,0);
};

class MuonTruthProxy: public TLorentzVector
{
 public:
  MuonTruthProxy();
  MuonTruthProxy(const xAOD::TruthParticle* muontruth);
  inline const xAOD::TruthParticle* muontruth() const {return m_muontruth;}

 private:
  const xAOD::TruthParticle* m_muontruth;

 public:
  ClassDef(MuonTruthProxy,0);
};


#endif
