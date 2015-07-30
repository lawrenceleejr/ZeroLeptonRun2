#ifndef ZeroLeptonRun2_ZeroLeptonCRWT_H_
#define ZeroLeptonRun2_ZeroLeptonCRWT_H_

#include "cafe/Processor.h"
#include "ZeroLeptonRun2/ZeroLeptonNTVars.h"
#include "ZeroLeptonRun2/PhysObjProxyFiller.h"
#include "ZeroLeptonRun2/PhysObjProxyFillerTruth.h"
#include "ZeroLeptonRun2/ZeroLeptonUtils.h"
#include "ZeroLeptonRun2/ZeroLeptonRunPeriod.h"
#include "ZeroLeptonRun2/ZeroLeptonCutVal.h"
#include "ZeroLeptonRun2/PhysObjProxyUtils.h"
#include "ZeroLeptonRun2/Counter.h"

class TTree;
class TLorentzVector;
class TVector2;

#include <string>
#include <map>


class ZeroLeptonCRWT : public cafe::Processor {
public:
  ZeroLeptonCRWT(const char *name);
  ~ZeroLeptonCRWT();
  void begin();
  void finish();
  bool processEvent(xAOD::TEvent& event);

private:
  TTree* bookTree(const std::string& name);
  TTree* getTree(const std::string& name);
  void FillCRWTVars(NTCRWTVars& crwtvars, const TLorentzVector& lepton, const TVector2& met, int lepsign);


  TTree* m_tree;
  std::string m_stringRegion;
  bool m_doSmallNtuple;
  bool m_fillTRJigsawVars;
  bool m_IsData;
  bool m_IsTruth;
  bool m_IsSignal;
  bool m_DoSystematics;
  ZeroLeptonRunPeriod m_period;
  bool m_isVR;
  bool m_isMuonChannel;
  bool m_isElectronChannel;

  NTVars m_ntv;
  NTReclusteringVars m_RTntv;
  NTExtraVars m_extrantv;
  NTRJigsawVars m_rjigsawntv;
  NTTheoryVars m_theoryntv;
  NTISRVars m_isrntv;
  NTCRWTVars m_crwtntv;

  std::string m_suffix;
  PhysObjProxyFiller* m_physobjsFiller;
  PhysObjProxyFillerTruth* m_physobjsFillerTruth;
  ZeroLeptonCutVal m_cutVal;
  PhysObjProxyUtils m_proxyUtils;
  ZeroLeptonUtils m_ZLUtils;

  Counter* m_counter;
  CounterRepository m_counterRepository;

  // TODO : an unordered_map would be faster
  std::map<std::string,TTree*> m_treeRepository;

  ZeroLeptonDerivationTag m_derivationTag;
public:
  ClassDef(ZeroLeptonCRWT,0);
};

#endif // ZeroLeptonRun2_ZeroLeptonCRWT_H_
