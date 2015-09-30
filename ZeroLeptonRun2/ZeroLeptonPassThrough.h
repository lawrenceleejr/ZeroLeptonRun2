#ifndef ZeroLeptonRun2_ZeroLeptonPassThrough_H_
#define ZeroLeptonRun2_ZeroLeptonPassThrough_H_

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

#include <string>
#include <map>


class ZeroLeptonPassThrough : public cafe::Processor {
public:
  ZeroLeptonPassThrough(const char *name);
  ~ZeroLeptonPassThrough();
  void begin();
  void finish();
  bool processEvent(xAOD::TEvent& event);

private:
  TTree* bookTree(const std::string& name);
  TTree* getTree(const std::string& name);

  TTree* m_tree;
  std::string m_stringRegion;
  bool m_doSmallNtuple;
  bool m_fillTRJigsawVars;
  bool m_fillReclusteringVars;
  bool m_doRecl;
  bool m_IsData;
  bool m_IsTruth;
  bool m_IsSignal;
  bool m_DoSystematics;
  int  m_maxJetCut;
  ZeroLeptonRunPeriod m_period;

  NTVars m_ntv;
  NTExtraVars m_extrantv;
  NTRJigsawVars m_rjigsawntv;
  NTRJigsawVars m_rjigsawntv_hlt;
  NTReclusteringVars m_RTntv;
  NTTheoryVars m_theoryntv;
  NTISRVars m_isrntv;
  bool   m_buildTriggerJetAndMET;

  std::string m_suffix;
  std::string m_suffixRecl;
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
  ClassDef(ZeroLeptonPassThrough,0);
};

#endif // ZeroLeptonRun2_ZeroLeptonPassThrough_H_
