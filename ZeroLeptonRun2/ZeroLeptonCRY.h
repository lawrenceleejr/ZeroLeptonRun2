#ifndef ZeroLeptonRun2_ZeroLeptonCRY_H_
#define ZeroLeptonRun2_ZeroLeptonCRY_H_

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


class ZeroLeptonCRY : public cafe::Processor {
public:
  ZeroLeptonCRY(const char *name);
  ~ZeroLeptonCRY();
  void begin();
  void finish();
  bool processEvent(xAOD::TEvent& event);

private:
  TTree* bookTree(const std::string& name);
  TTree* getTree(const std::string& name);

  TTree* m_tree;
  std::string m_stringRegion;
  bool m_doSmallNtuple;
  bool m_IsData;
  bool m_IsTruth;
  bool m_IsSignal;
  bool m_UseSystematics;
  ZeroLeptonRunPeriod m_period;

  NTVars m_ntv;
  NTExtraVars m_extrantv;
  NTReclusteringVars m_RTntv;
  NTCRYVars m_cryntv;

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
  ClassDef(ZeroLeptonCRY,0);
};

#endif // ZeroLeptonRun2_ZeroLeptonCRY_H_
