#ifndef ZeroLeptonRun2_ZeroLeptonCRZ_H_
#define ZeroLeptonRun2_ZeroLeptonCRZ_H_

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


class ZeroLeptonCRZ : public cafe::Processor {
public:
  ZeroLeptonCRZ(const char *name);
  ~ZeroLeptonCRZ();
  void begin();
  void finish();
  bool processEvent(xAOD::TEvent& event);

private:
  TTree* bookTree(const std::string& name);
  TTree* getTree(const std::string& name);
  void FillCRZVars(NTCRZVars& crzvars, std::vector<TLorentzVector>& lepton, const TVector2& met, std::vector<int> lepsigns);


  TTree* m_tree;
  std::string m_stringRegion;
  bool m_fillReclusteringVars;
  bool m_doRecl;
  bool m_doSmallNtuple;
  bool m_fillTRJigsawVars;
  bool m_IsData;
  bool m_IsSignal;
  bool m_IsTruth;
  bool m_DoSystematics;
  ZeroLeptonRunPeriod m_period;
  bool m_isMuonChannel;
  bool m_isElectronChannel;
  float m_PtLepton2;

  NTVars m_ntv;
  NTExtraVars m_extrantv;
  NTRJigsawVars m_rjigsawntv;
  NTReclusteringVars m_RTntv;
  NTTheoryVars m_theoryntv;
  NTCRZVars m_crzntv;

  std::string m_suffix;
  std::string m_suffixRecl;
  std::string m_suffixSyst;
  PhysObjProxyFiller* m_physobjsFiller;
  PhysObjProxyFillerTruth* m_physobjsFillerTruth;
  ZeroLeptonCutVal m_cutVal;
  PhysObjProxyUtils m_proxyUtils;
  ZeroLeptonUtils m_ZLUtils;

  Counter* m_counter;
  CounterRepository m_counterRepository;

  // TODO : an unordered_map would be faster
  std::map<std::string,TTree*> m_treeRepository;

public:
  ClassDef(ZeroLeptonCRZ,0);
};

#endif // ZeroLeptonRun2_ZeroLeptonCRZ_H_
