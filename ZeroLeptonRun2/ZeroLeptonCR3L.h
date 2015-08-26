#ifndef ZeroLeptonRun2_ZeroLeptonCR3L_H_
#define ZeroLeptonRun2_ZeroLeptonCR3L_H_

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


class ZeroLeptonCR3L : public cafe::Processor {
public:
  ZeroLeptonCR3L(const char *name);
  ~ZeroLeptonCR3L();
  void begin();
  void finish();
  bool processEvent(xAOD::TEvent& event);

private:
  TTree* bookTree(const std::string& name);
  TTree* getTree(const std::string& name);
  void FillCR3LVars(NTCR3LVars& cr3lvars, std::vector<TLorentzVector>& lepton, const TVector2& met, std::vector<int> lepsigns, int lepfromW, float InvMassLepPair, std::vector<float> vptvarcone20, std::vector<float> vptvarcone30, std::vector<float> vetcone20, bool isTruth);


  TTree* m_tree;
  std::string m_stringRegion;
  bool m_doSmallNtuple;
  bool m_fillTRJigsawVars;
  bool m_fillReclusteringVars;
  bool m_IsData;
  bool m_IsSignal;
  bool m_IsTruth;
  bool m_isVR;
  bool m_doRecl;
  bool m_DoSystematics;
  ZeroLeptonRunPeriod m_period;
  bool m_isMuonChannel;
  bool m_isElectronChannel;

  NTVars m_ntv;
  NTExtraVars m_extrantv;
  NTRJigsawVars m_rjigsawntv;
  NTReclusteringVars m_RTntv;
  NTTheoryVars m_theoryntv;
  NTCR3LVars m_cr3lntv;

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
  ClassDef(ZeroLeptonCR3L,0);
};

#endif // ZeroLeptonRun2_ZeroLeptonCR3L_H_
