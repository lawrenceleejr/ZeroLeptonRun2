#ifndef ZeroLeptonRun2_ZeroLeptonCRY_H_
#define ZeroLeptonRun2_ZeroLeptonCRY_H_

#include "cafe/Processor.h"
#include "ZeroLeptonRun2/ZeroLeptonNTVars.h"
#include "ZeroLeptonRun2/PhysObjProxies.h"
#include "ZeroLeptonRun2/PhysObjProxyFiller.h"
#include "ZeroLeptonRun2/PhysObjProxyFillerTruth.h"
#include "ZeroLeptonRun2/ZeroLeptonUtils.h"
#include "ZeroLeptonRun2/ZeroLeptonRunPeriod.h"
#include "ZeroLeptonRun2/ZeroLeptonCutVal.h"
#include "ZeroLeptonRun2/PhysObjProxyUtils.h"
#include "ZeroLeptonRun2/CleaningHelper.h"
#include "ZeroLeptonRun2/Counter.h"

#include "AsgTools/AsgMetadataTool.h"
#include "AsgTools/ToolHandle.h"

class TTree;
class IAsgPhotonIsEMSelector;

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
  void FillNTCRYVars(NTCRYVars& cryntv, const TLorentzVector& photon, TVector2& origmisset,int tight, int loose, float etcone20,float ptvarcone20,
		     float ptcone20,float etcone40, float ptvarcone40,float ptcone40,
		     int truthtype, int truthorigin,int isEMtight, int isSignal);

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
  ZeroLeptonRunPeriod m_period;

  NTVars m_ntv;
  NTExtraVars m_extrantv;
  NTRJigsawVars m_rjigsawntv;
  NTReclusteringVars m_RTntv;
  NTCRYVars m_cryntv;

  std::string m_suffix;
  std::string m_suffixRecl;
  std::string m_suffixSyst;
  PhysObjProxyFiller* m_physobjsFiller;
  PhysObjProxyFillerTruth* m_physobjsFillerTruth;
  CleaningHelper m_cleaningHelper;
  ZeroLeptonCutVal m_cutVal;
  PhysObjProxyUtils m_proxyUtils;
  ZeroLeptonUtils m_ZLUtils;

  Counter* m_counter;
  CounterRepository m_counterRepository;

  ToolHandle<IAsgPhotonIsEMSelector>     m_photonSelIsEM;

  // TODO : an unordered_map would be faster
  std::map<std::string,TTree*> m_treeRepository;

public:
  ClassDef(ZeroLeptonCRY,0);
};

#endif // ZeroLeptonRun2_ZeroLeptonCRY_H_
