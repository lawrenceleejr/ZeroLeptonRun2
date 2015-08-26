#ifndef PhysObjProxyFiller_h
#define PhysObjProxyFiller_h

class JetProxy;
class ElectronProxy;
class MuonProxy;
class TauProxy;
class JetReclusteringTool;

#include <vector>
#include <string>

class PhysObjProxyFiller
{
 public:
  PhysObjProxyFiller(float jetPtCut, float elPtCut, float muonPtCut, const std::string suffix, bool doRecl, const std::string suffixRecl);

  // fill with good and bad jets, overlap removal already performed
  void FillJetProxies(std::vector<JetProxy>& good_jets,
		      std::vector<JetProxy>& bad_jets,
		      std::vector<JetProxy>& b_jets) const;

  void FillJetReclProxies(std::vector<JetProxy>& good_jets_recl,
			  std::vector<float>& vD2) const;


  // fill with baseline and signal electron after overlap removal
  void FillElectronProxies(std::vector<ElectronProxy>& baseline_electrons,
			   std::vector<ElectronProxy>& isolated_baseline_electrons,
			   std::vector<ElectronProxy>& isolated_signal_electrons);

  // fill with baseline and signal muon after overlap removal
  void FillMuonProxies(std::vector<MuonProxy>& baseline_muons,
		       std::vector<MuonProxy>& isolated_baseline_muons,
		       std::vector<MuonProxy>& isolated_signal_muons);

  // fill with baseline and signal tau
  void FillTauProxies(std::vector<TauProxy>& baseline_taus,
		      std::vector<TauProxy>& signal_taus);

  // change suffix
  void setSuffix(const std::string& suffix) {m_suffix = suffix;}

 private:
  float m_jetPtCut;
  float m_elPtCut;
  float m_muonPtCut;
  bool m_doRecl;
  std::string m_suffix;
  std::string m_suffixRecl;
  JetReclusteringTool* m_jetReclusteringTool;
};



#endif
