#ifndef PhysObjProxyFiller_h
#define PhysObjProxyFiller_h

class JetProxy;
class ElectronProxy;
class MuonProxy;
class TauProxy;

#include <vector>
#include <string>

class PhysObjProxyFiller
{
 public:
  PhysObjProxyFiller(float jetPtCut, float elPtCut, float muonPtCut, const std::string suffix);

  // fill with good and bad jets, overlap removal already performed
  void FillJetProxies(std::vector<JetProxy>& good_jets,
		      std::vector<JetProxy>& bad_jets,
		      std::vector<JetProxy>& b_jets) const;

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
 private:
  float m_jetPtCut;
  float m_elPtCut;
  float m_muonPtCut;
  std::string m_suffix;
};



#endif
