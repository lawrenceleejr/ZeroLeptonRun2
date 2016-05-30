#ifndef PhysObjProxyFiller_h
#define PhysObjProxyFiller_h

class JetProxy;
class ElectronProxy;
class MuonProxy;
class TauProxy;
class PhotonProxy;
class JetReclusteringTool;

#include <vector>
#include <map>
#include <string>

class PhysObjProxyFiller
{
 public:
  PhysObjProxyFiller(float jetPtCut, float elPtCut, float muonPtCut, float phPtCut, const std::string suffix, bool doRecl, const std::string suffixRecl, const std::string suffixSyst, const std::string suffixSmear="");

  // fill with good and bad jets, overlap removal already performed
  void FillJetProxies(std::vector<JetProxy>& good_jets,
		      std::vector<JetProxy>& bad_jets,
		      std::vector<JetProxy>& b_jets) const;

  void FillFatJetProxies(std::vector<JetProxy>& good_fat_jets,
			 std::vector<JetProxy>& bad_fat_jets,
			 std::vector<float>& vD2_fat,
			 std::vector<bool>& visWmedium_fat) const;

  void FillJetReclProxies(std::vector<JetProxy>& good_jets_recl,
			  std::vector<float>& vD2);


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

  // fill with baseline and signal photon after overlap removal 
  void FillPhotonProxies(std::vector<PhotonProxy>& baseline_photons,
			 std::vector<PhotonProxy>& isolated_baseline_photons,
			 std::vector<PhotonProxy>& isolated_signal_photons);

  // change suffix
  void setSuffix(const std::string& suffix) {m_suffix = suffix;}
  void setSuffixSyst(const std::string& suffixSyst) {m_suffixSyst = suffixSyst;}

 private:
  float m_jetPtCut;
  float m_elPtCut;
  float m_muonPtCut;
  float m_phPtCut;
  bool m_doRecl;
  std::string m_suffix;
  std::string m_suffixRecl;
  std::string m_suffixSyst;
  std::string m_suffixSmear;
  JetReclusteringTool* m_jetReclusteringTool;
  std::map<std::string, JetReclusteringTool*> m_jrtMap;

  void setJRT();

};



#endif
