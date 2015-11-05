#ifndef PhysObjProxyFillerTruth_h
#define PhysObjProxyFillerTruth_h

class JetProxy;
class ElectronTruthProxy;
class MuonTruthProxy;
class PhotonTruthProxy;

#include <vector>
#include <string>

class PhysObjProxyFillerTruth
{
 public:
  PhysObjProxyFillerTruth(float jetPtCut, float elPtCut, float muonPtCut, float phPtCut, const std::string suffix);

  // fill with good and bad jets, overlap removal already performed
  void FillJetProxies(std::vector<JetProxy>& good_jets,
		      std::vector<JetProxy>& b_jets) const;

  // fill with baseline and signal electron after overlap removal
  void FillElectronProxies(std::vector<ElectronTruthProxy>& baseline_electrons,
			   std::vector<ElectronTruthProxy>& isolated_baseline_electrons,
			   std::vector<ElectronTruthProxy>& isolated_signal_electrons);

  // fill with baseline and signal muon after overlap removal
  void FillMuonProxies(std::vector<MuonTruthProxy>& baseline_muons,
		       std::vector<MuonTruthProxy>& isolated_baseline_muons,
		       std::vector<MuonTruthProxy>& isolated_signal_muons);

  // fill with baseline and signal photon after overlap removal
  void FillPhotonProxies(std::vector<PhotonTruthProxy>& baseline_photons,
			 std::vector<PhotonTruthProxy>& isolated_baseline_photons,
			 std::vector<PhotonTruthProxy>& isolated_signal_photons);


 private:
  float m_jetPtCut;
  float m_elPtCut;
  float m_muonPtCut;
  float m_phPtCut;
  std::string m_suffix;
};



#endif
