#ifndef LeptonEfficiency_h
#define LeptonEfficiency_h

#include <string>
#include <iostream>
#include <cmath>
#include <cstdio>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>


class LeptonEfficiency {

 public:
  
  // constructor
  LeptonEfficiency();
  
  // destructor
  ~LeptonEfficiency();

  void initEfficiencyFile(const std::string inputFileName);
  void setHistName(const std::string histName){m_histName = histName;}
  void setPtMax(const Double_t ptMax){m_ptMax = ptMax;}
  void setPtMin(const Double_t ptMin){m_ptMin = ptMin;}
  void setEtaMax(const Double_t etaMax){m_etaMax = etaMax;}

  Double_t getEfficiency(const Double_t pt, const Double_t eta, const int sys=0);

 private:
  TFile *m_inputFile;
  std::string m_histName; 
  Bool_t m_histSet;
  Double_t m_ptMax;
  Double_t m_ptMin;
  Double_t m_etaMax;
  
};

Double_t Nfake(const Double_t realEff, const Double_t fakeEff, const Double_t Ntight, const Double_t Nloose);
Double_t Nfakefake(const Double_t realEff0, const Double_t realEff1, const Double_t fakeEff0, const Double_t fakeEff1,
		   const Double_t Ntt, const Double_t Ntl, const Double_t Nlt, const Double_t Nll);
Double_t Nrealfake(const Double_t realEff0, const Double_t realEff1, const Double_t fakeEff0, const Double_t fakeEff1,
		   const Double_t Ntt, const Double_t Ntl, const Double_t Nlt, const Double_t Nll);
Double_t Nfakereal(const Double_t realEff0, const Double_t realEff1, const Double_t fakeEff0, const Double_t fakeEff1,
		   const Double_t Ntt, const Double_t Ntl, const Double_t Nlt, const Double_t Nll);
Double_t Nrealreal(const Double_t realEff0, const Double_t realEff1, const Double_t fakeEff0, const Double_t fakeEff1,
		   const Double_t Ntt, const Double_t Ntl, const Double_t Nlt, const Double_t Nll);
Double_t Nrrf_rfr_frr(const Double_t r1, const Double_t r2, const Double_t r3, const Double_t f1, const Double_t f2, const Double_t f3,
		      const Double_t Nttt, const Double_t Nttl, const Double_t Ntlt, const Double_t Nltt);

#endif // LeptonEfficiency_h
