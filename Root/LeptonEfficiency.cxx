#include "ZeroLeptonRun2/LeptonEfficiency.h"
#include <TMath.h>

LeptonEfficiency::LeptonEfficiency()
{
  m_ptMax = 1000000.; // MeV
  m_ptMin = 10000.;
  m_etaMax = 2.47;
  m_histSet = kFALSE; 
}


LeptonEfficiency::~LeptonEfficiency()
{
  if( m_histSet ){
    m_inputFile->Close();
    delete m_inputFile;
  }
}


void LeptonEfficiency::initEfficiencyFile(const std::string inputFileName)
{
  if( !m_histSet ){
    m_inputFile = new TFile( inputFileName.c_str(), "READ");
    m_histSet = kTRUE;
    std::cout << "Input file " << inputFileName.c_str() << " is set" << std::endl;
  }else{
    std::cerr << "Input file was already set!" << std::endl;
  }
}

Double_t LeptonEfficiency::getEfficiency(const Double_t pt, const Double_t eta, const int sys)
{
  Double_t etaLocal = eta;
  if( eta<-m_etaMax ){
    etaLocal = -m_etaMax+0.05;
  }else if( eta>m_etaMax ){
    etaLocal = m_etaMax-0.05;
  }

  Double_t ptLocal = pt;
  if( pt>m_ptMax ) ptLocal = m_ptMax - 1.;
  if( pt<0. ) return 0.;
  if( pt<m_ptMin ) ptLocal = m_ptMin+0.5;
  
  if( m_histSet ){
    TH2F *hist = (TH2F*)m_inputFile->Get(m_histName.c_str());
    Double_t eff = hist->GetBinContent( hist->FindFixBin( ptLocal*0.001, etaLocal ) );
    if(sys!=0){
      Double_t err = hist->GetBinError( hist->FindFixBin( ptLocal*0.001, etaLocal ) );
      if(sys==1){
	eff += err;
      }else if(sys==2){
	eff -= err;
      }
    }
    delete hist;
    return eff;
  }else{
    std::cerr << "Efficiency histogram is not set" << std::endl;
    return 0.;
  }
}



Double_t Nfake(const Double_t realEff, const Double_t fakeEff, const Double_t Ntight, const Double_t Nloose)
{
  if(realEff!=fakeEff){
    return ((realEff-1)*Ntight + realEff*Nloose)*fakeEff/(realEff-fakeEff);
  }else{
    return (Ntight+Nloose)*0.5*fakeEff;
  }
}

Double_t Nfakefake(const Double_t realEff0, const Double_t realEff1, const Double_t fakeEff0, const Double_t fakeEff1,
		   const Double_t Ntt, const Double_t Ntl, const Double_t Nlt, const Double_t Nll)
{
  if(realEff0!=fakeEff0 && realEff1!=fakeEff1){
    return fakeEff0*fakeEff1/(realEff0-fakeEff0)/(realEff1-fakeEff1)*((1.-realEff0)*(1.-realEff1)*Ntt+(realEff0-1.)*realEff1*Ntl+realEff0*(realEff1-1.)*Nlt+realEff0*realEff1*Nll); 
  }else{
    return 0.;
  }
}

Double_t Nrealfake(const Double_t realEff0, const Double_t realEff1, const Double_t fakeEff0, const Double_t fakeEff1,
		   const Double_t Ntt, const Double_t Ntl, const Double_t Nlt, const Double_t Nll)
{
  if(realEff0!=fakeEff0 && realEff1!=fakeEff1){
    return realEff0*fakeEff1/(realEff0-fakeEff0)/(realEff1-fakeEff1)*((fakeEff0-1.)*(1-realEff1)*Ntt+(1.-fakeEff0)*realEff1*Ntl+fakeEff0*(1-realEff1)*Nlt-fakeEff0*realEff1*Nll);
  }else{
    return 0.;
  }  
}

Double_t Nfakereal(const Double_t realEff0, const Double_t realEff1, const Double_t fakeEff0, const Double_t fakeEff1,
		   const Double_t Ntt, const Double_t Ntl, const Double_t Nlt, const Double_t Nll)
{
  if(realEff0!=fakeEff0 && realEff1!=fakeEff1){
    return fakeEff0*realEff1/(realEff0-fakeEff0)/(realEff1-fakeEff1)*((realEff0-1.)*(1-fakeEff1)*Ntt+(1.-realEff0)*fakeEff1*Ntl+realEff0*(1.-fakeEff1)*Nlt-realEff0*fakeEff1*Nll);
  }else{
    return 0.;
  }
}

Double_t Nrealreal(const Double_t realEff0, const Double_t realEff1, const Double_t fakeEff0, const Double_t fakeEff1,
		   const Double_t Ntt, const Double_t Ntl, const Double_t Nlt, const Double_t Nll)
{
  if(realEff0!=fakeEff0 && realEff1!=fakeEff1){
    return realEff0*realEff1/(realEff0-fakeEff0)/(realEff1-fakeEff1)*((1-fakeEff0)*(1-fakeEff1)*Ntt+(fakeEff0-1.)*fakeEff1*Ntl+fakeEff0*(fakeEff1-1.)*Nlt+fakeEff0*fakeEff1*Nll);
  }else{
    return 0.;
  }
}

Double_t Nrrf_rfr_frr(const Double_t r1, const Double_t r2, const Double_t r3, const Double_t f1, const Double_t f2, const Double_t f3,
		      const Double_t Nttt, const Double_t Nttl, const Double_t Ntlt, const Double_t Nltt)
{
  Double_t Nfake = 0.;
  if(r3!=f3){
    Nfake += (-(1.-r3)*Nttt + r3*Nttl)*f3/(r3-f3);
  }
  if(r2!=f3){
    Nfake += (-(1-r2)*Nttt + r2*Ntlt)*f2/(r2-f2);
  }
  if(r1!=f1){
    Nfake += (-(1-r1)*Nttt + r1*Nltt)*f1/(r1-f1);
  }

  return Nfake;
}
