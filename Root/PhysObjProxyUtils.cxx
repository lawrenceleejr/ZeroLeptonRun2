#include "ZeroLeptonRun2/PhysObjProxyUtils.h"
#include "ZeroLeptonRun2/PhysObjProxies.h"
#include "ZeroLeptonRun2/ZeroLeptonNTVars.h"
#include "xAODJet/Jet.h"
#include "TLorentzVector.h"
//#include "OxbridgeKinetics/OKTwoParents.h"

#include "ZeroLeptonRun2/Sphericity.h"

#include <cmath>
#include <iostream>

#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"


#include "RestFrames/RestFrame.hh"
#include "RestFrames/RFrame.hh"
#include "RestFrames/RLabFrame.hh"
#include "RestFrames/RDecayFrame.hh"
#include "RestFrames/RVisibleFrame.hh"
#include "RestFrames/RInvisibleFrame.hh"
#include "RestFrames/RSelfAssemblingFrame.hh"
#include "RestFrames/InvisibleMassJigsaw.hh"
#include "RestFrames/InvisibleRapidityJigsaw.hh"
#include "RestFrames/ContraBoostInvariantJigsaw.hh"
#include "RestFrames/MinimizeMassesCombinatoricJigsaw.hh"
#include "RestFrames/InvisibleGroup.hh"
#include "RestFrames/CombinatoricGroup.hh"
#include "RestFrames/FramePlot.hh"

void PhysObjProxyUtils::EnergyWeightedTime(const std::vector<JetProxy>& jets, 
					   std::vector<float>& time) const
{
  time.resize(5,-999.f);
  
  double denom = 0.;
  double num = 0.;
  std::size_t maxi = jets.size();
  if ( maxi > 6 ) maxi = 6;
  for ( std::size_t i = 0; i < maxi; i++ ) {
    if ( jets[i].jet() ) {  // skip leptons 
      denom = denom + jets[i].E();
      float time = -99999.f;
      jets[i].jet()->getAttribute(xAOD::JetAttribute::Timing,time);
      num = num + jets[i].E() * time;
    }    
    if ( i == 1 ) time[0] = num/denom; // 2jets
    else if ( i == 2 ) time[1] = num/denom; // 3jets
    else if ( i == 3 ) time[2] = num/denom; // 4jets
    else if ( i == 4 ) time[3] = num/denom; //5jets
    else if ( i == 5 ) time[4] = num/denom; //6jets
  }

  /*
  std::cout << " jet timing ";
  for ( std::size_t i = 0; i <5; ++i ) std::cout << " " << time[i];
  std::cout << std::endl;
  */
}

double PhysObjProxyUtils::SmallestdPhi(const std::vector<JetProxy>& jets, double met_phi) const
{
  if ( jets.size() < 2 ) return 999.;
  double dphi1 = std::acos(std::cos(jets[0].Phi() - met_phi));
  double dphi2 = std::acos(std::cos(jets[1].Phi() - met_phi));
  double dphi3 = 999.;
  if ( jets.size() > 2 && jets[2].Pt() > 40000. ) { 
    dphi3= std::acos(std::cos(jets[2].Phi() - met_phi));
  }
  double min1 = std::min(dphi1,dphi2);
  return std::min(min1,dphi3);
}

double PhysObjProxyUtils::SmallestRemainingdPhi(const std::vector<JetProxy>& jets, double met_phi) const
{
  double remainingDPhi=999;
  double dphiMin=999;
  unsigned int jetcount = 0;
  for ( std::vector<JetProxy>::const_iterator itjet = jets.begin();
	  itjet != jets.end(); ++itjet )
    {      
      jetcount++;
      if ( jetcount>3 && itjet->Pt()>40000 ) {
	remainingDPhi = std::acos(std::cos(itjet->Phi() - met_phi));
	dphiMin = std::min(remainingDPhi,dphiMin);
      }
    }
  return dphiMin;
}


double PhysObjProxyUtils::Meff(const std::vector<JetProxy>& jets, size_t njets, double met, double jetPtCut, double extraJetPtCut)
{
  double meff=met;
  if ( jets.size() < njets ) njets=jets.size();
  for(size_t i=0; i<njets; i++) 
  {
    if ( i<=3 && jets[i].Pt() > jetPtCut ) meff += jets[i].Pt();
    if ( i>3  && jets[i].Pt() > extraJetPtCut ) meff += jets[i].Pt();
  }
  return meff;
}


void PhysObjProxyUtils::ComputeSphericity(const std::vector<JetProxy>& jets, double & Sp, double & ST, double & Ap)
{
  
  int njet = jets.size(); //you can set Number of jets for calculation
  Sp=-1;
  ST=-1;
  Ap=-1;
  
  vector<TLorentzVector> v_tlv;
 
  //prepare vector<TLorentzVector> of jets to use
  for(size_t ijet=0; ijet<jets.size(); ijet++) 
    {
      
      //      if(jets[ijet].Pt()<ptcut)break;
      TLorentzVector jet;
      jet.SetPtEtaPhiM(jets[ijet].Pt(),
		       jets[ijet].Eta(),
		       jets[ijet].Phi(),
		       jets[ijet].M());
      v_tlv.push_back(jet);
    }
  
    if(v_tlv.size() < (size_t)njet || v_tlv.size()==0)return ;


    Sphericity sp; //construct
    sp.SetTLV(v_tlv, njet);
    sp.GetSphericity(Sp, ST, Ap);


};

void PhysObjProxyUtils::RJigsawVariables( const std::vector<JetProxy>& jets,
  std::map<TString,double> RJigsawVars
  ){}


void PhysObjProxyUtils::RazorVariables(const std::vector<JetProxy>& jets, 
					 Double_t metx,
					 Double_t mety, 
					 double &gaminvRp1 ,
					 double &shatR ,
					 double &mdeltaR ,
					 double &cosptR ,
					 double &Minv2 ,
					 double &Einv ,
					 double & gamma_R,
					 double &dphi_BETA_R , 
					 double &dphi_leg1_leg2 , 
					 double &costhetaR ,
					 double &dphi_BETA_Rp1_BETA_R,
					 double &gamma_Rp1,
					 double &Eleg1,
					 double &Eleg2, 
					 double &costhetaRp1)
{
  if ( jets.size() < 2 ) {
    gaminvRp1 = -999.;
    shatR = -999.;
    mdeltaR = -999.;
    cosptR = -999.;
    Minv2 = -999.;
    Einv = -999.;
    gamma_R= -999.;
    dphi_BETA_R = -999.; 
    dphi_leg1_leg2 = -999.; 
    costhetaR = -999.;
    dphi_BETA_Rp1_BETA_R= -999.;
    gamma_Rp1= -999.;
    Eleg1= -999.;
    Eleg2= -999.; 
    costhetaRp1 = -999.;
    return;
  }


  //=============================================================
  // Step 1: make megajet
  //=============================================================
  //This code is adapted from a code from CMS
  //https://twiki.cern.ch/twiki/bin/view/CMSPublic/RazorLikelihoodHowTo  
 
  //To minimize the change in the code, the vector of JetProxy is converted to a vector a TLorentzVector
  std::vector<TLorentzVector> myjets;
  for(size_t ijet=0; ijet<jets.size(); ijet++) 
    {
      TLorentzVector jet;
      jet.SetPtEtaPhiM(jets[ijet].Pt(),
		       jets[ijet].Eta(),
		       jets[ijet].Phi(),
		       jets[ijet].M());
      myjets.push_back(jet);
    }



  //build megajets
  //  vector<TLorentzVector> mynewjets;
  TLorentzVector J1, J2;
  //  bool foundGood = false;
  size_t N_comb = 1;
  //for(size_t i = 0; i < myjets.size(); i++)    
  for(size_t i = 0; i < myjets.size() && i<15; i++)//code very slow if there are many jets
    {
      N_comb *= 2;
    }
  double M_min = 9999999999999999999999999.0;
  int j_count;
  for(size_t i = 1; i < N_comb-1; i++)
  //  for(size_t i = 1; i+1 < N_comb; i++)
    {
      TLorentzVector  j_temp1, j_temp2;
      int itemp = i;
      j_count = N_comb/2;
      int count = 0;
      
      while(j_count > 0)
	{
	  
	  TLorentzVector TLorentzJets_count = myjets[count];      
	  if(itemp/j_count == 1)
	    {
	      j_temp1 += TLorentzJets_count;
	    } 
	  else 
	    {
	      j_temp2 += TLorentzJets_count;
	    }
	  
	      itemp -= j_count*(itemp/j_count);
	      j_count /= 2;
	      count++;
	}
      
      double M_temp = j_temp1.M2()+j_temp2.M2();    
      // cout << j_temp1.M2()<< " " << j_temp2.M2()<< " "<< M_temp << " "  << endl;
      
      // smallest mass
      if(M_temp < M_min)
	{
	  M_min = M_temp;
	  J1 = j_temp1;
	  J2 = j_temp2;
	}        
    }  
  
  if(J2.Pt() > J1.Pt())
    {
      TLorentzVector temp = J1;
      J1 = J2;
      J2 = temp;
    }
  //  mynewjets.push_back(J1);
  //  mynewjets.push_back(J2);
  

  //=============================================================
  // Step 2: compute superrazor variables
  //=============================================================
  //based on code provided privately by L. Lee
    
  TVector3 MET(metx, mety, 0.0);
  
  J1.SetVectM(J1.Vect(),0.0);
  J2.SetVectM(J2.Vect(),0.0);
  
  TVector3 vBETA_z = (1./(J1.E()+J2.E()))*(J1+J2).Vect();
  vBETA_z.SetX(0.0);
  vBETA_z.SetY(0.0);
  

  //transformation from lab frame to approximate rest frame along beam-axis
  J1.Boost(-vBETA_z);
  J2.Boost(-vBETA_z);

  TVector3 pT_CM = (J1+J2).Vect() + MET;
  pT_CM.SetZ(0.0); //should be redundant...
  
  Minv2 = (J1+J2).M2();
  Einv = sqrt(MET.Mag2()+Minv2);
  
  //////////////////////
  // definition of shatR
  //////////////////////
  TLorentzVector J1J2 = J1+J2;

  shatR = sqrt( ((J1J2).E()+Einv)*((J1J2).E()+Einv) - pT_CM.Mag2() );
  
  TVector3 vBETA_R = (1./sqrt(pT_CM.Mag2() + shatR*shatR))*pT_CM;
  gamma_R = 1./sqrt(1.-vBETA_R.Mag2());
  
  
  //transformation from lab frame to R frame
  J1.Boost(-vBETA_R);
  J2.Boost(-vBETA_R);
  

  dphi_BETA_R = ((J1J2).Vect()).DeltaPhi(vBETA_R);
  dphi_leg1_leg2 = J1.Vect().DeltaPhi(J2.Vect());
  costhetaR =  fabs((J1J2).Vect().Dot(vBETA_R)/((J1J2).Vect().Mag()*vBETA_R.Mag()));


  /////////////
  //
  // R-frame
  //
  /////////////
  
  TVector3 vBETA_Rp1 = (1./(J1.E()+J2.E()))*(J1.Vect() - J2.Vect());
  
  ////////////////////////
  // definition of gaminvRp1
  ////////////////////////
  gaminvRp1 = sqrt(1.-vBETA_Rp1.Mag2() );
  gamma_Rp1 = 1/sqrt(1.-vBETA_Rp1.Mag2() );
  dphi_BETA_Rp1_BETA_R = vBETA_Rp1.DeltaPhi(vBETA_R);

  //transformation from R frame to R+1 frames
  J1.Boost(-vBETA_Rp1);
  J2.Boost(vBETA_Rp1);
  //////////////
  //
  // R+1-frames
  //
  //////////////

  ///////////////////////
  // definition of mdeltaR
  ////////////////////////
  mdeltaR = J1.E()+J2.E();

  ///////////////////////
  // definition of cosptR
  ////////////////////////
  cosptR = pT_CM.Mag()/sqrt(pT_CM.Mag2()+mdeltaR * mdeltaR);


  Eleg1 = J1.E();
  Eleg2 = J2.E();
  costhetaRp1 = fabs(J1.Vect().Dot(vBETA_Rp1)/(J1.Vect().Mag()*vBETA_Rp1.Mag()));





 }

/*
double PhysObjProxyUtils::MT2(const std::vector<JetProxy>& jets,const TVector2& MissingET) const
{
  double mT2=0;
  if (jets.size()<2) return -1;
  
  std::vector<TLorentzVector> jets_tmp; 
  for (size_t i = 0 ; i <jets.size(); i++ ) {
    jets_tmp.push_back(TLorentzVector(*dynamic_cast<const TLorentzVector*>(&(jets[i]))));
  }

  OxbridgeKinetics::OKTwoParents ok2p;
  TVector2 ptmiss(MissingET.X(),MissingET.Y());
  double Minvis = 0.;
  ok2p.addVis(jets_tmp[0],1);
  ok2p.addVis(jets_tmp[1],2);
  ok2p.setPtMiss(ptmiss);
  ok2p.setMinvis(Minvis);
  
  mT2 = ok2p.calcM2T();
    
  return mT2;    
}
*/

void PhysObjProxyUtils::GetAlphaISRVar(const std::vector<JetProxy>& jets, double met, std::vector<double>& alpha_vec) const
{
  alpha_vec.clear();
  for (size_t ijet=0; ijet<jets.size() ; ijet++) {
    if ( jets[ijet].Pt() < 50000.) continue;
    double alpha = std::min(jets[ijet].Pt(),met)/std::max(jets[ijet].Pt(),met);
    alpha_vec.push_back(alpha);
  }
  return; 
}

void PhysObjProxyUtils::GetMinPtDistinctionISR(const std::vector<JetProxy>& jets, std::vector<double>& minPtDistinction_vec) const
{
  minPtDistinction_vec.clear();
  for (size_t isrcandidate=0; isrcandidate<jets.size() ; isrcandidate++) {
    if ( jets[isrcandidate].Pt() < 50000.) continue;
    double minDist=99.;
    for (size_t jjet=0; jjet<jets.size() ; jjet++) {
      if (jjet==isrcandidate) continue;
      double max2jet_pt=std::max(jets[jjet].Pt(),jets[isrcandidate].Pt());
      double min2jet_pt=std::min(jets[jjet].Pt(),jets[isrcandidate].Pt());
      double fraction=max2jet_pt/min2jet_pt;
      minDist=std::min(minDist,fraction);
    }
    minPtDistinction_vec.push_back(minDist);
  }
  return; 

}

void PhysObjProxyUtils::GetMinDeltaFraction(const std::vector<JetProxy>& jets, std::vector<double>& minDeltaFrac_vec) const
{
  minDeltaFrac_vec.clear();
  for (size_t isrcandidate=0; isrcandidate<jets.size() ; isrcandidate++) {
    if (jets[isrcandidate].Pt() < 50000.) continue;
    double delta_isrcand = jets[isrcandidate].M()/jets[isrcandidate].Pt();
    double minDeltaFrac=99.;
    for (size_t jjet=0; jjet<jets.size() ; jjet++) {
      if (jjet==isrcandidate) continue;
      double delta_jjet = jets[jjet].M()/jets[jjet].Pt();
      double mindelta = std::min(delta_isrcand,delta_jjet);
      double maxdelta = std::max(delta_isrcand,delta_jjet);
      double fraction = maxdelta/mindelta;
      minDeltaFrac = std::min(minDeltaFrac,fraction);
    }
    minDeltaFrac_vec.push_back(minDeltaFrac);
  }
  return; 
}


void PhysObjProxyUtils::GetMinRapidityGap(const std::vector<JetProxy>& jets, std::vector<double>& minRapGap_vec) const
{
  minRapGap_vec.clear();
  for (size_t isrcandidate=0; isrcandidate<jets.size() ; isrcandidate++) {
    if (jets[isrcandidate].Pt() < 50000.) continue;
    double minRapGap=99.;
    for (size_t jjet=0; jjet<jets.size() ; jjet++) {
      if (jjet==isrcandidate) continue;
      double diffrap = std::abs(jets[isrcandidate].Eta()-jets[jjet].Eta());
      minRapGap = std::min(minRapGap,diffrap);
    }
    minRapGap_vec.push_back(minRapGap);
  }
  return; 
}

void PhysObjProxyUtils::GetMaxRapidityOtherJets(const std::vector<JetProxy>& jets, std::vector<double>& maxRapOtherJets_vec) const
{
  maxRapOtherJets_vec.clear();
  for (size_t isrcandidate=0; isrcandidate<jets.size() ; isrcandidate++) {
    if (jets[isrcandidate].Pt() < 50000.) continue;
    double maxRap=-99.;
    for (size_t jjet=0; jjet<jets.size() ; jjet++) {
      if (jjet==isrcandidate) continue;
      maxRap = std::max(maxRap,jets[jjet].Eta());
    }
    maxRapOtherJets_vec.push_back(maxRap);
  }
  return; 
}

void PhysObjProxyUtils::GetdPhiJetMet(const std::vector<JetProxy>& jets, double met_phi, std::vector<double>& dPhiJetMet_vec) const
{
  dPhiJetMet_vec.clear();
  for (size_t isrcandidate=0; isrcandidate<jets.size() ; isrcandidate++) {
    if (jets[isrcandidate].Pt() < 50000.) continue;
    double dphi = std::abs(jets[isrcandidate].Phi()-met_phi);
    dPhiJetMet_vec.push_back(dphi);
  }
  return; 
}

void PhysObjProxyUtils::GetISRJet(const std::vector<JetProxy>& jets,
				  std::vector<size_t>& isr_jet_indices,
				  double met,
				  double phi_met,
				  std::string signal, 
				  bool usealpha) const
{
  using std::cout;
  using std::endl;

  std::vector<size_t> tmp_isr_jet_indices;
  std::vector<bool> maxpt_crit_vector;
  std::vector<bool> maxdelta_crit_vector;
  std::vector<bool> rapid_crit_vector; 
  if (!(signal=="gluino" || signal=="squark")) { 
    cout << " [GetISRJet]: a signal type other than 'squark' or 'gluino' is requested. Please use one of those two. " << endl;
    cout << " Exiting without tagging an ISR jet " << endl; 
  } 


  for (size_t isrcandidate=0; isrcandidate<jets.size() ; isrcandidate++) {
    if (jets[isrcandidate].Pt() < 50000.) continue;  // Only consider jets with Pt > 50 GeV
    double delta_isrjet = jets[isrcandidate].M()/jets[isrcandidate].Pt();
    bool pretag_crit=true;
    bool maxpt_crit=true;
    bool rapid_crit=true;
    bool maxdelta_crit=true;
    for (size_t jjet=0; jjet<jets.size() ; jjet++) { 
      if (jjet==isrcandidate) continue;
      
      double max2jet_pt = std::max(jets[jjet].Pt(),jets[isrcandidate].Pt());
      double min2jet_pt = std::min(jets[jjet].Pt(),jets[isrcandidate].Pt());
      double diffrap = std::abs(jets[isrcandidate].Eta()-jets[jjet].Eta());
      double delta_jjet = jets[jjet].M()/jets[jjet].Pt();
      double max2jet_delta = std::max(delta_isrjet,delta_jjet);
      double min2jet_delta = std::min(delta_isrjet,delta_jjet);
      
      // pretag criteria: 

      if (max2jet_pt/min2jet_pt<=2.0) maxpt_crit = false;
      if (diffrap<=1.0) rapid_crit = false;
      if (max2jet_delta/min2jet_delta<=1.5) maxdelta_crit = false;
    } 
    if (signal=="gluino") { 
      if (maxpt_crit || rapid_crit || maxdelta_crit) pretag_crit = true;
    }
    if (signal=="squark") { 
      if ( rapid_crit) pretag_crit = true;
    }
    if (!pretag_crit) continue;
      
    // Only isr candidate jets which have survived pre-tag conditions
    bool jjetrapidity_crit = true;
    bool isr_rapidity_diff_crit = true;
    for (size_t jjet = 0; jjet<jets.size() ; jjet++) { 
      if (jjet==isrcandidate) continue;
      if (std::abs(jets[jjet].Eta())>=2.0) jjetrapidity_crit = false;
      if (std::abs(jets[isrcandidate].Eta()-jets[jjet].Eta())<=0.5) isr_rapidity_diff_crit = false;
    } 
    bool isrjetrapidity_crit = false;
    if (std::abs(jets[isrcandidate].Eta())>1.0) isrjetrapidity_crit = true;
    bool deltaphi_crit = false;
    double dphi = std::abs(jets[isrcandidate].Phi()-phi_met);
    if (dphi>=TMath::Pi()) dphi = 2*TMath::Pi()-dphi;
    if (dphi>2.0) deltaphi_crit = true;
    bool alpha_crit = false;
    double alpha = std::min(jets[isrcandidate].Pt(),met)/std::max(jets[isrcandidate].Pt(),met);
    if (alpha>0.4) alpha_crit = true;
    
    // tagging:
    if (signal=="gluino") { 
      if (jjetrapidity_crit && isrjetrapidity_crit && isr_rapidity_diff_crit && deltaphi_crit && ( ( alpha_crit && usealpha ) || (!usealpha))) {
	tmp_isr_jet_indices.push_back(isrcandidate);
	
	maxpt_crit_vector.push_back(maxpt_crit);
	maxdelta_crit_vector.push_back(maxdelta_crit);
	rapid_crit_vector.push_back(rapid_crit);
      }
    } 
    if (signal=="squark") { 
      if (jjetrapidity_crit && isrjetrapidity_crit && deltaphi_crit && ( ( alpha_crit && usealpha ) || (!usealpha))) {
	tmp_isr_jet_indices.push_back(isrcandidate);
	rapid_crit_vector.push_back(rapid_crit);
      }
    } 
  }


  if (tmp_isr_jet_indices.size() > 1) { 
    if (signal=="gluino") {
      std::vector<int> which_criteria;
      std::vector<size_t> whichjet_crit1;
      std::vector<size_t> whichjet_crit2;
      std::vector<size_t> whichjet_crit3;
      for (size_t ijet=0; ijet <tmp_isr_jet_indices.size() ; ijet++) { 
	if (maxpt_crit_vector[ijet]==true) whichjet_crit1.push_back(tmp_isr_jet_indices[ijet]);
	if (maxpt_crit_vector[ijet]==false && rapid_crit_vector[ijet]==true)  whichjet_crit1.push_back(tmp_isr_jet_indices[ijet]);
	if (maxpt_crit_vector[ijet]==false && rapid_crit_vector[ijet]==false && maxdelta_crit_vector[ijet]==true)  whichjet_crit1.push_back(tmp_isr_jet_indices[ijet]);
      }
      if (whichjet_crit1.size()==1) isr_jet_indices.push_back(whichjet_crit1[0]);
      else if (whichjet_crit1.size()==0 && whichjet_crit2.size()==1 ) isr_jet_indices.push_back(whichjet_crit2[0]);
      else if (whichjet_crit1.size()>1 ) {
	double maxpt = 0;
	size_t thisjet = 0;
	for (size_t ijet=0; ijet < whichjet_crit1.size() ; ijet++){
	  if (jets[whichjet_crit1[ijet]].Pt() > maxpt){
	    maxpt =  jets[whichjet_crit1[ijet]].Pt();
	    thisjet = ijet;
	  }
	}
	isr_jet_indices.push_back(whichjet_crit1[thisjet]);

	return;
      } 

      else if (whichjet_crit1.size()==0 && whichjet_crit2.size()>1 )   return;
      else if (whichjet_crit1.size()==0 && whichjet_crit2.size()==0 && whichjet_crit3.size()>1 ) { 
	double maxdelta = 0;
	size_t thisjet = 0;
	for (size_t ijet=0; ijet < whichjet_crit3.size() ; ijet++){
	  if (jets[whichjet_crit3[ijet]].M()/jets[whichjet_crit3[ijet]].Pt()> maxdelta) {
	    maxdelta = jets[whichjet_crit3[ijet]].M()/jets[whichjet_crit3[ijet]].Pt();
	    thisjet = ijet;
	  }
	}
	isr_jet_indices.push_back(whichjet_crit3[thisjet]);

	return;
      }
    }
    
    else if (signal=="squark") { 
      // if more than 1 jet passes these selections: none of them is tagged as an ISR jet! 
      
      return;
    }
  }
  else if (tmp_isr_jet_indices.size() == 1) { 
    isr_jet_indices = tmp_isr_jet_indices;
  }
  return;
}


bool PhysObjProxyUtils::CosmicMuon(const std::vector<MuonProxy>& muons) const
{
  for ( size_t i = 0; i < muons.size(); ++i ) {
    if ( muons[i].isCosmic() ) return true;
  }
  return false;
}


bool PhysObjProxyUtils::isbadMETmuon(const std::vector<MuonProxy>& muons,
				     float MET, const TVector2& MissingET) const
{
  bool isbadmetmuon=false; 
  TVector2 MissingETMuon(0,0);
  for ( size_t iMu = 0; iMu < muons.size(); ++iMu ) {
    MissingETMuon -= muons[iMu].Vect().XYvector();
  }
  double METMuon = MissingETMuon.Mod();
  
  double MET_muon_ratio = METMuon/MET*std::cos(MissingETMuon.Phi()-MissingET.Phi()) ;
  if(MET_muon_ratio>=0.5) isbadmetmuon=true; 

  return isbadmetmuon; 
}


bool PhysObjProxyUtils::badTileVeto(const std::vector<JetProxy>& jets, const TVector2& MissingET) const
{
  bool isDeadTile=false;   
  for ( std::vector<JetProxy>::const_iterator itjet = jets.begin();
	itjet != jets.end(); itjet++ ) {     
    if ( ! (*itjet).jet() ) continue;
    double jet_pt = (*itjet).Pt();
    if(jet_pt<40000.) continue;
    double jet_phi = (*itjet).Phi();
    double  jet_BCH_CORR_JET = 0.;
    (*itjet).jet()->getAttribute(xAOD::JetAttribute::BchCorrJet,jet_BCH_CORR_JET);
    double MET_phi = MissingET.Phi();
    if(std::acos(std::cos( jet_phi-MET_phi )  )<0.3 && jet_BCH_CORR_JET>0.05) isDeadTile=true;
  }
  return isDeadTile ;
}

bool PhysObjProxyUtils::chfVeto(const std::vector<JetProxy>& jets) const
{
  bool shouldbecleaned=false;
  for ( std::size_t i = 0; i < std::min((std::size_t)2,jets.size()); ++i ) {
    if ( !jets[i].jet() ) continue; // skip lepton disguised as jet

    std::vector<float> sumPtTrk;
    jets[i].jet()->getAttribute(xAOD::JetAttribute::SumPtTrkPt500,sumPtTrk); // FIXME or SumPtTrkPt1000 ??
    float chf = sumPtTrk[0]/jets[i].Pt();

    float emf = 0.;
    jets[i].jet()->getAttribute(xAOD::JetAttribute::EMFrac,emf);

    if ( jets[i].Pt()>100000. && chf<0.02 && std::abs(jets[i].Eta())<2. ) shouldbecleaned=true;

    if ( jets[i].Pt()>100000. && chf<0.05 && std::abs(jets[i].Eta())<2. && emf>0.9) shouldbecleaned=true;
  }
  return shouldbecleaned;
}

bool PhysObjProxyUtils::chfTileVeto(const std::vector<JetProxy>& jets) const
{
  bool shouldbecleaned=false;
  for ( std::size_t i = 0; i < std::min((std::size_t)2,jets.size()); ++i ) {
    bool isIn=false;
    if ( jets[i].Eta()<-0.6 && jets[i].Eta()>-1.0 && 
	 jets[i].Phi()<-0.6 && jets[i].Phi()>-1.0 ) isIn=true;

    std::vector<float> sumPtTrk;
    jets[i].jet()->getAttribute(xAOD::JetAttribute::SumPtTrkPt500,sumPtTrk); // FIXME or SumPtTrkPt1000 ??
    float chf = sumPtTrk[0]/jets[i].Pt();

    float emf = 0.;
    jets[i].jet()->getAttribute(xAOD::JetAttribute::EMFrac,emf);

    if ( isIn && jets[i].Pt()>100000. &&  chf<0.20 && emf<0.3 ) shouldbecleaned=true;
  }
  return shouldbecleaned;
}


void PhysObjProxyUtils::FillNTExtraVars(NTExtraVars& extrantv,
					double MET_Track, 
					double MET_Track_phi,
					double mT2,
					double mT2_noISR,
					double gaminvRp1,
					double shatR,
					double mdeltaR,
					double cosptR,
					double gamma_R,
					double dphi_BETA_R,
					double  dphi_leg1_leg2,
					double  costhetaR,
					double dphi_BETA_Rp1_BETA_R,
					double gamma_Rp1,
					double costhetaRp1,
					double Ap)
{
  extrantv.Reset();
  extrantv.mettrack = MET_Track;
  extrantv.mettrack_phi = MET_Track_phi;
  extrantv.mT2=mT2;
  extrantv.mT2_noISR=mT2_noISR;  
  extrantv. gaminvRp1=gaminvRp1;
  extrantv. shatR=shatR;
  extrantv. mdeltaR=mdeltaR;
  extrantv. cosptR=cosptR;
  extrantv. gamma_R=gamma_R;
  extrantv. dphi_BETA_R=dphi_BETA_R;
  extrantv.  dphi_leg1_leg2=dphi_leg1_leg2;
  extrantv.  costhetaR=costhetaR;
  extrantv. dphi_BETA_Rp1_BETA_R=dphi_BETA_Rp1_BETA_R;
  extrantv. gamma_Rp1=gamma_Rp1;
  extrantv. costhetaRp1=costhetaRp1;
  extrantv. Ap=Ap;
}




void PhysObjProxyUtils::FillNTVars(NTVars& ntv, 
				   unsigned int RunNumber, 
				   unsigned int EventNumber,
				   unsigned int veto, 
				   float weight, 
				   std::vector<float>& normWeight, 
				   std::vector<float>& pileupWeight, 
				   float genWeight, 
				   float ttbarWeightHT,
				   float ttbarWeightPt2, 
				   float ttbarAvgPt,
				   float WZweight,
				   std::vector<float>& bTagWeight,
				   std::vector<float>& cTagWeight,
				   int nBJet,
				   int nCJet, 
				   double MissingEt,
				   double METPhi, 
				   double* Meff, 
				   double meffincl,
				   double minDphi, 
				   double RemainingminDPhi,
				   const std::vector<JetProxy>& good_jets,
				   int hardproc,
				   unsigned int cleaning,
				   float timing,
				   const std::vector<float>& jetSmearSystW,
				   const std::vector<float>* flaggedtau, 
				   float tauMt,
				   float SherpaBugMET,
				   float metLHTOPOx,
				   float metLHTOPOy,
				   bool isTruth)
{
  ntv.Reset();

  ntv.RunNumber = RunNumber;
  ntv.EventNumber = EventNumber;
  ntv.veto = veto;
  ntv.weight = weight;
  ntv.pileupWeight     = (pileupWeight.size() >= 1) ? pileupWeight.at(0) : 1.;
  ntv.pileupWeightUp   = (pileupWeight.size() >= 2) ? pileupWeight.at(1) : 1.;
  ntv.pileupWeightDown = (pileupWeight.size() >= 3) ? pileupWeight.at(2) : 1.;
  ntv.genWeight = genWeight;
  ntv.ttbarWeightHT = ttbarWeightHT;
  ntv.ttbarWeightPt2 = ttbarWeightPt2;
  ntv.ttbarAvgPt = ttbarAvgPt;
  ntv.WZweight = WZweight;
  ntv.bTagWeight     = (bTagWeight.size() >= 1) ? bTagWeight.at(0) : 1.;
  ntv.bTagWeightBUp   = (bTagWeight.size() >= 5) ? bTagWeight.at(4) : 1.;
  ntv.bTagWeightBDown = (bTagWeight.size() >= 2) ? bTagWeight.at(1) : 1.;
  ntv.bTagWeightCUp   = (bTagWeight.size() >= 6) ? bTagWeight.at(5) : 1.;
  ntv.bTagWeightCDown = (bTagWeight.size() >= 3) ? bTagWeight.at(2) : 1.;
  ntv.bTagWeightLUp   = (bTagWeight.size() >= 7) ? bTagWeight.at(6) : 1.;
  ntv.bTagWeightLDown = (bTagWeight.size() >= 4) ? bTagWeight.at(3) : 1.;
  ntv.cTagWeight     = (cTagWeight.size() >= 1) ? cTagWeight.at(0) : 1.;
  ntv.cTagWeightBUp   = (cTagWeight.size() >= 5) ? cTagWeight.at(4) : 1.;
  ntv.cTagWeightBDown = (cTagWeight.size() >= 2) ? cTagWeight.at(1) : 1.;
  ntv.cTagWeightCUp   = (cTagWeight.size() >= 6) ? cTagWeight.at(5) : 1.;
  ntv.cTagWeightCDown = (cTagWeight.size() >= 3) ? cTagWeight.at(2) : 1.;
  ntv.cTagWeightLUp   = (cTagWeight.size() >= 7) ? cTagWeight.at(6) : 1.;
  ntv.cTagWeightLDown = (cTagWeight.size() >= 4) ? cTagWeight.at(3) : 1.;  
  ntv.nBJet = nBJet;
  ntv.nCJet = nCJet;  
  ntv.MET = MissingEt;
  ntv.METPhi = METPhi;
  ntv.deltaPhi = minDphi;
  ntv.deltaPhiRemaining=RemainingminDPhi;
  
  ntv.MeffIncl= meffincl;
  ntv.normWeight=(normWeight.size() >= 1) ? normWeight.at(0) : 1.;
  ntv.normWeightUp=(normWeight.size() >= 2) ? normWeight.at(1) : 1.;
  ntv.normWeightDown=(normWeight.size() >= 3) ? normWeight.at(2) : 1.;

  ntv.hardproc=hardproc;
  ntv.cleaning=cleaning;
  ntv.timing=timing;

  ntv.SherpaBugMET = SherpaBugMET;

  // pdgid of incoming partons
  // FIXME
  /*
  if (!m_isData) {
    if (m_data->mcevt_pdf_id1 && m_data->mcevt_pdf_id1->size()) 
      ntv.pdfId1 = m_data->mcevt_pdf_id1->at(0);
    if (m_data->mcevt_pdf_id2 && m_data->mcevt_pdf_id2->size()) 
      ntv.pdfId2 = m_data->mcevt_pdf_id2->at(0);
  }
  */


  ntv.Njet = 0; 
  float metNOCHCORRCELLx = MissingEt*std::cos(METPhi);
  float metNOCHCORRCELLy = MissingEt*std::sin(METPhi);
  float metLHTOPONOCHCORRCELLx = metLHTOPOx;
  float metLHTOPONOCHCORRCELLy = metLHTOPOy;
  TLorentzVector jet1TLV;
  TLorentzVector jet2TLV;
  for ( size_t jet0=0; jet0<good_jets.size(); jet0++) 
  {
    const JetProxy& thisjet = good_jets[jet0];
    float pt = thisjet.Pt();
    float eta = (pt > 0.) ? thisjet.Eta() : 0.;
    float phi = (pt > 0.) ? thisjet.Phi() : 0.;
    float m = (pt > 0.) ? thisjet.M() : 0.;

    float emf = 0.;
    float chf = 1.;
    float jetbtag = -1.;  // MV1 not available on xAOD
    if ( thisjet.jet() && !isTruth ) {
      thisjet.jet()->getAttribute(xAOD::JetAttribute::EMFrac,emf);
      std::vector<float> sumPtTrk;
      thisjet.jet()->getAttribute(xAOD::JetAttribute::SumPtTrkPt500,sumPtTrk); // FIXME or SumPtTrkPt1000 ??
      chf = (pt > 0.) ? sumPtTrk[0]/pt : 0;
    }
    int jetflav = 0;
    float jettagU = -1.;
    float jettagB = -1.;
    float jettagC = -1.;
    double  jet_BCH_CORR_CELL = 0.;
    float jetFracSamplingMax = 0.;
    int jetFracSamplingMaxIndex = 0.;
    if ( thisjet.jet() ) {
      thisjet.jet()->getAttribute(xAOD::JetAttribute::BchCorrCell,jet_BCH_CORR_CELL);
      thisjet.jet()->getAttribute(xAOD::JetAttribute::FracSamplingMax,jetFracSamplingMax);
      thisjet.jet()->getAttribute(xAOD::JetAttribute::FracSamplingMaxIndex,jetFracSamplingMaxIndex);
      metNOCHCORRCELLx += thisjet.jet()->px()*jet_BCH_CORR_CELL;
      metNOCHCORRCELLy += thisjet.jet()->py()*jet_BCH_CORR_CELL;
      metLHTOPONOCHCORRCELLx += thisjet.jet()->px()*jet_BCH_CORR_CELL;
      metLHTOPONOCHCORRCELLy += thisjet.jet()->py()*jet_BCH_CORR_CELL;
    }

    if (pt > 40000.) {
      ntv.jetPt.push_back(pt);
      ntv.jetEta.push_back(eta);
      ntv.jetPhi.push_back(phi);
      ntv.jetM.push_back(m);
      ntv.jetBTag.push_back(jetbtag);
      ntv.jetFlav.push_back(jetflav);
      ntv.jetTagU.push_back(jettagU);
      ntv.jetTagB.push_back(jettagB);
      ntv.jetTagC.push_back(jettagC);
      ntv.jetBCH_CORR_CELL.push_back(jet_BCH_CORR_CELL);
      ntv.jetFracSamplingMax.push_back(jetFracSamplingMax);
      ntv.jetFracSamplingMaxIndex.push_back(jetFracSamplingMaxIndex);
      ntv.Njet++;
    }
    switch (jet0)
    {
    case 0:
      jet1TLV.SetPtEtaPhiM(pt,eta,phi,m);
      ntv.emfjet0 = emf;
      ntv.chfjet0 = chf;
      break;
    case 1:
      jet2TLV.SetPtEtaPhiM(pt,eta,phi,m);
      ntv.emfjet1 = emf;
      ntv.chfjet1 = chf;
      break;
    }
  }
  ntv.metNOCHCORRCELL = std::sqrt(metNOCHCORRCELLx*metNOCHCORRCELLx+
				  metNOCHCORRCELLy*metNOCHCORRCELLy);
  ntv.metLHTOPONOCHCORRCELL = std::sqrt(metLHTOPONOCHCORRCELLx*metLHTOPONOCHCORRCELLx+
					metLHTOPONOCHCORRCELLy*metLHTOPONOCHCORRCELLy);
  ntv.metLHTOPO = std::sqrt(metLHTOPOx*metLHTOPOx+
			    metLHTOPOy*metLHTOPOy);
  
  ntv.jetSmearSystW = jetSmearSystW;
}


void PhysObjProxyUtils::FillNTReclusteringVars(NTReclusteringVars& RTntv, 
				   const std::vector<JetProxy>& good_jets)
{
  RTntv.Reset(); 
  
  //Reclustering:
  PhysObjProxyUtils::ReclJets myRT;
  const float fCut=0.1;  //trimming cut 
  myRT=Recluster(good_jets, 40000., fCut, 1.0);
  RTntv.RTjets10SubJetIndeces = myRT.recl_jets_subInds;
  RTntv.RTjetM = myRT.recl_jets_M;
  RTntv.RTjetPt = myRT.recl_jets_Pt;
  RTntv.RTjetEta = myRT.recl_jets_Eta;
  RTntv.RTjetPhi = myRT.recl_jets_Phi;

  int NWcandidates= 0;
  for (unsigned int iRT=0; iRT< RTntv.RTjetM.size(); ++iRT){
    if(RTntv.RTjetM[iRT]>60000. && RTntv.RTjetM[iRT]<100000.) NWcandidates++; 
  }
  RTntv.NWcandidates= NWcandidates;
  
}

PhysObjProxyUtils::ReclJets PhysObjProxyUtils::Recluster(const std::vector<JetProxy>& small_jets, double PTcut, double fcut, double jetRad){
  vector<int> NumOfSubJets;
  vector<fastjet::PseudoJet> particles = ObjsToPJ(small_jets);
  fastjet::JetDefinition fJetDef(fastjet::antikt_algorithm, jetRad,fastjet::E_scheme, fastjet::Best);
  fastjet::ClusterSequence hardClustSeq(particles, fJetDef);
  vector<fastjet::PseudoJet> StandardJets = fastjet::sorted_by_pt(hardClustSeq.inclusive_jets(PTcut));
  vector<TLorentzVector> RCjets;
  for (unsigned int i=0; i<StandardJets.size(); ++i){
    TLorentzVector sub = TLorentzVector();   
    sub.SetPtEtaPhiE(StandardJets[i].pt(), StandardJets[i].eta(), StandardJets[i].phi(), StandardJets[i].e());
    vector<fastjet::PseudoJet> constituents = StandardJets[i].constituents();
    for(unsigned int iCons = 0; iCons < constituents.size(); ++iCons){
      //Do something with the small radius jets.
    }
    RCjets.push_back(sub);
  }   
  
  vector<TLorentzVector> RTjets; 
  RTjets.clear();
  vector<vector<int> > RTjets_small_jets_inds;
  vector<float> RTjets_M;

  //Now for my trimming on the re-clustered jets 
  for (unsigned int i=0; i<StandardJets.size(); ++i){
    int NumSubJets=0; 
    TLorentzVector trimmedjet = TLorentzVector(); 
    vector<fastjet::PseudoJet> constituents = StandardJets[i].constituents();
    //float emfrac_recalculated=0;
    vector<int> this_jet_subinds;
    for(unsigned int iCons = 0; iCons < constituents.size(); ++iCons){
      TLorentzVector subjet = TLorentzVector();
      subjet.SetPtEtaPhiE(constituents[iCons].pt(), constituents[iCons].eta(), constituents[iCons].phi(), constituents[iCons].e());    
      //float emfrac= FindEmFrac(small_jets, subjet.Pt());
      if (subjet.Pt() > fcut*RCjets[i].Pt()){     
        //emfrac_recalculated+= subjet.E()*emfrac;
	trimmedjet+=subjet;
        NumSubJets+=1;
        for (unsigned int iSub = 0; iSub < small_jets.size(); ++iSub){
            if (fabs(small_jets[iSub].Pt()-subjet.Pt())>0.00001) continue;
            if (fabs(small_jets[iSub].Phi()-subjet.Phi())>0.01) continue;
            if (fabs(small_jets[iSub].Eta()-subjet.Eta())>0.01) continue;            
            this_jet_subinds.push_back(iSub);      
        }
      }
    }
    //emfrac_recalculated= emfrac_recalculated/trimmedjet.E();
    //if (emfrac_recalculated<0.99) {
      NumOfSubJets.push_back(NumSubJets);
      RTjets.push_back(trimmedjet);
      RTjets_small_jets_inds.push_back(this_jet_subinds);
    //}  
  }
  
  //finally sort RTjets according to pt
  vector<int> sorted_indexes= GetSortedJetIndexes(RTjets);
  PhysObjProxyUtils::ReclJets myRecl;
  for (unsigned int i=0; i<sorted_indexes.size(); ++i){
    myRecl.recl_jets_tlv.push_back(RTjets[sorted_indexes[i]]);
    myRecl.recl_jets_subInds.push_back(RTjets_small_jets_inds[sorted_indexes[i]]);
    myRecl.recl_jets_Pt.push_back(RTjets[sorted_indexes[i]].Pt());
    myRecl.recl_jets_Eta.push_back(RTjets[sorted_indexes[i]].Eta());
    myRecl.recl_jets_Phi.push_back(RTjets[sorted_indexes[i]].Phi());
    myRecl.recl_jets_M.push_back(RTjets[sorted_indexes[i]].M());
    myRecl.sub_jets.push_back(NumOfSubJets[sorted_indexes[i]]);
    if(i==0){
      //fill leading jet constituents
      vector<fastjet::PseudoJet> constituents = StandardJets[sorted_indexes[i]].constituents();
      for(unsigned int iCons = 0; iCons < constituents.size(); ++iCons){
        TLorentzVector subjet = TLorentzVector();
        subjet.SetPtEtaPhiE(constituents[iCons].pt(), constituents[iCons].eta(), constituents[iCons].phi(), constituents[iCons].e());    
	if (subjet.Pt() > fcut*StandardJets[sorted_indexes[i]].pt()){     
	  myRecl.recljet1_smalljets_tlv.push_back(subjet);
	}
      }	
      //extract two latest protojets
      fastjet::PseudoJet parent1, parent2;
      if (!hardClustSeq.has_parents(StandardJets[sorted_indexes[i]],parent1,parent2)) {
          myRecl.recljet1_massdrop= 0;
	 continue;
      }
      vector<fastjet::PseudoJet> constituents1 = parent1.constituents();
      TLorentzVector parentTrim1(0,0,0,0),parentTrim2(0,0,0,0);
      for(unsigned int iCons = 0; iCons < constituents1.size(); ++iCons){
        TLorentzVector subjet = TLorentzVector();
        subjet.SetPtEtaPhiE(constituents1[iCons].pt(), constituents1[iCons].eta(),
	constituents1[iCons].phi(), constituents1[iCons].e());    
	if (subjet.Pt() > fcut*StandardJets[sorted_indexes[i]].pt()){     
	  parentTrim1+=subjet;
	}
      }
      vector<fastjet::PseudoJet> constituents2 = parent2.constituents();
      for(unsigned int iCons = 0; iCons < constituents2.size(); ++iCons){
        TLorentzVector subjet = TLorentzVector();
        subjet.SetPtEtaPhiE(constituents2[iCons].pt(), constituents2[iCons].eta(),
	constituents2[iCons].phi(), constituents2[iCons].e());    
	if (subjet.Pt() > fcut*StandardJets[sorted_indexes[i]].pt()){     
	  parentTrim2+=subjet;
	}
      }
      myRecl.recljet1_massdrop= TMath::Max(parentTrim1.M(),parentTrim2.M())/RTjets[sorted_indexes[i]].M();
    }
  }
  
  return myRecl;
}   

std::vector<int> PhysObjProxyUtils::GetSortedJetIndexes(const std::vector<TLorentzVector> jets)
{
   float vec_pt[500];
   int vec_index[500];
  
   for (UInt_t iJet=0; iJet< jets.size(); ++iJet){
       vec_pt[iJet]= jets[iJet].Pt();
       vec_index[iJet]= iJet;
   }
  
   // now obtain list of sorted indexes
   int sorted_index[500];
   TMath::Sort((int)jets.size(),vec_pt,sorted_index);
 
   vector<int> sorted;
   if(!sorted.empty()) sorted.clear();
   for (UInt_t iJet=0; iJet< jets.size(); ++iJet){
       sorted.push_back(sorted_index[iJet]);
   }
   return sorted;
  
}


std::vector<fastjet::PseudoJet> PhysObjProxyUtils::ObjsToPJ(const std::vector<JetProxy>& good_jets){
  vector<fastjet::PseudoJet> out;
  for(unsigned int i = 0; i < good_jets.size(); i++){
      const JetProxy& thisjet = good_jets[i];
      fastjet::PseudoJet newjet(thisjet.Px(), thisjet.Py(), thisjet.Pz(), thisjet.E());                                                                                     
      out.push_back(newjet);
  }
  return out;
}
