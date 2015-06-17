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

void PhysObjProxyUtils::RJigsawInit(){

  LAB_alt = new RestFrames::RLabFrame("LAB","lab");

  S_alt = new RestFrames::RSelfAssemblingFrame("CM","CM");
  V_alt = new RestFrames::RVisibleFrame("V_alt","Vis");
  I_alt = new RestFrames::RInvisibleFrame("I_alt","Iinv");
  INV_alt = new RestFrames::InvisibleGroup ("INV_alt","Invisible State Jigsaws");
  VIS_alt = new RestFrames::CombinatoricGroup("VIS_alt","Visible Object Jigsaws");


  MinMass_alt = new RestFrames::InvisibleMassJigsaw("MINMASS_JIGSAW_ALT", "Invisible system mass Jigsaw");
  Rapidity_alt = new RestFrames::InvisibleRapidityJigsaw("RAPIDITY_JIGSAW_ALT", "Invisible system rapidity Jigsaw");


  LAB = new RestFrames::RLabFrame("LAB","lab");
  SS = new RestFrames::RDecayFrame("SS","SS");
  S1 = new RestFrames::RSelfAssemblingFrame("S1","#tilde{S}_{a}");
  S2 = new RestFrames::RSelfAssemblingFrame("S2","#tilde{S}_{b}");
  V1 = new RestFrames::RVisibleFrame("V1","V_{a}");
  V2 = new RestFrames::RVisibleFrame("V2","V_{b}");
  I1 = new RestFrames::RInvisibleFrame("I1","I_{a}");
  I2 = new RestFrames::RInvisibleFrame("I2","I_{b}");
  INV = new RestFrames::InvisibleGroup("INV","Invisible State Jigsaws");
  VIS = new RestFrames::CombinatoricGroup("VIS","Visible Object Jigsaws");

  MinMassJigsaw = new RestFrames::InvisibleMassJigsaw("MINMASS_JIGSAW", "Invisible system mass Jigsaw");
  RapidityJigsaw = new RestFrames::InvisibleRapidityJigsaw("RAPIDITY_JIGSAW", "Invisible system rapidity Jigsaw");
  ContraBoostJigsaw = new RestFrames::ContraBoostInvariantJigsaw("CB_JIGSAW","Contraboost invariant Jigsaw");
  HemiJigsaw = new RestFrames::MinimizeMassesCombinatoricJigsaw("HEM_JIGSAW","Minimize m _{V_{a,b}} Jigsaw");

  LAB_R = new RestFrames::RLabFrame("LAB_R","LAB");
  GG_R = new RestFrames::RDecayFrame("GG_R","#tilde{g}#tilde{g}");
  Ga_R = new RestFrames::RDecayFrame("Ga_R","#tilde{g}_{a}");
  Gb_R = new RestFrames::RDecayFrame("Gb_R","#tilde{g}_{b}");
  Ca_R = new RestFrames::RDecayFrame("Ca_R","C_{a}");
  Cb_R = new RestFrames::RDecayFrame("Cb_R","C_{b}");
  V1a_R = new RestFrames::RVisibleFrame("V1a_R","j_{1a}");
  V2a_R = new RestFrames::RVisibleFrame("V2a_R","j_{2a}");
  Xa_R = new RestFrames::RInvisibleFrame("Xa_R","#tilde{#chi}_{a}");
  V1b_R = new RestFrames::RVisibleFrame("V1b_R","j_{1b}");
  V2b_R = new RestFrames::RVisibleFrame("V2b_R","j_{2b}");
  Xb_R = new RestFrames::RInvisibleFrame("Xb_R","#tilde{#chi}_{b}");
  INV_R = new RestFrames::InvisibleGroup ("INV_R","WIMP Jigsaws");
  VIS_R = new RestFrames::CombinatoricGroup("VIS","Visible Object Jigsaws");
  MinMassJigsaw_R = new RestFrames::InvisibleMassJigsaw("MINMASS_R", "Invisible system mass Jigsaw");
  RapidityJigsaw_R = new RestFrames::InvisibleRapidityJigsaw("RAPIDITY_R", "Invisible system rapidity Jigsaw");
  ContraBoostJigsaw_R = new RestFrames::ContraBoostInvariantJigsaw("CONTRA_R","Contraboost invariant Jigsaw");
  HemiJigsaw_R = new RestFrames::MinimizeMassesCombinatoricJigsaw ("HEM_JIGSAW_R","Minimize m _{V_{a,b}} Jigsaw");
  CaHemiJigsaw_R = new RestFrames::MinimizeMassesCombinatoricJigsaw("CbHEM_JIGSAW_R","Minimize m _{C_{a}} Jigsaw");
  CbHemiJigsaw_R = new RestFrames::MinimizeMassesCombinatoricJigsaw("CaHEM_JIGSAW_R","Minimize m _{C_{b}} Jigsaw");


  INV_alt->AddFrame(I_alt);
  VIS_alt->AddFrame(V_alt);
  VIS_alt->SetNElementsForFrame(V_alt,1,false);

  LAB_alt->SetChildFrame(S_alt);
  S_alt->AddChildFrame(V_alt);
  S_alt->AddChildFrame(I_alt);

  LAB_alt->InitializeTree(); 

// Will just set invisible mass to zero
  INV_alt->AddJigsaw(MinMass_alt);

// will set rapidity to zero
  INV_alt->AddJigsaw(Rapidity_alt);
  Rapidity_alt->AddVisibleFrame((LAB_alt->GetListVisibleFrames()));

  LAB_alt->InitializeAnalysis(); 

  //
  //
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  // SQUARK TREE /////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  //
  //

  // The invisible group is all of the weakly interacting particles
  INV->AddFrame(I1);
  INV->AddFrame(I2);

  // the combinatoric group is the list of visible objects
  // that go into our hemispheres 
  VIS->AddFrame(V1);
  VIS->SetNElementsForFrame(V1,1,false);
  VIS->AddFrame(V2);
  VIS->SetNElementsForFrame(V2,1,false);

  LAB->SetChildFrame(SS);

  SS->AddChildFrame(S1);
  SS->AddChildFrame(S2);

  S1->AddChildFrame(V1);
  S1->AddChildFrame(I1);
  S2->AddChildFrame(V2);
  S2->AddChildFrame(I2);

  std::cout << "Is consistent tree topology? " << LAB->InitializeTree() << std::endl; 

  INV->AddJigsaw(MinMassJigsaw);

  INV->AddJigsaw(RapidityJigsaw);
  RapidityJigsaw->AddVisibleFrame((LAB->GetListVisibleFrames()));

  INV->AddJigsaw(ContraBoostJigsaw);
  ContraBoostJigsaw->AddVisibleFrame((S1->GetListVisibleFrames()), 0);
  ContraBoostJigsaw->AddVisibleFrame((S2->GetListVisibleFrames()), 1);
  ContraBoostJigsaw->AddInvisibleFrame((S1->GetListInvisibleFrames()), 0);
  ContraBoostJigsaw->AddInvisibleFrame((S2->GetListInvisibleFrames()), 1);

  VIS->AddJigsaw(HemiJigsaw);
  HemiJigsaw->AddFrame(V1,0);
  HemiJigsaw->AddFrame(V2,1);

  std::cout << "Is consistent analysis tree? : " << LAB->InitializeAnalysis() << std::endl; 

  //
  //
  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  // GLUINO TREE /////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  //
  //


  // Set up 'signal-like' analysis tree
  LAB_R->SetChildFrame(GG_R);
  GG_R->AddChildFrame(Ga_R);
  GG_R->AddChildFrame(Gb_R);
  Ga_R->AddChildFrame(V1a_R);
  Ga_R->AddChildFrame(Ca_R);
  Ca_R->AddChildFrame(V2a_R);
  Ca_R->AddChildFrame(Xa_R);
  Gb_R->AddChildFrame(V1b_R);
  Gb_R->AddChildFrame(Cb_R);
  Cb_R->AddChildFrame(V2b_R);
  Cb_R->AddChildFrame(Xb_R);


  if(!LAB_R->InitializeTree()) cout << "Problem with signal-like reconstruction tree" << endl; 


  INV_R->AddFrame(Xa_R);
  INV_R->AddFrame(Xb_R);
  // visible frames in first decay step must always have at least one element
  VIS_R->AddFrame(V1a_R);
  VIS_R->AddFrame(V1b_R);
  VIS_R->SetNElementsForFrame(V1a_R,1,false);
  VIS_R->SetNElementsForFrame(V1b_R,1,false);
  // visible frames in second decay step can have zero elements
  VIS_R->AddFrame(V2a_R);
  VIS_R->AddFrame(V2b_R);
  VIS_R->SetNElementsForFrame(V2a_R,0,false);
  VIS_R->SetNElementsForFrame(V2b_R,0,false);

  INV_R->AddJigsaw(MinMassJigsaw_R);
  INV_R->AddJigsaw(RapidityJigsaw_R);
  RapidityJigsaw_R->AddVisibleFrame((LAB_R->GetListVisibleFrames()));
  INV_R->AddJigsaw(ContraBoostJigsaw_R);
  ContraBoostJigsaw_R->AddVisibleFrame((Ga_R->GetListVisibleFrames()), 0);
  ContraBoostJigsaw_R->AddVisibleFrame((Gb_R->GetListVisibleFrames()), 1);
  ContraBoostJigsaw_R->AddInvisibleFrame((Ga_R->GetListInvisibleFrames()), 0);
  ContraBoostJigsaw_R->AddInvisibleFrame((Gb_R->GetListInvisibleFrames()), 1);
  VIS_R->AddJigsaw(HemiJigsaw_R);
  HemiJigsaw_R->AddFrame(V1a_R,0);
  HemiJigsaw_R->AddFrame(V1b_R,1);
  HemiJigsaw_R->AddFrame(V2a_R,0);
  HemiJigsaw_R->AddFrame(V2b_R,1);
  VIS_R->AddJigsaw(CaHemiJigsaw_R);
  CaHemiJigsaw_R->AddFrame(V1a_R,0);
  CaHemiJigsaw_R->AddFrame(V2a_R,1);
  CaHemiJigsaw_R->AddFrame(Xa_R,1);
  VIS_R->AddJigsaw(CbHemiJigsaw_R);
  CbHemiJigsaw_R->AddFrame(V1b_R,0);
  CbHemiJigsaw_R->AddFrame(V2b_R,1);
  CbHemiJigsaw_R->AddFrame(Xb_R,1);

  if(!LAB_R->InitializeAnalysis()) cout << "Problem with signal-tree jigsaws" << endl;

  return;

}


void PhysObjProxyUtils::CalculateRJigsawVariables(const std::vector<JetProxy>& jets, 
  Double_t metx,
  Double_t mety,
  std::map<TString,float>& RJigsawVariables){



  LAB->ClearEvent();
  LAB_R->ClearEvent();
  LAB_alt->ClearEvent();


  vector<RestFrames::GroupElementID> jetID_R;                    // ID for tracking jets in tree

  std::cout << "number of jets is " << jets.size() << std::endl;

  if(jets.size() < 2){
    RJigsawVariables = std::map<TString, float>();
    return;
  } 

  // Still need to add jets to frames ///////////////
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

  for(size_t ijet=0; ijet<jets.size(); ijet++) 
    {
      VIS->AddLabFrameFourVector( myjets[ijet] );
      jetID_R.push_back( VIS_R->AddLabFrameFourVector( myjets[ijet] )  );
      TLorentzVector temp = myjets[ijet];
      temp.SetPtEtaPhiM(temp.Pt(),0.,temp.Phi(),temp.M());
      VIS_alt->AddLabFrameFourVector( temp );
    }



  TVector3 MET_TV3;

  MET_TV3.SetZ(0.);
  MET_TV3.SetX(metx);
  MET_TV3.SetY(mety);


  INV->SetLabFrameThreeVector(MET_TV3);
  LAB->AnalyzeEvent();

  INV_alt->SetLabFrameThreeVector(MET_TV3);
  LAB_alt->AnalyzeEvent();


  // if(goodJets->size()>3){
  INV_R->SetLabFrameThreeVector(MET_TV3);
  LAB_R->AnalyzeEvent();



  ////////////////////////////////////////////////////////////////////////////////
  // 1st order squark vars


  RJigsawVariables[ "RJVars_SS_Mass"           ] = SS->GetMass();
  RJigsawVariables[ "RJVars_SS_InvGamma"       ] = 1./SS->GetGammaInParentFrame();
  RJigsawVariables[ "RJVars_SS_dPhiBetaR"      ] = SS->GetDeltaPhiBoostVisible();
  RJigsawVariables[ "RJVars_SS_dPhiVis"        ] = SS->GetDeltaPhiVisible();
  RJigsawVariables[ "RJVars_SS_CosTheta"       ] = SS->GetCosDecayAngle();
  RJigsawVariables[ "RJVars_SS_dPhiDecayAngle" ] = SS->GetDeltaPhiDecayAngle();
  RJigsawVariables[ "RJVars_SS_VisShape"       ] = SS->GetVisibleShape();
  RJigsawVariables[ "RJVars_SS_MDeltaR"        ] = SS->GetVisibleShape() * SS->GetMass() ;
  RJigsawVariables[ "RJVars_S1_Mass"           ] = S1->GetMass();
  RJigsawVariables[ "RJVars_S1_CosTheta"       ] = S1->GetCosDecayAngle();
  RJigsawVariables[ "RJVars_S2_Mass"           ] = S2->GetMass();
  RJigsawVariables[ "RJVars_S2_CosTheta"       ] = S2->GetCosDecayAngle();
  RJigsawVariables[ "RJVars_I1_Depth"          ] = S1->GetFrameDepth(I1);
  RJigsawVariables[ "RJVars_I2_Depth"          ] = S2->GetFrameDepth(I2);
  RJigsawVariables[ "RJVars_V1_N"              ] = VIS->GetNElementsInFrame(V1);
  RJigsawVariables[ "RJVars_V2_N"              ] = VIS->GetNElementsInFrame(V2);

  // end
  ////////////////////////////////////////////////////////////////////////////////





  ////////////////////////////////////////////////////////////////////////////////
  // 2nd order "gluino-like" vars

  RestFrames::RDecayFrame* G[2];
  RestFrames::RDecayFrame* C[2];
  RestFrames::RVisibleFrame* VS[2];
  RestFrames::RVisibleFrame* VC[2];
  RestFrames::RInvisibleFrame* X[2];
  // Randomize the two hemispheres
  int flip = (m_random.Rndm() > 0.5);
  G[flip] = Ga_R;
  G[(flip+1)%2] = Gb_R;
  C[flip] = Ca_R;
  C[(flip+1)%2] = Cb_R;
  VS[flip] = V1a_R;
  VS[(flip+1)%2] = V1b_R;
  VC[flip] = V2a_R;
  VC[(flip+1)%2] = V2b_R;
  X[flip] = Xa_R;
  X[(flip+1)%2] = Xb_R;


  double NV[2];
  double jet1PT[2];
  double jet2PT[2];


  for(int i = 0; i < 2; i++){

    NV[i] =  VIS_R->GetNElementsInFrame(VS[i]);
    NV[i] += VIS_R->GetNElementsInFrame(VC[i]);

    int N = jetID_R.size();
    // std::cout << "In SklimmerAnalysis:  N Jets " << N << std::endl;

    double pTmax[2]; pTmax[0] = -1.; pTmax[1] = -1.;
    for(int j = 0; j < N; j++){
      const RestFrames::RestFrame* frame = VIS_R->GetFrame(jetID_R[j]);
      if(VS[i]->IsSame(frame) || VC[i]->IsSame(frame)){
        double pT = VIS_R->GetLabFrameFourVector(jetID_R[j]).Pt();
        // std::cout << "In SklimmerAnalysis: ijet pT " << pT << std::endl;

        if(pT > pTmax[0]){
          pTmax[1] = pTmax[0];
          pTmax[0] = pT;
        } else {
          if(pT > pTmax[1]) pTmax[1] = pT;
        }
      }
    }

    jet1PT[i] = pTmax[0];
    jet2PT[i] = pTmax[1];
    std::cout << "In SklimmerAnalysis: " << jet1PT[i] << " " << jet2PT[i] << std::endl;
    std::cout << "In SklimmerAnalysis: " << NV[i] << std::endl;



    if(NV[i] > 1){
      RJigsawVariables[Form("RJVars_C_%d_CosTheta",i)     ] = C[i]->GetCosDecayAngle();
      RJigsawVariables[Form("RJVars_G_%d_dPhiGC",i)       ] = G[i]->GetDeltaPhiDecayPlanes(C[i]);
      RJigsawVariables[Form("RJVars_G_%d_MassRatioGC",i)  ] = (C[i]->GetMass()-X[i]->GetMass())/(G[i]->GetMass()-X[i]->GetMass());
    } else {
      RJigsawVariables[Form("RJVars_C_%d_CosTheta",i)     ] = -10.;
      RJigsawVariables[Form("RJVars_G_%d_dPhiGC",i)       ] = -10.;
      RJigsawVariables[Form("RJVars_G_%d_MassRatioGC",i)  ] = -10.;
    }

    RJigsawVariables[ Form("RJVars_G_%d_CosTheta",i)    ] = G[i]->GetCosDecayAngle();
    RJigsawVariables[ Form("RJVars_G_%d_Jet1_pT",i)     ] = jet1PT[i];
    RJigsawVariables[ Form("RJVars_G_%d_Jet2_pT",i)     ] = jet2PT[i];

  }



  // Some new Gluino variables....

  TLorentzVector vV1 = G[0]->GetVisibleFourVector(G[0]);
  TLorentzVector vV2 = G[1]->GetVisibleFourVector(G[1]);
  float MG = (vV1.M2()-vV2.M2())/(2.*(vV1.E()-vV2.E()));

  float PG = G[0]->GetMomentum(GG_R);
  float MGG = 2.*sqrt(PG*PG + MG*MG);
  float gaminvGG = 2.*MG/MGG;
  //float gaminv = 1./SS->GetGammaInParentFrame();
  float gaminv = SS->GetVisibleShape();
  float beta = sqrt(1.- gaminv*gaminv);
  float betaGG = sqrt(1.- gaminvGG*gaminvGG);

  //*** velocity difference between 'massive' and 'mass-less'
  float DeltaBetaGG = -(betaGG-beta)/(1.-betaGG*beta);

  //*** delta phi between GG visible decay products and GG decay axis
  float dphiVG = GG_R->GetDeltaPhiDecayVisible();


  RJigsawVariables[ "RJVars_MG"          ] = MG;
  RJigsawVariables[ "RJVars_DeltaBetaGG" ] = DeltaBetaGG;
  RJigsawVariables[ "RJVars_dphiVG"      ] = dphiVG;


  // Gluino-level variables end
  ////////////////////////////////////////////////////////////////////////////////




  ////////////////////////////////////////////////////////////////////////////////
  // QCD Variables


  // dphiR and Rptshat (formerly cosPT)
  // for QCD rejection
  double dphiR = SS->GetDeltaPhiBoostVisible();
  double PTCM = SS->GetFourVector(LAB).Pt();
  double Rptshat = PTCM / (PTCM + SS->GetMass()/4.);

  // QCD rejection using the 'background tree'
  // MET 'sibling' in background tree auxillary calculations
  TLorentzVector Psib = I_alt->GetSiblingFrame()->GetFourVector(LAB_alt);
  TLorentzVector Pmet = I_alt->GetFourVector(LAB_alt);
  double Psib_dot_METhat = max(0., Psib.Vect().Dot(MET_TV3.Unit()));
  double Mpar2 = Psib.E()*Psib.E()-Psib.Vect().Dot(MET_TV3.Unit())*Psib.Vect().Dot(MET_TV3.Unit());
  double Msib2 = Psib.M2();
  double MB2 = 2.*(Pmet.E()*Psib.E()-MET_TV3.Dot(Psib.Vect()));
  TVector3 boostPsibM = (Pmet+Psib).BoostVector();


  // QCD rejection variables from 'background tree'
  double DepthBKG = S_alt->GetFrameDepth(I_alt);
  int Nsib = I_alt->GetSiblingFrame()->GetNDescendants();
  double cosBKG = I_alt->GetParentFrame()->GetCosDecayAngle();
  double dphiMsib = fabs(MET_TV3.DeltaPhi(Psib.Vect()));
  double RpsibM = Psib_dot_METhat / (Psib_dot_METhat + MET_TV3.Mag());
  double RmsibM = 1. / ( MB2/(Mpar2-Msib2) + 1.);
  Psib.Boost(-boostPsibM);
  double cosPsibM = -1.*Psib.Vect().Unit().Dot(boostPsibM.Unit());
  cosPsibM = (1.-cosPsibM)/2.;
  double DeltaQCD1 = (cosPsibM-RpsibM)/(cosPsibM+RpsibM);
  double DeltaQCD2 = (cosPsibM-RmsibM)/(cosPsibM+RmsibM);

  RJigsawVariables[ "RJVars_QCD_dPhiR"    ] = dphiR;
  RJigsawVariables[ "RJVars_QCD_Rpt"      ] = Rptshat;
  RJigsawVariables[ "RJVars_QCD_Rmsib"    ] = RmsibM;
  RJigsawVariables[ "RJVars_QCD_Delta2"   ]  = DeltaQCD2;
  RJigsawVariables[ "RJVars_QCD_Rpsib"    ] = RpsibM;
  RJigsawVariables[ "RJVars_QCD_Delta1"   ]  = DeltaQCD1;

  // end
  ////////////////////////////////////////////////////////////////////////////////



  return;

}


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
  extrantv.mettrack = MET_Track * 0.001;
  extrantv.mettrack_phi = MET_Track_phi;
  extrantv.mT2=mT2 * 0.001;
  extrantv.mT2_noISR=mT2_noISR * 0.001;  
  extrantv.gaminvRp1=gaminvRp1;
  extrantv.shatR=shatR;
  extrantv.mdeltaR=mdeltaR;
  extrantv.cosptR=cosptR;
  extrantv.gamma_R=gamma_R;
  extrantv.dphi_BETA_R=dphi_BETA_R;
  extrantv.dphi_leg1_leg2=dphi_leg1_leg2;
  extrantv.costhetaR=costhetaR;
  extrantv.dphi_BETA_Rp1_BETA_R=dphi_BETA_Rp1_BETA_R;
  extrantv.gamma_Rp1=gamma_Rp1;
  extrantv.costhetaRp1=costhetaRp1;
  extrantv.Ap=Ap;
}

void PhysObjProxyUtils::FillNTRJigsawVars(NTRJigsawVars& rjigsawntv,
              std::map<TString,float> & RJigsawVariables
          )
{
  // rjigsawntv.Reset();

  std::cout << "In filling function----------------" << std::endl;
  std::cout << RJigsawVariables["RJVars_G_0_CosTheta"] << " -----------------------" << std::endl;

  rjigsawntv.RJVars_SS_Mass           = RJigsawVariables["RJVars_SS_Mass"] ;
  rjigsawntv.RJVars_SS_InvGamma       = RJigsawVariables["RJVars_SS_InvGamma"] ;
  rjigsawntv.RJVars_SS_dPhiBetaR      = RJigsawVariables["RJVars_SS_dPhiBetaR"] ;
  rjigsawntv.RJVars_SS_dPhiVis        = RJigsawVariables["RJVars_SS_dPhiVis"] ;
  rjigsawntv.RJVars_SS_CosTheta       = RJigsawVariables["RJVars_SS_CosTheta"] ;
  rjigsawntv.RJVars_SS_dPhiDecayAngle = RJigsawVariables["RJVars_SS_dPhiDecayAngle"] ;
  rjigsawntv.RJVars_SS_VisShape       = RJigsawVariables["RJVars_SS_VisShape"] ;
  rjigsawntv.RJVars_SS_MDeltaR        = RJigsawVariables["RJVars_SS_MDeltaR"] ;
  rjigsawntv.RJVars_S1_Mass           = RJigsawVariables["RJVars_S1_Mass"] ;
  rjigsawntv.RJVars_S1_CosTheta       = RJigsawVariables["RJVars_S1_CosTheta"] ;
  rjigsawntv.RJVars_S2_Mass           = RJigsawVariables["RJVars_S2_Mass"] ;
  rjigsawntv.RJVars_S2_CosTheta       = RJigsawVariables["RJVars_S2_CosTheta"] ;
  rjigsawntv.RJVars_I1_Depth          = RJigsawVariables["RJVars_I1_Depth"] ;
  rjigsawntv.RJVars_I2_Depth          = RJigsawVariables["RJVars_I2_Depth"] ;
  rjigsawntv.RJVars_V1_N              = RJigsawVariables["RJVars_V1_N"] ;
  rjigsawntv.RJVars_V2_N              = RJigsawVariables["RJVars_V2_N"] ;
  rjigsawntv.RJVars_MG                = RJigsawVariables["RJVars_MG"] ;
  rjigsawntv.RJVars_DeltaBetaGG       = RJigsawVariables["RJVars_DeltaBetaGG"] ;
  rjigsawntv.RJVars_dphiVG            = RJigsawVariables["RJVars_dphiVG"] ;
  rjigsawntv.RJVars_G_0_CosTheta      = RJigsawVariables["RJVars_G_0_CosTheta"] ;
  rjigsawntv.RJVars_C_0_CosTheta      = RJigsawVariables["RJVars_C_0_CosTheta"] ;
  rjigsawntv.RJVars_G_0_dPhiGC        = RJigsawVariables["RJVars_G_0_dPhiGC"] ;
  rjigsawntv.RJVars_G_0_MassRatioGC   = RJigsawVariables["RJVars_G_0_MassRatioGC"] ;
  rjigsawntv.RJVars_G_0_Jet1_pT       = RJigsawVariables["RJVars_G_0_Jet1_pT"] ;
  rjigsawntv.RJVars_G_0_Jet2_pT       = RJigsawVariables["RJVars_G_0_Jet2_pT"] ;
  rjigsawntv.RJVars_G_1_CosTheta      = RJigsawVariables["RJVars_G_1_CosTheta"] ;
  rjigsawntv.RJVars_C_1_CosTheta      = RJigsawVariables["RJVars_C_1_CosTheta"] ;
  rjigsawntv.RJVars_G_1_dPhiGC        = RJigsawVariables["RJVars_G_1_dPhiGC"] ;
  rjigsawntv.RJVars_G_1_MassRatioGC   = RJigsawVariables["RJVars_G_1_MassRatioGC"] ;
  rjigsawntv.RJVars_G_1_Jet1_pT       = RJigsawVariables["RJVars_G_1_Jet1_pT"] ;
  rjigsawntv.RJVars_G_1_Jet2_pT       = RJigsawVariables["RJVars_G_1_Jet2_pT"] ;
  rjigsawntv.RJVars_QCD_dPhiR         = RJigsawVariables["RJVars_QCD_dPhiR"] ;
  rjigsawntv.RJVars_QCD_Rpt           = RJigsawVariables["RJVars_QCD_Rpt"] ;
  rjigsawntv.RJVars_QCD_Rmsib         = RJigsawVariables["RJVars_QCD_Rmsib"] ;
  rjigsawntv.RJVars_QCD_Rpsib         = RJigsawVariables["RJVars_QCD_Rpsib"] ;
  rjigsawntv.RJVars_QCD_Delta1        = RJigsawVariables["RJVars_QCD_Delta1"] ;
  rjigsawntv.RJVars_QCD_Delta2        = RJigsawVariables["RJVars_QCD_Delta2"] ;
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
				   bool isTruth,
				   std::vector<TauProxy> baseline_taus,
				   std::vector<TauProxy> signal_taus)
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
  ntv.ttbarAvgPt = ttbarAvgPt * 0.001;
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
  ntv.MET = MissingEt * 0.001;
  ntv.METPhi = METPhi;
  ntv.deltaPhi = minDphi;
  ntv.deltaPhiRemaining=RemainingminDPhi;
  
  ntv.MeffIncl= meffincl * 0.001;
  ntv.normWeight=(normWeight.size() >= 1) ? normWeight.at(0) : 1.;
  ntv.normWeightUp=(normWeight.size() >= 2) ? normWeight.at(1) : 1.;
  ntv.normWeightDown=(normWeight.size() >= 3) ? normWeight.at(2) : 1.;

  ntv.hardproc=hardproc;
  ntv.cleaning=cleaning;
  ntv.timing=timing;

  ntv.SherpaBugMET = SherpaBugMET * 0.001;

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
      ntv.jetPt.push_back(pt * 0.001);
      ntv.jetEta.push_back(eta);
      ntv.jetPhi.push_back(phi);
      ntv.jetM.push_back(m * 0.001);
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
				  metNOCHCORRCELLy*metNOCHCORRCELLy) * 0.001;
  ntv.metLHTOPONOCHCORRCELL = std::sqrt(metLHTOPONOCHCORRCELLx*metLHTOPONOCHCORRCELLx+
					metLHTOPONOCHCORRCELLy*metLHTOPONOCHCORRCELLy) * 0.001;
  ntv.metLHTOPO = std::sqrt(metLHTOPOx*metLHTOPOx+
			    metLHTOPOy*metLHTOPOy) * 0.001;
  
  ntv.jetSmearSystW = jetSmearSystW;

  ntv.tauN = signal_taus.size();
  ntv.tauLooseN = baseline_taus.size();
  //std::cout<<"tauN="<<ntv.tauN<<" tauLooseN="<<ntv.tauLooseN<<std::endl;
  for ( size_t tau0=0; tau0<baseline_taus.size(); tau0++) 
  {
    const TauProxy& thistau = baseline_taus[tau0];
    ntv.tauPt.push_back(thistau.Pt() * 0.001);
    ntv.tauEta.push_back(thistau.Eta());
    ntv.tauPhi.push_back(thistau.Phi());
    float sf, sfStatUp, sfStatDown, sfSystUp, sfSystDown;
    thistau.getSF(sf,sfStatUp,sfStatDown,sfSystUp,sfSystDown);
    ntv.tauLooseSF.push_back(sf);
    ntv.tauLooseSFStatUp.push_back(sfStatUp);
    ntv.tauLooseSFStatDown.push_back(sfStatDown);
    ntv.tauLooseSFSystUp.push_back(sfSystUp);
    ntv.tauLooseSFSystDown.push_back(sfSystDown);
    //std::cout<<" tauLooseSF = "<<sf<<"+"<<sfStatUp<<"-"<<sfStatDown<<"+"<<sfSystUp<<"-"<<sfSystDown<<std::endl;
  }
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

  std::vector<float> pts =  myRT.recl_jets_Pt;
  for ( size_t i = 0; i < pts.size(); i++ ) {pts[i] *= 0.001;}
  std::vector<float> masses =  myRT.recl_jets_M;
  for ( size_t i = 0; i < masses.size(); i++ ) {masses[i] *= 0.001;}

  RTntv.RTjetM = masses;
  RTntv.RTjetPt = pts;
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
  
   for (UInt_t iJet=0; iJet< jets.size(); ++iJet){
       vec_pt[iJet]= jets[iJet].Pt();
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
