#include "ZeroLeptonRun2/FilterUpdateMerge.h"
#include "ZeroLeptonRun2/ZeroLeptonNTVars.h"

#include "SUSYTools/SUSYCrossSection.h"
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"

#include <stdexcept>
#include <iostream>

FilterUpdateMerge::FilterUpdateMerge(SUSY::CrossSectionDB* xsecDB): 
  m_xsecDB(xsecDB),
  addExtraVars(false),
  addCRZVars(false),
  addCRWTVars(false),
  addCRYVars(false),
  addTheoryVars(false),
  addISRVars(false),
  addWVars(false)
{
  if ( !m_xsecDB ) throw std::invalid_argument("FilterUpdateMerge: need a valid xsecDB pointer !");
}

void  FilterUpdateMerge::process(TTree* outTree, const std::string& inTreeName, 
				 const std::vector<std::string>& inFiles, 
				 bool isSignal, bool doXSecNormalisation, 
				 bool doFiltering, bool doExtraVars, 
				 bool doCRWTVars, bool doCRZVars, 
				 bool doCRYVars)
{
  // book output tuple variables
  NTVars outVars;
  bookNTVars(outTree,outVars,false);
  
  NTReclusteringVars outRTVars;
  bookNTReclusteringVars(outTree,outRTVars);
  
  NTExtraVars inExtraVars, outExtraVars;
  if ( doExtraVars ) bookNTExtraVars(outTree,outExtraVars);

  NTCRWTVars inCRWTVars, outCRWTVars;
  if ( doCRWTVars ) bookNTCRWTVars(outTree,outCRWTVars);

  NTCRZVars inCRZVars, outCRZVars;
  if ( doCRZVars ) bookNTCRZVars(outTree,outCRZVars);

  NTCRYVars outCRYVars;
  if ( doCRYVars ) bookNTCRYVars(outTree,outCRYVars);

  // read, update, filter, copy
  unsigned int previousRun = 0;
  float rel_uncertainty = 0.;
  float sumW = 0.;
  float xsec = 0.;
  NTVarsRead inVars;
  NTReclusteringVarsRead inRTVars;
  NTCRYVarsRead inCRYVars;
  for ( std::size_t i = 0; i < inFiles.size(); ++i ){
    TFile* f = TFile::Open(inFiles[i].c_str());
    if ( !f || f->IsZombie() ) throw std::runtime_error("Could not open "+inFiles[i]);
    TTree* tree = dynamic_cast<TTree*>(f->Get(inTreeName.c_str()));
    if ( !tree  ) throw std::runtime_error("Could not find a tree named SRAllNT in "+inFiles[i]);
    inVars.setAddresses(tree);
    inRTVars.setAddresses(tree);
    if ( doExtraVars ) tree->GetBranch("NTExtraVars")->SetAddress(&inExtraVars.mettrack);
    if ( doCRWTVars ) tree->GetBranch("NTCRWTVars")->SetAddress(&inCRWTVars.lep1Pt);
    if ( doCRZVars ) tree->GetBranch("NTCRZVars")->SetAddress(&inCRZVars.lep1Pt);
    if ( doCRYVars) inCRYVars.setAddresses(tree);

    // loop over entries
    for ( size_t j = 0; j < tree->GetEntries(); ++j ) {
      tree->GetEntry(j);
      outVars = inVars.ntv;
      outRTVars = inRTVars.RTntv;
      outCRYVars = inCRYVars.ntv;
      if ( doExtraVars) outExtraVars = inExtraVars;
      if ( doCRWTVars) outCRWTVars = inCRWTVars;
      if ( doCRZVars) outCRZVars = inCRZVars;

      if ( doFiltering  && !acceptEvent(inVars.ntv) ) continue;

      if ( doXSecNormalisation ) {
	if ( isSignal ) { 
	  // cross-section to be recomputed for every event as the hard process changes
	  rel_uncertainty = m_xsecDB->rel_uncertainty(inVars.ntv.RunNumber,inVars.ntv.hardproc);
	  sumW = m_xsecDB->sumweight(inVars.ntv.RunNumber,inVars.ntv.hardproc);
	  xsec = m_xsecDB->xsectTimesEff(inVars.ntv.RunNumber,inVars.ntv.hardproc);
	}
	else {
	  if ( inVars.ntv.RunNumber != previousRun ) {
	    rel_uncertainty = m_xsecDB->rel_uncertainty(inVars.ntv.RunNumber);
	    xsec = m_xsecDB->xsectTimesEff(inVars.ntv.RunNumber);
	    if ( ( inVars.ntv.RunNumber >= 147910 and inVars.ntv.RunNumber <= 147917 ) ||
		 ( inVars.ntv.RunNumber >= 361020 and inVars.ntv.RunNumber <= 361032 )   )  {
	      std::cout << " JZxW sample, normalize with # of evts instead of sum weights " << std::endl;
	      sumW = m_xsecDB->process(inVars.ntv.RunNumber).stat();
	    }
	    else {
	      sumW = m_xsecDB->sumweight(inVars.ntv.RunNumber);
	    }
	    previousRun = inVars.ntv.RunNumber;
	  }
	}
        outVars.normWeight     = (sumW != 0.) ? xsec/sumW : 0.;
        outVars.normWeightUp   = outVars.normWeight*(1 + rel_uncertainty);
        outVars.normWeightDown = outVars.normWeight*(1 - rel_uncertainty);
      }

      //std::cout << inVars.ntv.jetPt.size() << " " <<  outVars.jetPt.size() << std::endl;
      outTree->Fill();
    }
    f->Close();
    delete f;
  }

}

bool FilterUpdateMerge::acceptEvent(const NTVars& ntv) const
{
  // Custom function which can be used to filter events
  if (ntv.MeffIncl < 700000.) return false;
  if (ntv.jetPt.size() < 2) return false;
  if (ntv.jetPt[1] < 60000.) return false;
  if (ntv.MET < 160000.) return false;
  if (ntv.jetPt[0] < 130000.) return false;
  
  return true;
}
