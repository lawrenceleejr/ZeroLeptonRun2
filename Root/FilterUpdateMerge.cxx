#include "ZeroLeptonRun2/FilterUpdateMerge.h"
#include "ZeroLeptonRun2/ZeroLeptonNTVars.h"

#include "SUSYTools/SUSYCrossSection.h"
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"

#include <stdexcept>
#include <iostream>
#include <chrono>
#include <thread>

FilterUpdateMerge::FilterUpdateMerge(SUSY::CrossSectionDB* xsecDB):
  m_xsecDB(xsecDB),
  addExtraVars(false),
  addRJigsawVars(false),
  addCRZVars(false),
  addCRWTVars(false),
  addCRYVars(false),
  addTheoryVars(false),
  addISRVars(false),
  addWVars(false)
{
  std::cout << "Printing xsecDB " << std::endl;

  for(std::map<SUSY::CrossSectionDB::Key,
	       SUSY::CrossSectionDB::Process>::const_iterator dbEntry = xsecDB->begin();
       dbEntry != xsecDB->end();
       ++dbEntry) {
    auto key     = dbEntry->first;
    auto process = dbEntry->second;


    std::cout <<     "process              : " << process.name() << std::endl;
    std::cout <<     "process.ID()         : "          <<process.ID()            << std::endl;//const { return m_id                      ;}
    std::cout <<     "process.name()       : "         <<process.name()          << std::endl;//const { return m_name                      ;}
    std::cout <<     "process.xsect()      : "         <<process.xsect()         << std::endl   ;//const { return m_cross_section                      ;}
    std::cout <<     "process.kfactor()    : "       <<process.kfactor()       << std::endl     ;//const { return m_kfactor                      ;}
    std::cout <<     "process.efficiency() : "    <<process.efficiency()    << std::endl        ;//const { return m_efficiency                      ;}
    std::cout <<     "process.relunc()     : "    <<process.relunc()        << std::endl    ;//const { return m_relunc                      ;}
    std::cout <<     "process.sumweight()  : "     <<process.sumweight()     << std::endl       ;//const { return m_sumweight                      ;}
    std::cout <<     "process.stat()       : "          <<process.stat()          << std::endl  ;//const { return m_stat;}
    std::cout << std::endl;

  }

  if ( !m_xsecDB ) throw std::invalid_argument("FilterUpdateMerge: need a valid xsecDB pointer !");
}

void  FilterUpdateMerge::process(TTree* outTree, const std::string& inTreeName,
				 const std::vector<std::string>& inFiles,
				 bool isSignal, bool doXSecNormalisation,
				 bool doFiltering, bool doExtraVars,
				 bool doRJigsawVars, bool doCRWTVars,
				 bool doCRZVars, bool doCRYVars)
{
  // book output tuple variables
  NTVars outVars;
  bookNTVars(outTree,outVars,false);

  NTReclusteringVars outRTVars;
  bookNTReclusteringVars(outTree,outRTVars);

  NTExtraVars inExtraVars, outExtraVars;
  if ( doExtraVars ) bookNTExtraVars(outTree,outExtraVars);

  NTRJigsawVars inRJigsawVars, outRJigsawVars;
  if ( doRJigsawVars ) bookNTRJigsawVars(outTree,outRJigsawVars);

  NTCRWTVars inCRWTVars, outCRWTVars;
  if ( doCRWTVars ) bookNTCRWTVars(outTree,outCRWTVars);

  NTCRZVars inCRZVars, outCRZVars;
  if ( doCRZVars ) bookNTCRZVars(outTree,outCRZVars);

  NTCRYVars inCRYVars, outCRYVars;
  if ( doCRYVars ) bookNTCRYVars(outTree,outCRYVars);

  // read, update, filter, copy
  unsigned int previousRun = 0;
  float rel_uncertainty = 0.;
  float sumW = 0.;
  float xsec = 0.;
  NTVarsRead inVars;
  NTReclusteringVarsRead inRTVars;
  //NTCRYVarsRead inCRYVars;
  auto checkMod = [](int numToMod, int modulo){return (numToMod%modulo==0);};

  for ( std::size_t i = 0; i < inFiles.size(); ++i ){
    if(checkMod(i, 100))  std::cout << "Ran over " << i << " files."  << std::endl;
    //TFile* f = TFile::Open(inFiles[i].c_str());
    TFile* f = OpenWithRetries(inFiles[i]);
    if ( !f || f->IsZombie() ) throw std::runtime_error("Could not open "+inFiles[i]);
    TTree* tree = dynamic_cast<TTree*>(f->Get(inTreeName.c_str()));
    if ( !tree  ) throw std::runtime_error("Could not find a tree"+inTreeName+inFiles[i]);
    inVars.setAddresses(tree);
    inRTVars.setAddresses(tree);
    if ( doExtraVars ) tree->GetBranch("NTExtraVars")->SetAddress(&inExtraVars.mettrack);
    if ( doRJigsawVars ) tree->GetBranch("NTRJigsawVars")->SetAddress(&inRJigsawVars.RJVars_PP_Mass);
    if ( doCRWTVars ) tree->GetBranch("NTCRWTVars")->SetAddress(&inCRWTVars.lep1Pt);
    if ( doCRZVars ) tree->GetBranch("NTCRZVars")->SetAddress(&inCRZVars.lep1Pt);
    //if ( doCRYVars) inCRYVars.setAddresses(tree);
    if ( doCRYVars) tree->GetBranch("NTCRYVars")->SetAddress(&inCRYVars.phPt);

    // loop over entries
    for ( size_t j = 0; j < static_cast<size_t>(tree->GetEntries()); ++j ) {
      tree->GetEntry(j);

      // Filter first to skip the unnecessary next block of copies 

      if ( doFiltering  && !acceptEvent(inVars.ntv) ) continue;

      // This block is probably somewhat slow. Could investigate speeding up.

      outVars = inVars.ntv;
      outRTVars = inRTVars.RTntv;
      //outCRYVars = inCRYVars.ntv;
      if ( doExtraVars) outExtraVars = inExtraVars;
      if ( doRJigsawVars) outRJigsawVars = inRJigsawVars;
      if ( doCRWTVars) outCRWTVars = inCRWTVars;
      if ( doCRZVars) outCRZVars = inCRZVars;
      if ( doCRYVars) outCRYVars = inCRYVars;


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

	if(sumW <= 0 ) {
	  std::cout << "SumW <= 0!!! Error in : " << std::endl;
	  std::cout << __FILE__ << " in function : " << __PRETTY_FUNCTION__ << " at line :"  << __LINE__ << std::endl;
	}
	if(xsec <= 0 ) {
	  std::cout << "xsec <= 0!!! Error in : " << std::endl;
	  std::cout << __FILE__ << " in function : " << __PRETTY_FUNCTION__ << " at line :"  << __LINE__ << std::endl;
	}

	outVars.normWeight     = (sumW != 0.) ? xsec/sumW : 0.;
        outVars.normWeightUp   = outVars.normWeight*(1 + rel_uncertainty);
        outVars.normWeightDown = outVars.normWeight*(1 - rel_uncertainty);
	if(checkMod(j, 10000)) {
	  std::cout << "xsec  : " << xsec << std::endl;
	  std::cout << "SumW  : " << sumW << std::endl;
	  std::cout << "xsec/sumW :" << xsec/sumW << std::endl;
	  std::cout << "outVars : " << outVars.normWeight << std::endl;
	}

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
  if (ntv.MeffIncl < 800.) return false;
  //if (ntv.jetPt.size() < 2) return false;
  //if (ntv.jetPt[1] < 60000.) return false;
  if (ntv.MET < 150.) return false;
  if (ntv.jetPt[0] < 100.) return false;

  return true;
}

TFile* FilterUpdateMerge::OpenWithRetries(const std::string& name)
{
  static std::vector<int> waitTimeInSec = {10,30,120,300,1200};
  // if there is a : in the name, it's a file on a server so it may be worth
  // retrying opening the file.
  bool canRetry = name.find(':') !=  std::string::npos;
  if ( canRetry ) {
    for ( auto delay :  waitTimeInSec ){
      TFile* f = TFile::Open(name.c_str(),"READ");
      if (f && !f->IsZombie()) {
	return f;
      }
      else {
	std::cout << "Could not open file, will retry after " << delay << " seconds" << std::endl;
	std::this_thread::sleep_for (std::chrono::seconds(delay));
      }
    }
    return 0;
  }
  else {
    return TFile::Open(name.c_str(),"READ");
  }
}

