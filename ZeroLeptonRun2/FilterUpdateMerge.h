#ifndef FilterUpdateMerge_h
#define FilterUpdateMerge_h

namespace SUSY{
  class CrossSectionDB;
}
class TTree;
class NTVars;

#include <vector>
#include <string>

// helper class python that has trouble to handle vector<>
class FilterUpdateMergeFileList
{
public:
  FilterUpdateMergeFileList(): m_files() {}
  void add(const std::string& filename) { m_files.push_back(filename);}
  const std::vector<std::string>& files() const {return m_files;}
private:
  std::vector<std::string> m_files;
};

//--------------------------------------------------------------------------
// A tool to merge small-tuples produce by the factory code, with the option
// to filter event or recalculate some variables.
//--------------------------------------------------------------------------
class FilterUpdateMerge
{
 public:
  FilterUpdateMerge(SUSY::CrossSectionDB* xsecDB);

  // merge input tree inTreeName from file inFiles to the outTree, possibly
  // applying event filtering and event reweighting according to cross-section 
  void process(TTree* outTree, const std::string& inTreeName, 
	       const std::vector<std::string>& inFiles, 
	       bool isSignal, bool doXSecReweighting, bool doFiltering,
<<<<<<< HEAD
	       bool doExtraVars, bool doCRWTVars, bool doCRZVars, 
	       bool doCRYVars);
=======
	       bool doExtraVars, bool doRJigsawVars);
>>>>>>> Works for some number of the variables right now. going to fill out other variables now

  // trick for python that has trouble to handle vector<>
  void process(TTree* outTree, const std::string& inTreeName, 
	       const FilterUpdateMergeFileList& inFiles, 
	       bool isSignal, bool doXSecNormalisation, bool doFiltering,
<<<<<<< HEAD
	       bool doExtraVars, bool doCRWTVars, bool doCRZVars, 
	       bool doCRYVars)
  {
    process(outTree, inTreeName, inFiles.files(), isSignal, doXSecNormalisation, doFiltering, doExtraVars, doCRWTVars, doCRZVars, doCRYVars);
=======
	       bool doExtraVars, bool doRJigsawVars)
  {
    process(outTree, inTreeName, inFiles.files(), isSignal, doXSecNormalisation, doFiltering,doExtraVars,doRJigsawVars);
>>>>>>> Works for some number of the variables right now. going to fill out other variables now
  }

 private:
  bool acceptEvent(const NTVars& ntv) const;

  SUSY::CrossSectionDB* m_xsecDB;
  bool addExtraVars;
  bool addRJigsawVars;
  bool addCRZVars;
  bool addCRWTVars;
  bool addCRYVars;
  bool addTheoryVars;
  bool addISRVars;
  bool addWVars;
};



#endif
