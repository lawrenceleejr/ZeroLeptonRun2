#ifndef ZeroLeptonRun2_CleaningHelper_H_
#define ZeroLeptonRun2_CleaningHelper_H_

// -------------------------------------------------------------------
//  Cleaninghelper:
//    helper class for MCEventVeto to work around the limitations of
//    xAOD EDM vs CINT
// -------------------------------------------------------------------
#include <map>
#include <string>
#include <bitset>

class CleaningHelper
{
 public:
  CleaningHelper() {
    //the list of potential cleaning cuts
    //if we need another cleaning, add it here
    cleaning["badJetVeto"]             = false;
    cleaning["badMetMuonVeto"]         = false;
    cleaning["badTileVeto"]            = false;
    cleaning["negativeEnergyCellVeto"] = false;
    cleaning["leadingJetTimingVeto"]   = false;
    cleaning["chfTileVeto"]            = false;
    cleaning["chfVeto"]                = false;
    cleaning["metTSTCleaningVeto"]     = false;
    cleaning["badMuonVeto"]            = false;
    cleaning["cosmicMuonVeto"]         = false;
  };

  unsigned long finalCleaning(){
    assert(32 > cleaning.size());//require that we don't have too many cleaning cuts
    std::bitset< 32 > cleaningBitset;

    int counter = 0;
    for(std::map<std::string,bool>::const_iterator nameBoolPair = cleaning.begin();
       nameBoolPair != cleaning.end();
       ++nameBoolPair) {
      auto name    = nameBoolPair->first;
      auto passCut = nameBoolPair->second;

      cleaningBitset[counter] = passCut;
      ++counter;
    }

    return cleaningBitset.to_ulong();
  }

  std::map<std::string, bool> cleaning;

};


#endif
