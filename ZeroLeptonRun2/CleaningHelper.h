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
    cleaning["badJetVeto"] = 0;
    cleaning["badMetMuon"] = 0;
    cleaning["badTileVeto"] = 0;
    cleaning["negativeEnergyCellVeto"] = 0;
    cleaning["leadingJetTiming"] = 0;
    cleaning["chfTileVeto"] = 0;
    cleaning["chfVeto"] = 0;
    cleaning["failMetCleaning"] = 0;
    cleaning["badMuon"] = 0;
    cleaning["cosmicMuon"] = 0;
  };

  unsigned long finalCleaning(){
    std::bitset< cleaning.size() > cleaningBitset;

    int counter = 0;
    for(std::map<std::string,bool>::const_iterator nameBoolPair = cleaning->begin();
       nameBoolPair != cleaning->end();
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
