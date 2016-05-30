#ifndef ZeroLeptonRun2_TileTriggerCheck_H_
#define ZeroLeptonRun2_TileTriggerCheck_H_

//--------------------------------------------------------------------------
//  Implement potential trigger tower satturation effect and mistriggering
//  fom https://twiki.cern.ch/twiki/bin/view/AtlasProtected/SusyBgForumRun2CheckList
// To make sure your analysis is not sensitive to L1Calo timing issue due 
// to saturation in Tile, plot the energy deposited in Tile for the jets of a 
// few signal points close to your exclusion limit 
// (jetem_E * (1 - f_em - f_hec)). If this gets up to around 700 GeV, please 
// contact the Background Forum conveners. 
//
//--------------------------------------------------------------------------


#include "cafe/Processor.h"
class TH1F;

class TileTriggerCheck : public cafe::Processor {
public:
  TileTriggerCheck(const char *name);
  void begin();
  bool processEvent(xAOD::TEvent& event);

private:
  TH1F* m_TileEnergy;

public:
  ClassDef(TileTriggerCheck,0);
};

#endif // ZeroLeptonRun2_TileTriggerCheck_H_
