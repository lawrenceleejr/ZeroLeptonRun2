#ifndef QCDSEEDSELECTION_H
#define QCDSEEDSELECTION_H

class QCDSeedSelection {

 public:
  bool passTrigger_2011Data(int RunNumber, bool EF_j180_a4_EFFS, bool EF_j135_a4_EFFS, bool EF_j100_a4_EFFS, bool EF_j75_a4_EFFS, bool EF_j55_a4_EFFS, bool EF_j240_a4tc_EFFS,  bool EF_j180_a4tc_EFFS,  bool EF_j135_a4tc_EFFS,  bool EF_j100_a4tc_EFFS,  bool EF_j75_a4tc_EFFS,  bool EF_j55_a4tc_EFFS, int nJets, float leadingJetPt, float& trigger_weight);
  bool passTrigger_2012Data(int RunNumber, bool EF_j460_a4tchad, bool EF_j360_a4tchad, bool EF_j280_a4tchad, bool EF_j220_a4tchad, bool EF_j180_a4tchad, bool EF_j145_a4tchad,bool EF_j110_a4tchad,bool EF_j80_a4tchad, bool EF_j55_a4tchad,int nJets, float leadingJetPt, float& trigger_weight);
  bool pass2B35JetTrigger_2012Data(int RunNumber, bool EF_2b35_j145_j35, bool EF_2b35_j145_j100 ,bool EF_2b35_j110_2j35 ,bool EF_2b35_3j35, int nJets, float leadingJetPt, float& trigger_weight);

    bool passBJetTrigger_2012Data(int RunNumber,bool EF_b360, bool EF_b280, bool EF_b220, bool EF_b180,bool EF_b110, bool EF_b80, bool EF_b55 , int nJets, float leadingJetPt, float& trigger_weight);
  bool selectSeedEvent(double MissingET, double sumET, double etMissSigCut);

  static const double GeV;
};
#endif
