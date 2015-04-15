#ifndef ZeroLeptonRun2_ZeroLeptonRunPeriod_h
#define ZeroLeptonRun2_ZeroLeptonRunPeriod_h

#include <string>

enum ZeroLeptonRunPeriod {p7tev,p8tev,p13tev,INVALID};

inline ZeroLeptonRunPeriod periodFromString(const std::string& period)
{
  if ( period == "p7tev" ) return p7tev;
  else if ( period == "p8tev" ) return p8tev;
  else if ( period == "p13tev" ) return p13tev;
  return INVALID;
}

enum ZeroLeptonDerivationTag {NotADerivation, p1872, INVALID_Derivation};

inline ZeroLeptonDerivationTag derivationTagFromString(const std::string& tag)
{
  if ( tag == "" || tag =="None" || tag == "none" || tag == "NA") return NotADerivation;
  if ( tag == "p1872" ) return p1872;
  return INVALID_Derivation;
}



#endif
