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

enum ZeroLeptonDerivationTag {NotADerivation, p1872, p2353, p2363, p2372, p2375, p2377, p2384, p2419, INVALID_Derivation};

inline ZeroLeptonDerivationTag derivationTagFromString(const std::string& tag)
{
  if ( tag == "" || tag =="None" || tag == "none" || tag == "NA") return NotADerivation;
  if ( tag == "p1872" ) return p1872;
  if ( tag == "p2353" ) return p2353;
  if ( tag == "p2363" ) return p2363;
  if ( tag == "p2372" ) return p2372;
  if ( tag == "p2375" ) return p2375;
  if ( tag == "p2377" ) return p2377;
  if ( tag == "p2384" ) return p2384;
  if ( tag == "p2419" ) return p2419;
  return INVALID_Derivation;
}



#endif
