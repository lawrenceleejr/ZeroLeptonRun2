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

#endif
