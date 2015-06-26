/*
 *  This file is part of Healpix_cxx.
 *
 *  Healpix_cxx is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  Healpix_cxx is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Healpix_cxx; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix, see http://healpix.jpl.nasa.gov
 */

/*
 *  Healpix_cxx is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*! \file cxxutils.h
 *  Various convenience functions used by the Planck LevelS package.
 *
 *  Copyright (C) 2002, 2003, 2004, 2005, 2006, 2007 Max-Planck-Society
 *  \author Martin Reinecke \author Reinhard Hell
 * adapted by Sophie Henrot-Versille for ATLAS 
 */

#ifndef ZeroLepton_CxxUTILS_H
#define ZeroLepton_CxxUTILS_H

#include <algorithm>
#include <string>
#include <map>
#include <cmath>

template<typename T> inline const char *type2typename ()
  { return "unknown type"; }
template<> inline const char *type2typename<signed char> ()
  { return "signed char"; }
template<> inline const char *type2typename<unsigned char> ()
  { return "unsigned char"; }
template<> inline const char *type2typename<short> ()
  { return "short"; }
template<> inline const char *type2typename<unsigned short> ()
  { return "unsigned short"; }
template<> inline const char *type2typename<int> ()
  { return "int"; }
template<> inline const char *type2typename<unsigned int> ()
  { return "unsigned int"; }
template<> inline const char *type2typename<long> ()
  { return "long"; }
template<> inline const char *type2typename<unsigned long> ()
  { return "unsigned long"; }
template<> inline const char *type2typename<long long> ()
  { return "long long"; }
template<> inline const char *type2typename<unsigned long long> ()
  { return "unsigned long long"; }
template<> inline const char *type2typename<float> ()
  { return "float"; }
template<> inline const char *type2typename<double> ()
  { return "double"; }
template<> inline const char *type2typename<bool> ()
  { return "bool"; }
template<> inline const char *type2typename<std::string> ()
  { return "std::string"; }

//! Returns the string \a orig without leading and trailing whitespace.
std::string trim (const std::string &orig);

//! Returns a string containing the text representation of \a x.
/*! Care is taken that no information is lost in the conversion. */
template<typename T> std::string dataToString(const T &x);
template<> std::string dataToString (const bool &x);
template<> std::string dataToString (const std::string &x);
template<> std::string dataToString (const float &x);
template<> std::string dataToString (const double &x);

/*! Returns a string containing the text representation of \a x, padded
    with leading zeroes to \a width characters. */
std::string intToString(int x, int width);

//! Reads a value of a given datatype from a string
template<typename T> void stringToData (const std::string &x, T &value);
template<> void stringToData (const std::string &x, std::string &value);
template<> void stringToData (const std::string &x, bool &value);

//! Reads a value of a given datatype from a string
template<typename T> inline T stringToData (const std::string &x)
  { T result; stringToData(x,result); return result; }

//! Parses the file \a filename and returns the key/value pairs in \a dict.
void parse_file (const std::string &filename,
  std::map<std::string,std::string> &dict);

//! Case-insensitive string comparison
/*! Returns \a true, if \a a and \a b differ only in capitalisation,
    else \a false. */
bool equal_nocase (const std::string &a, const std::string &b);

//! Returns lowercase version of \a input.
std::string tolower(const std::string &input);

/*! \} */

#endif
