/*
 *  Class for parsing parameter files
 *
 *  Copyright (C) 2003, 2004, 2005, 2008 Max-Planck-Society
 *  Authors: Martin Reinecke, Reinhard Hell
 * Adapted for SUSY/ATLAS Sophie Henrot-Versille
 */

#ifndef ZeroLepton_ParamFile_H
#define ZeroLepton_ParamFile_H

#include <map>
#include <set>
#include <string>
#include <iostream>
#include "ZeroLeptonRun2/ZeroLeptonCxxUtils.h"
#include <stdexcept>

class ZeroLeptonParamFile
  {
  private:
    typedef std::map<std::string,std::string> params_type;
    params_type params;
    mutable std::set<std::string> read_params;
    bool verbose;

    std::string get_valstr(const std::string &key) const
      {
      params_type::const_iterator loc=params.find(key);
      if (loc!=params.end()) return loc->second;
      throw std::logic_error("key not found: "+key);
      return "";
      }

    bool param_unread (const std::string &key) const
      { return (read_params.find(key)==read_params.end()); }

  public:
    ZeroLeptonParamFile (const std::string &filename, bool verbose_=false)
      : verbose(verbose_)
      { parse_file (filename, params); }

    ZeroLeptonParamFile (const params_type &par)
      : params (par), verbose(false)
      {}

    ~ZeroLeptonParamFile ()
      {
      if (verbose)
        for (params_type::const_iterator loc=params.begin();
             loc!=params.end(); ++loc)
          if (param_unread(loc->first))
            std::cout << "Parser warning: unused parameter "
                      << loc->first << std::endl;
      }

    bool param_present(const std::string &key) const
      { return (params.find(key)!=params.end()); }

    template<typename T> T find (const std::string &key) const
      {
      T result;
      stringToData(get_valstr(key),result);
      if (verbose && param_unread(key))
        std::cout << "Parser: " << key << " = " << dataToString(result)
                  << std::endl;
      read_params.insert(key);
      return result;
      }
    template<typename T> T find
      (const std::string &key, const T &deflt)
      {
      if (param_present(key)) return find<T>(key);
      if (verbose && param_unread(key))
        std::cout << "Parser: " << key << " = " << dataToString(deflt)
                  << " <default>" << std::endl;
      params[key]=dataToString(deflt);
      read_params.insert(key);
      return deflt;
      }

    const params_type &getParams() const
      { return params; }

    template<typename T> void findParam
      (const std::string &key, T &value) const
      { value = find<T>(key); }
  };

#endif
