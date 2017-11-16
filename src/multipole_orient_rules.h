/*************************************************************************
 * multipole_orient_rules.h - Handles multipole orient rules.            *
 *                                                                       *
 * Copyright (C) 2014-2015 by Laboratory of Medicinal Chemistry of       *
 * Department of Chemistry, Lomonosov Moscow State University            *
 *                                                                       *
 * Programmed by Oleg. I. Titov                                          *
 *                                                                       *
 * This file is part of the Electrostatic Tools project.                 *
 * For more information, see <http://qsar.chem.msu.ru/elec_tools/>       *
 *                                                                       *
 * Based on work:                                                        *
 * O.I. Titov, D.A. Shulga, V.A. Palyulin, N.S. Zefirov                  *
 * Description of halogen bonding on the basis of multicenter multipole  *
 * expansion // Doklady Chemistry, 2013, vol. 450, p. 139                *
 *                                                                       *
 * This program is free software; you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation version 2 of the License.                *
 *                                                                       *
 * This program is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 * GNU General Public License for more details.                          *
 ************************************************************************/

#ifndef MULTIPOLES_ORIENT_RULES
#define MULTIPOLES_ORIENT_RULES

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "mult_tools_defines.h"
#include "comment_skipper.h"

struct MultipoleOrientRule
{
public:
  MultipoleOrientRule() { }
  MultipoleOrientRule(char scenario, const std::string& SMARTS) : SMARTS(SMARTS), scenario(scenario) { }

  std::string readMe(std::istream& strm, std::string header = "")
  {
    std::string pattern;
    if ( header.empty() ) strm >> header;
    if ( header != "rule:" ) throw std::domain_error("Reading orient rules: Expected \"rule:\" but found " + header);
    char c;
    strm >> c;
    if (c != 'z' && c != 'a' && c != 'b' && c != 'c' && c != 'd' && c != 'e' && c != 'f') 
       throw std::domain_error(std::string("Reading orient rules: Unknown rule type: ") + c);
    strm >> pattern;
    if (pattern.at(0) != '\"' || pattern.at(pattern.length() - 1) != '\"') 
       throw std::domain_error(std::string("Reading orient rules: Bad SMARTS: ") + pattern);
    scenario = c;
    SMARTS = pattern.substr(1, pattern.length() - 2);
    return std::string();
  }

  bool readMe(const std::string& filename )
  {
    std::ifstream input(filename.c_str());
    if (!input.good()) return false;
    CommentSkipper skipper;
    std::istringstream ifile(skipper.skip(input));
    if (!readMe(ifile).empty()) return false;
  }

  void printMe(std::ostream& strm) const
  {
    strm << "rule:\t" << scenario << "\t\"" << SMARTS << "\"";
  }

  bool operator == (const MultipoleOrientRule& other) const
  {
    return SMARTS == other.SMARTS; // ignore scenario
  }

  std::string SMARTS;
  char scenario;
};

#ifndef SWIG
DECLARE_IO_STREAM_OPERATORS(MultipoleOrientRule)
#endif // SWIG

#ifdef SWIG
__STR__(MultipoleOrientRule)
#endif // SWIG

class MultipoleOrientRules
{
public:

  int size() const { return rules.size(); }

  const MultipoleOrientRule& rule(const int i) const { return rules[i]; }

  std::string readMe(std::istream& strm, std::string header = "")
  {
    if (header.empty()) strm >> header;
    if (header != "Orient-rules:") return header;

    clear();
    while (strm.good())
    {
      std::string buff;
      MultipoleOrientRule r;
      strm >> buff;
      if (buff == "rule:") r.readMe(strm, buff); else return buff;
      rules.push_back(r);
    }
    return std::string();
  }

  bool readMe(const std::string& filename )
  {
    std::ifstream input(filename.c_str());
    if (!input.good()) return false;
    CommentSkipper skipper;
    std::istringstream ifile(skipper.skip(input));
    if (!readMe(ifile).empty()) return false;
  }
  
  void clear()
  {
    rules.clear();
  }

  void printMe(std::ostream& strm) const
  {
    strm << "Orient-rules:" << std::endl;
    for (int i = 0; i < rules.size(); ++i)
    {
      strm << rules[i] << std::endl;
    }
  }
  bool hasRule(const MultipoleOrientRule& r) const
  {
    return rules.end() != std::find(rules.begin(), rules.end(), r);
  }
  void addRule(const MultipoleOrientRule& r)
  {
    rules.push_back(r);
  }
  const MultipoleOrientRules& operator += (const MultipoleOrientRules& other)
  {
    for (int i = 0; i < other.size(); ++i)
    {
      MultipoleOrientRule r = other.rule(i);
      if ( ! hasRule(r) ) addRule(r);
    }
  }

private:
  std::vector< MultipoleOrientRule > rules;
};

#ifndef SWIG
DECLARE_IO_STREAM_OPERATORS(MultipoleOrientRules)
#endif // SWIG

#ifdef SWIG
__STR__(MultipoleOrientRules)
#endif // SWIG

#endif
