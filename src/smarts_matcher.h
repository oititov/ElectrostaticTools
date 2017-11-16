/*************************************************************************
 * smarts_matcher.h - Performs SMARTS matches with several SMARTS        *
 *                    patterns at once                                   *
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

#ifndef SMARTS_MATCHER_H
#define SMARTS_MATCHER_H

#include <vector>
#include <string>

#include <openbabel/mol.h>
#include <openbabel/parsmart.h>

class SMARTSMatcher
{
public:
  SMARTSMatcher() {}
  void addPattern(const std::string& p) 
  { 
    OpenBabel::OBSmartsPattern patt; 
    patt.Init(p); 
    patterns.push_back(patt); 
  }
  void addPatterns(const std::vector<std::string> p)
  {
    for (std::vector<std::string>::const_iterator i = p.begin(), ie = p.end(); i != ie; ++i) addPattern(*i);
  }
  std::vector<int> getMatches(OpenBabel::OBMol& mol)
  {
    std::vector< int > matches;
    for (int i = 0; i < patterns.size(); ++i)
    {
      patterns[i].Match(mol);
      std::vector< std::vector< int > > hits = patterns[i].GetUMapList();
      for (int j = 0; j < hits.size(); ++j) matches.push_back(hits[j][0]);
    }
    std::sort(matches.begin(), matches.end());
    matches.erase(std::unique(matches.begin(), matches.end()), matches.end());
    return matches;
  }
  void clear() { patterns.clear(); }
private:
  std::vector< OpenBabel::OBSmartsPattern > patterns;
};

#endif // SMARTS_MATCHER_H