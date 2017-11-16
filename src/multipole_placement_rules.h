/*************************************************************************
 * multipole_placement_rules.h - Handles multipole placement rules.      *
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

#ifndef MULTIPOLES_PLACEMENT_RULES_HEADER
#define MULTIPOLES_PLACEMENT_RULES_HEADER

#include <iostream>
#include <string>
#include <cctype>
#include <stdexcept>

#include "vector3.h"
#include "multipoles.h"
#include "comment_skipper.h"
#include "mult_tools_defines.h"

class MultipoleFlags
{
public:
  MultipoleFlags() : value(0), dipole_dir(Vector3(0.0, 0.0, 0.0)) { }

  MultipoleFlags(const std::string s) : value(0), dipole_dir(Vector3(0.0, 0.0, 0.0))
  {
    if (!s.empty())
    {
      std::istringstream strm(s);
      readMe(strm);
    }
  }

  MultipoleFlags(unsigned int v) : value(v), dipole_dir(Vector3(0.0, 0.0, 0.0)) { }

  MultipoleFlags(const MultipoleFlags& other) : value(other.value), dipole_dir(other.dipole_dir) { }

  bool test(MULTIPOLES v) const { return value & v; }

  void clear() { value = NONE; dipole_dir = Vector3(0.0, 0.0, 0.0); }

  void set(MULTIPOLES v) { value |= v; }

  void unset(MULTIPOLES v) { value &= ~v; }

  const Vector3& dipoleDirection() const { return dipole_dir; }

#define SWITCH_PREVC(D,Q)            \
        switch (prevc)               \
        {                            \
          case 'd' :                 \
            set(DIPOLE_DIRECTED);    \
            dipole_dir = D;          \
            break;                   \
          case 'q' :                 \
            set(Q);     \
            break;                   \
          default :                  \
            throw std::domain_error(std::string("Reading placement rules: Unknown multipole: \"") + \
                  prevc + "\" followed by \"" + c + "\""); \
        }                            
  std::string readMe(std::istream& strm, std::string header = "")
  {
    std::string buff;
    if (strm.good()) strm >> buff;
    std::istringstream sstrm(buff);
    char c = '\0', prevc = '\0';
    bool expect_vector = false;
    while(sstrm)
    {
      prevc = c;
      c = '\0';
      sstrm >> c;
      switch (c)
      {
        case 'm' :
          set(MONOPOLE);
          break;
        case 'd' :
          set(DIPOLE);
          break;
        case 'q' :
          set(QUADRUPOLE);
          break;
        case '*' :
        case 'z' :
          SWITCH_PREVC(Vector3(0.0, 0.0, 1.0), QUADRUPOLE_Z)
          break;
        case 'x' :
          SWITCH_PREVC(Vector3(1.0, 0.0, 0.0), QUADRUPOLE_X)
          break;
        case 'y' :
          SWITCH_PREVC(Vector3(0.0, 1.0, 0.0), QUADRUPOLE_Y)
          break;
        case 'v' :
          if (prevc == 'd') expect_vector = true; 
          else throw std::domain_error(std::string("Reading placement rules: Unknown rule: \"v\" followed by \"") + prevc + "\"");
          break;
        case 'n' :
          value = NONE;
          break;
        case '\0':
          break;
        default :
          throw std::domain_error(std::string("Reading placement rules: Unknown multipole: \"") + c + "\"");
      }
    }
    if (expect_vector) 
    {
      Vector3 v; 
      std::string result = v.readMe(strm);
      if (result.empty()) 
      {
        v.normalize();
        dipole_dir = v;
        set(DIPOLE_DIRECTED);
      } else throw std::domain_error("Reading placement rules: Vector direction not found");
    }
    return std::string();
  }

  void printMe(std::ostream& strm) const
  {
    char dipole_flag = 'v';
    if ( (dipole_dir - Vector3(1.0, 0.0, 0.0)).length_squared() < 0.0001 ) dipole_flag = 'x'; else
    if ( (dipole_dir - Vector3(0.0, 1.0, 0.0)).length_squared() < 0.0001 ) dipole_flag = 'y'; else
    if ( (dipole_dir - Vector3(0.0, 0.0, 1.0)).length_squared() < 0.0001 ) dipole_flag = 'z'; 

    if (test(MONOPOLE)) strm << "m";
    if (test(DIPOLE)) 
    {
      strm << "d";
      if (test(DIPOLE_DIRECTED)) strm << dipole_flag;
    }
    if (test(QUADRUPOLE)) strm << "q";
    {
      if (test(QUADRUPOLE_X)) strm << "x";
      if (test(QUADRUPOLE_Y)) strm << "y";
      if (test(QUADRUPOLE_Z)) strm << "z";
    }
    if (test(DIPOLE) && test(DIPOLE_DIRECTED) && dipole_flag == 'v') strm << " " << dipole_dir;
  }

private:
  unsigned int value;
  Vector3 dipole_dir;
};

DECLARE_IO_STREAM_OPERATORS(MultipoleFlags)

class MultipolePlacementRule
{
public:
  MultipolePlacementRule(const std::string& s, MultipoleFlags r) : rule(r), SMARTS(s) { }
  
  MultipoleFlags rule;
  std::string SMARTS;
};

class MultipolePlacementRules
{
public:

  std::string readMe(std::istream& strm, std::string header = "")
  {
    if (header.empty()) strm >> header;
    if (header != "Placement-rules:") return header;
    clear();
    while (strm.good())
    {
      std::string type;
      strm >> type;
      if ((type == "atom:") || (type == "group:"))
      {
        std::string pattern;
        strm >> pattern;
        if (pattern.at(0) != '\"' || pattern.at(pattern.length() - 1) != '\"') throw std::domain_error("Reading placement rules: Bad atom SMARTS: " + pattern);
        MultipoleFlags rule;
        strm >> rule;
        if (type == "atom:") atomRules.push_back(MultipolePlacementRule(pattern.substr(1, pattern.length() - 2), rule));
        else groupRules.push_back(MultipolePlacementRule(pattern.substr(1, pattern.length() - 2), rule));
      } else return type;
    }
    return std::string();
  }

  void clear()
  {
    atomRules.clear();
    groupRules.clear();
  }
  
  int numAtomRules() const { return atomRules.size(); }

  int numGroupRules() const { return groupRules.size(); }

  const MultipolePlacementRule& atomRule(int i) const { return atomRules[i]; }

  const MultipolePlacementRule& groupRule(int i) const { return groupRules[i]; }

  void printMe(std::ostream& strm) const
  {
    strm << "Placement-rules:" << std::endl;
    for (int i = 0; i < atomRules.size(); ++i)
      strm << "atom: \"" << atomRules[i].SMARTS << "\" " << atomRules[i].rule << std::endl;
    for (int i = 0; i < groupRules.size(); ++i)
      strm << "atom: \"" << groupRules[i].SMARTS << "\" " << groupRules[i].rule << std::endl;
  }

private:
  std::vector< MultipolePlacementRule > atomRules;
  std::vector< MultipolePlacementRule > groupRules;
};

DECLARE_IO_STREAM_OPERATORS(MultipolePlacementRules)

#endif // MULTIPOLES_PLACEMENT_RULES_HEADER

