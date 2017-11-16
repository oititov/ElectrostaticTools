/*************************************************************************
 * multipole_setter.h - Modifies atomic multipoles' values.              *
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

#ifndef MULTIPOLE_SETTER_H
#define MULTIPOLE_SETTER_H

#include <boost/algorithm/string.hpp>

#include <sstream>

#include <openbabel/parsmart.h>

#include "molecule.h" 

class MultipoleSetter
{
public: 
  MultipoleSetter() {}
  void setMolecule(GeneralMultipoledMolecule* m) { mol = m; }
  void placeMonopole(std::string& s, bool add = false)
  {
    std::istringstream ss(s);
    std::string SMARTS;
    ss >> SMARTS;    
    Monopole v = readMonopole(ss);
    OpenBabel::OBSmartsPattern patt;
    patt.Init(SMARTS);
    patt.Match(*mol);
    std::vector< std::vector< int > > hits = patt.GetUMapList();
    for (int i = 0; i < hits.size(); ++i)
    {
      OpenBabel::OBAtom *a = mol->GetAtom(hits[i][0]);
      a->SetPartialCharge(v.value());
      mol->SetPartialChargesPerceived();
      Multipoles *m = static_cast<Multipoles*>(a->GetData("Multipoles"));
      if (m)
      {
        if (add) m->setMonopole(m->monopole().value() + v.value()); else 
        m->setMonopole(v);
      } else
      {
        m = new Multipoles(Vector3(a->GetX(), a->GetY(), a->GetZ()));
        m->setMonopole(v);
        a->SetData(m);
      }
    }
  }
  void placeDipole(std::string& s, bool add = false)
  {
    std::istringstream ss(s);
    std::string SMARTS;
    ss >> SMARTS;    
    Dipole v = readDipole(ss);
    OpenBabel::OBSmartsPattern patt;
    patt.Init(SMARTS);
    patt.Match(*mol);
    std::vector< std::vector< int > > hits = patt.GetUMapList();
    for (int i = 0; i < hits.size(); ++i)
    {
      OpenBabel::OBAtom *a = mol->GetAtom(hits[i][0]);
      Multipoles *m = static_cast<Multipoles*>(a->GetData("Multipoles"));
      if (m)
      {
        if (add) m->setDipole(m->dipole().value() + v.value()); else
        m->setDipole(v);
      } else
      {
        m = new Multipoles(Vector3(a->GetX(), a->GetY(), a->GetZ()));
        m->setMonopole(Monopole(a->GetPartialCharge()));
        m->setDipole(v);
        a->SetData(m);
      }
    }
  }
  void placeQuadrupole(std::string& s, bool add = false)
  {
    std::istringstream ss(s);
    std::string SMARTS;
    ss >> SMARTS;    
    Quadrupole v = readQuadrupole(ss);
    OpenBabel::OBSmartsPattern patt;
    patt.Init(SMARTS);
    patt.Match(*mol);
    std::vector< std::vector< int > > hits = patt.GetUMapList();
    for (int i = 0; i < hits.size(); ++i)
    {
      OpenBabel::OBAtom *a = mol->GetAtom(hits[i][0]);
      Multipoles *m = static_cast<Multipoles*>(a->GetData("Multipoles"));
      if (m)
      {
        if (add) m->setQuadrupole(m->quadrupole().value() + v.value()); else
        m->setQuadrupole(v);
      } else
      {
        m = new Multipoles(Vector3(a->GetX(), a->GetY(), a->GetZ()));
        m->setMonopole(Monopole(a->GetPartialCharge()));
        m->setQuadrupole(v);
        a->SetData(m);
      }
    }
  }
  
private:
  std::vector<std::string> prepareInput(std::string& s)
  {
    boost::algorithm::trim(s);
    std::vector<std::string> strs;
    boost::split(strs, "string to split", boost::is_any_of("\t \"\'"));
    return strs;
  }
  std::vector<double> getDV(std::istringstream& ss)
  {
    std::vector<double> dv;
    while (ss)
    {
      double v;
      ss >> v;
      dv.push_back(v);
    } 
    dv.pop_back(); // we've read the last value twice
    return dv;
  }
  Monopole readMonopole(std::istringstream& ss)
  {
    std::vector<double> dv = getDV(ss);
    Monopole result;
    if (dv.size() == 1)
    {
      result = Monopole(dv[0]);
    } else throw std::invalid_argument("Incorrect number of parameters specified for charge addition.");
    return result;
  }
  Dipole readDipole(std::istringstream& ss)
  {
    std::vector<double> dv = getDV(ss);
    Dipole result;
    if (dv.size() == 1)
    {
      result = Dipole(dv[0]);
    } else if (dv.size() == 3)
    {
      result = Dipole(dv[0], dv[1], dv[2]);
    } else throw std::invalid_argument("Incorrect number of parameters specified for dipole addition.");
    return result;
  }
  Quadrupole readQuadrupole(std::istringstream& ss)
  {
    std::vector<double> dv = getDV(ss);
    Quadrupole result;
    if (dv.size() == 1)
    {
      result = Quadrupole(dv[0]);
    } else if (dv.size() == 3)
    {
      result = Quadrupole(dv[0], dv[1], dv[2]);
    } else if (dv.size() == 6)
    {
      result = Quadrupole(dv[0], dv[1], dv[2], dv[3], dv[4], dv[5]);
    } else throw std::invalid_argument("Incorrect number of parameters specified for quadrupole addition.");
    return result;
  }
  GeneralMultipoledMolecule *mol;
};

#endif // MULTIPOLE_SETTER_H
