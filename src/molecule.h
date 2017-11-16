/*************************************************************************
 * molecule.h - Handles multipoled molecule.                             *
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

#ifndef MOLECULE_H
#define MOLECULE_H

#ifndef SWIG
#include <stdexcept>
#include <string>

#include <openbabel/mol.h>
#endif // SWIG

#include "vector3.h"
#include "matrix3.h"
#include "multipoles.h"
#ifndef SWIG
#include "comment_skipper.h"
#include "multipole_orient_rules.h"
#include "mult_tools_defines.h"
#endif // SWIG

class MultipoledMolecule : public OpenBabel::OBMol
{
public:
  MultipoledMolecule() { };

  MultipoledMolecule(const MultipoledMolecule& other) : OpenBabel::OBMol(other)
  {
    using OpenBabel::OBMolAtomIter;
    FOR_ATOMS_OF_MOL(a, this)
    {
      Multipoles *m = static_cast<Multipoles *>(a->GetData("Multipoles"));
      if (m) GetAtom(a->GetIdx())->SetData(new Multipoles(*m));
    }
  }

  void correctMultPolPositions()
  {
    using::OpenBabel::OBMolAtomIter;
    Multipoles *m;
    FOR_ATOMS_OF_MOL(a, this)
    {
      if (m = static_cast<Multipoles *>(a->GetData("Multipoles")))
      {
        m->moveTo(Vector3(a->GetCoordinate()));
      }
    }
  }

  void Translate(const OpenBabel::vector3& v)
  {
    using OpenBabel::OBMolAtomIter;
    OpenBabel::OBMol::Translate(v);
    correctMultPolPositions();
  }

  virtual void Rotate(const double u[3][3])
  {
    OpenBabel::OBMol::Rotate(u);
    correctMultPolPositions();
    calcOrientMatrices();
    prepareMultipoles();
  }

  virtual void Rotate(const double u[9])
  {
    OpenBabel::OBMol::Rotate(u);
    correctMultPolPositions();
    calcOrientMatrices();
    prepareMultipoles();
  }
  
  virtual void Rotate(const Matrix3& rot)
  {
    OpenBabel::OBMol::Rotate(rot.transposed().data());
    correctMultPolPositions();
    calcOrientMatrices();
    prepareMultipoles();
  }

  void moveAtomTo(OpenBabel::OBAtom *a, const Vector3& vec)
  {
    a->SetVector(OpenBabel::vector3(vec[0], vec[1], vec[2]));
  }

  void addMultipoles(OpenBabel::OBAtom *a, const Multipoles& mult)
  {
    a->SetData(new Multipoles(mult));
  }

  double E_El(const MultipoledMolecule& other) const
  {
    double result = 0.0;
    for (unsigned i = 1, ie = NumAtoms(); i <= ie; ++i)
      for (unsigned j = 1, je = other.NumAtoms(); j <= je; ++j)
      {
        const Multipoles *ma, *mb;
        ma = other.GetMultipoles(j);
        mb = GetMultipoles(i);
        if (ma && mb)
        {
          result += ma->energy(*mb);
        }
      }
    return result;
  }

  double MEP(const Vector3& pos) const
  {
    double result = 0.0;
    for (unsigned i = 1, ie = NumAtoms(); i <= ie; ++i)
    {
      const Multipoles *ma;
      ma = GetMultipoles(i);
      if (ma)
        result += ma->potential(pos);
    }
    return result;
  }

  
  virtual void calcOrientMatrices() = 0;

  virtual void prepareMultipoles() { };

  virtual const Multipoles *GetMultipoles(const unsigned int idx) const = 0;
  virtual const Multipoles *GetRawMultipoles(const unsigned int idx) const = 0;

  bool hasMultipoles(const unsigned int idx) { return GetRawMultipoles(idx); }

  void printMe(const std::string& filename) const
  {
    std::ofstream output(filename.c_str());
    if (output.good()) printMe(output);
  }

  void printMe(std::ostream& strm = std::cout) const
  {
    do_printMe(strm);
  }

  bool readMe(const std::string& filename )
  {
    std::ifstream input(filename.c_str());
    if (!input.good()) return false;
    CommentSkipper skipper;
    std::istringstream ifile(skipper.skip(input));
    if (!readMe(ifile).empty()) return false;
  }

  std::string readMe(std::istream& strm, std::string header = "")
  {
    return do_readMe(strm, header);
  }

  void raiseDummies() // since OBSMARTS does not match "[#0]" we change all dummies to element #99
  {
    unsigned z;
    for (int i = 1; i <= NumAtoms(); ++i)
    {
      z = GetAtom(i)->GetAtomicNum();
      if (z == 0) GetAtom(i)->SetAtomicNum(99);
    }
  }
  void lowerDummies() // since OBSMARTS does not match "[#0]" we change all dummies to element #99
  {
    unsigned z;
    for (int i = 1; i <= NumAtoms(); ++i)
    {
      z = GetAtom(i)->GetAtomicNum();
      if (z == 99) GetAtom(i)->SetAtomicNum(0);
    }
  }
  
  void syncronizeCharges()
  {
    double total_charge = 0;
    for (int i = 1; i <= NumAtoms(); ++i)
    {
      Multipoles* m = static_cast<Multipoles*>(GetAtom(i)->GetData("Multipoles"));
      if (m) 
        if (m->hasMonopole()) 
        {
          GetAtom(i)->SetPartialCharge(m->monopole().value());
          total_charge += m->monopole().value();
        } else GetAtom(i)->SetPartialCharge(0); 
    }
    SetTotalCharge(total_charge);
    SetPartialChargesPerceived(); 
  }
  void reverseSyncronizeCharges()
  {
    for (int i = 1; i <= NumAtoms(); ++i)
    {
      OpenBabel::OBAtom *a = GetAtom(i);
      a->DeleteData("Multipoles");
      Multipoles* m = new Multipoles(Vector3(a->GetX(), a->GetY(), a->GetZ()));
      m->setMonopole(Monopole(a->GetPartialCharge()));
      a->SetData(m);
    }
    prepareMultipoles();
  }
  
protected:
  virtual std::string do_readMe(std::istream& strm, std::string header = "")
  {
    std::string buff;
    if ( header.empty() ) strm >> header;
    if ( header != "Molecule:" ) return header;
    BeginModify();
    strm >> buff;
    while (strm.good())
    {
      if (buff == "Atom:")
      {
        unsigned int z;
        Vector3 vec;
        std::string tmp;
        strm >> tmp;
        std::istringstream ss(tmp);
        ss >> z;
        strm >> vec;
        OpenBabel::OBAtom *a = NewAtom();
        a->SetAtomicNum(z);
        a->SetVector(OpenBabel::vector3(vec.x(), vec.y(), vec.z()));
        Multipoles mult;
        strm >> buff;
        if (buff == "Multipoles:")
        {
          buff = mult.readMe(strm, buff);
          a->SetData(new Multipoles(mult));
          if (mult.hasMonopole()) a->SetPartialCharge(mult.monopole().value());
        }
      }
      else if (buff == "Bond:")
      {
        unsigned int from, to;
        char c1, c2, order;

        strm >> from >> c1 >> to >> c2 >> order;
        if (c1 != '-' || c2 != ':' || !( order == 'a' || order == '1' || order == '2' || order == '3' )) throw std::domain_error("Reading molecule: Unknown bond order: " + order);
        OpenBabel::OBBond *b = NewBond();
        OpenBabel::OBAtom *beg = GetAtom(from), *end = GetAtom(to);
        if (order == 'a')
        {
          b->SetBondOrder(5); b->SetAromatic();
          beg->SetAromatic(); end->SetAromatic();
        } else b->SetBondOrder(order - '1' + 1);
        b->SetBegin(beg);
        b->SetEnd(end);
        beg->AddBond(b);
        end->AddBond(b);
        strm >> buff;
      } else { break; }
    }
    EndModify();
    syncronizeCharges();    
    return buff;
  }

  virtual void do_printMe(std::ostream& strm) const
  {
    Multipoles *m;
    strm << "Molecule:" << std::endl;
    for (int i = 1, ie = NumAtoms(); i <= ie; ++i)   // OpenBabel indexes atoms from 1
    {
      OpenBabel::OBAtom *a = GetAtom(i);
      m = static_cast<Multipoles *>(a->GetData("Multipoles"));
      OpenBabel::vector3 tv = a->GetVector();
      Vector3 v(tv[0], tv[1], tv[2]);
      unsigned z = a->GetAtomicNum();
      strm << "/* Atom " << i << " " << OpenBabel::etab.GetSymbol((z == 99) ? 0 : z) << " */" << std::endl;
      strm << "Atom: " << ( (z == 99) ? 0 : z ) << " " << v << std::endl;
      if (m) strm << *m;
      strm << std::endl;
    }
    for (int i = 0, ie = NumBonds(); i < ie; ++i)
    {
      OpenBabel::OBBond *b = GetBond(i);
      strm << "Bond: " << b->GetBeginAtomIdx() << " - " << b->GetEndAtomIdx() << " : ";
      if (b->IsAromatic()) strm << "a"; else strm << b->GetBondOrder();
      strm << std::endl;
    }
  }
};

#ifndef SWIG
DECLARE_IO_STREAM_OPERATORS(MultipoledMolecule)
#endif // SWIG

class GeneralMultipoledMolecule : public MultipoledMolecule
{
public:
  GeneralMultipoledMolecule() { }

  GeneralMultipoledMolecule(const GeneralMultipoledMolecule& other)
      : MultipoledMolecule(other),
        orientMap(other.orientMap),
        matrices(other.matrices),
        orientedMultipoles(other.orientedMultipoles),
        orientRules(other.orientRules)
  {
  };
  
  void EndModify()
  {
    MultipoledMolecule::EndModify();
    createOrientMap();
    calcOrientMatrices();
    prepareMultipoles();
  }

  void setOrientRules(const MultipoleOrientRules& rules)
  {
    orientRules = rules;
    createOrientMap();
    calcOrientMatrices();
    prepareMultipoles();
  }
  const MultipoleOrientRules& getOrientRules() const
  {
    return orientRules;
  }

  void createOrientMap()
  {
    orientMap.resize(NumAtoms() + 1);
    SMARTS_hits.resize(NumAtoms() + 1);
    for (unsigned i = 0; i < orientRules.size(); ++i)
    {
      MultipoleOrientRule rule = orientRules.rule(i);
      OpenBabel::OBSmartsPattern SMARTS;
      SMARTS.Init(rule.SMARTS);
      SMARTS.Match(*this);
      std::vector< std::vector< int > > hits = SMARTS.GetMapList();
      for (int unsigned k = 0; k < hits.size(); ++k)
      {
        int idx[4] = { -1, -1, -1, -1 };
        for (unsigned l = 0; l < hits[k].size() && l < 4u; ++l) idx[l] = hits[k][l];
        orientMap[idx[0]] = OrientMapData(rule.scenario, idx[1], idx[2], idx[3]);
        SMARTS_hits[idx[0]] = rule.SMARTS;
      }
    }
  }
  
  const std::vector<std::string>& getOrientMatches() const { return SMARTS_hits; }

  void calcOrientMatrices()
  {
    matrices.resize(NumAtoms() + 1);
    Vector3 Zaxis, Yaxis, Xaxis;
    for (unsigned i = 1; i < orientMap.size(); ++i)
    {
      switch (orientMap[i].type)
      {
        case 'a':
        {
          Zaxis = Vector3(GetAtom(i)->GetCoordinate()) - Vector3(GetAtom(orientMap[i].idx[0])->GetCoordinate());
          Zaxis.normalize();
          if (fabs(Zaxis.dot(Vector3(1.0, 0.0, 0.0))) < 1.0)
            Xaxis = Zaxis.vec(Vector3(1.0, 0.0, 0.0));
          else Xaxis = Zaxis.vec(Vector3(0.0, 1.0, 0.0));
          Yaxis = Xaxis.vec(Zaxis);
          break;
        }
        case 'b':
        {
          Zaxis = Vector3(GetAtom(i)->GetCoordinate()) - Vector3(GetAtom(orientMap[i].idx[0])->GetCoordinate());
          Xaxis = Zaxis.vec(Vector3(GetAtom(i)->GetCoordinate()) - Vector3(GetAtom(orientMap[i].idx[1])->GetCoordinate()));
          Yaxis = Xaxis.vec(Zaxis);
          break;
        }
        case 'c':
        {
          Vector3 b1 = Vector3(GetAtom(i)->GetCoordinate())
              - Vector3(GetAtom(orientMap[i].idx[0])->GetCoordinate());
          Vector3 b2 = Vector3(GetAtom(i)->GetCoordinate())
              - Vector3(GetAtom(orientMap[i].idx[1])->GetCoordinate());
          b1.normalize(); b2.normalize();
          Zaxis = b1 + b2;
          Xaxis = Zaxis.vec(Vector3(GetAtom(i)->GetCoordinate()) - Vector3(GetAtom(orientMap[i].idx[0])->GetCoordinate()));
          Yaxis = Xaxis.vec(Zaxis);
          break;
        }
        case 'f':
        {
          Vector3 b1 = Vector3(GetAtom(i)->GetCoordinate())
              - Vector3(GetAtom(orientMap[i].idx[0])->GetCoordinate());
          Vector3 b2 = Vector3(GetAtom(i)->GetCoordinate())
              - Vector3(GetAtom(orientMap[i].idx[1])->GetCoordinate());
          b1.normalize(); b2.normalize();
          Vector3 Zax = b1 + b2;
          Vector3 Xax = Zax.vec(Vector3(GetAtom(i)->GetCoordinate()) - Vector3(GetAtom(orientMap[i].idx[0])->GetCoordinate()));
          Vector3 Yax = Xax.vec(Zax);
          Zax.normalize(); Yax.normalize();
          Xaxis = Zax + Yax; Yaxis = (Zax - Yax).normalized(); Zaxis = Xax;
          break;
        }
        case 'd':
        {
          Vector3 cent(GetAtom(i)->GetCoordinate()), p1(GetAtom(orientMap[i].idx[0])->GetCoordinate()),
              p2(GetAtom(orientMap[i].idx[1])->GetCoordinate()), p3(GetAtom(orientMap[i].idx[2])->GetCoordinate());
          Zaxis = ( cent - p1 ).normalized() + ( cent - p2 ).normalized() + ( cent - p3 ).normalized();
          if (Zaxis.length() < 0.001)
          {
            Xaxis = cent - p1;
            Zaxis = Xaxis.vec(cent - p2);
          } else
          {
            Xaxis = Zaxis.vec(cent - p1);
          }
          Yaxis = Zaxis.vec(Xaxis);
          break;
        }
        case 'e':
        {
          Xaxis = Vector3(GetAtom(i)->GetCoordinate()) - Vector3(GetAtom(orientMap[i].idx[0])->GetCoordinate());
          Zaxis = Xaxis.vec(Vector3(GetAtom(i)->GetCoordinate()) - Vector3(GetAtom(orientMap[i].idx[1])->GetCoordinate()));
          Yaxis = Xaxis.vec(Zaxis);
          break;
        }
        default:
        case 'z':
        {
          Xaxis = Vector3(1.0, 0.0, 0.0);
          Yaxis = Vector3(0.0, 1.0, 0.0);
          Zaxis = Vector3(0.0, 0.0, 1.0);
          break;
        }
      } // switch
      Zaxis.normalize();
      Yaxis.normalize();
      Xaxis.normalize();
//      matrices[i] = Matrix3(Xaxis, Yaxis, Zaxis).transposed();
      matrices[i] = Matrix3(Xaxis, Yaxis, Zaxis);
    }
  }

  void prepareMultipoles()
  {
    for (unsigned i = 1, ie = NumAtoms() + 1; i < ie; ++i)
    {
      Multipoles *m = static_cast<Multipoles *>(GetAtom(i)->GetData("Multipoles"));
      if (m)
      {
        Multipoles m1(*m);
        m1.rotate(matrices[i]);
        orientedMultipoles[i] = m1;
      }
    }
  }

  virtual const Multipoles *GetMultipoles(const unsigned int idx) const
  {
    std::map< unsigned, Multipoles >::const_iterator i = orientedMultipoles.find(idx);
    if (i != orientedMultipoles.end())
    {
      return &( i->second );
    } else { return 0; }
  }

  virtual const Multipoles *GetRawMultipoles(const unsigned int idx) const
  {
    return static_cast<const Multipoles*>(GetAtom(idx)->GetData("Multipoles"));
  }

  const Matrix3& matrix(const unsigned int i) const
  {
    return matrices[i];
  }

protected:
  std::string do_readMe(std::istream& strm, std::string header = "")
  {
    if (header.empty()) strm >> header;
    MultipoleOrientRules oRules;
    if (header == "Molecule:")
    {
      header = MultipoledMolecule::do_readMe(strm, header);
      header = oRules.readMe(strm, header);
    } else if (header == "Orient-rules:")
    {
      header = oRules.readMe(strm, header);
      header = MultipoledMolecule::do_readMe(strm, header);
    }
    setOrientRules(oRules);
    return header;
  }
  void do_printMe(std::ostream& strm) const
  {
    MultipoledMolecule::do_printMe(strm);
    strm << std::endl;
    orientRules.printMe(strm);
  }

private:
  struct OrientMapData
  {
    char type;
    int idx[3];

    OrientMapData(char t = '\0', int a = 0, int b = 0, int c = 0)
    {
      type = t;
      idx[0] = a;
      idx[1] = b;
      idx[2] = c;
    }
  };

  // matrices to convert multipoles from local to global coordinate systems
  std::vector< OrientMapData > orientMap;
  std::vector< Matrix3 > matrices;  // Ox=X
  std::map< unsigned, Multipoles > orientedMultipoles;
  MultipoleOrientRules orientRules;
  // for orient rules debugging
  std::vector<std::string> SMARTS_hits;
};

#ifndef SWIG
DECLARE_IO_STREAM_OPERATORS(GeneralMultipoledMolecule)
#endif // SWIG

#ifdef SWIG
__STR__(GeneralMultipoledMolecule);
#endif // SWIG

#endif
