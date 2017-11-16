/*************************************************************************
 * mcc_creator.h - Converts multipoles to multipole charge clusters.     *
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

#ifndef MCC_CREATOR_H
#define MCC_CREATOR_H

#include <openbabel/parsmart.h>

#include <Eigen/Dense>

#include "molecule.h"

#include "smarts_matcher.h"

class MCCCreator 
{
public:
  MCCCreator() {}
  void setMolecule(GeneralMultipoledMolecule* m) { mol = m; }
  void ignoreDipole(bool f) { iDipole = f; }
  void ignoreQuadrupole(bool f) { iQuadrupole = f; }
  void setRadius(const double r) { rx = ry = rz = r/bohr; }
  void setRadius(const double x, const double y, const double z) { rx = x/bohr; ry = y/bohr; rz = z/bohr; }
  void setEps(const double e) { eps = e; }
  void installMCCs(const std::vector<std::string>& patterns)
  {
    if (iDipole && iQuadrupole) return; // nothing to do
    
    SMARTSMatcher matcher;
    matcher.addPatterns(patterns);
    std::vector< int > mcc_matches = matcher.getMatches(*mol);
    
    if (mcc_matches.empty()) return; // nothing to do
    
    // since OBMol looses charges when adding atoms back the up
    std::vector<double> charges;

    mol->createOrientMap();
    mol->calcOrientMatrices();
    charges.reserve(mol->NumAtoms()+1);
    charges.push_back(0); // move the beginning of the array to 1
    for (int i = 1; i <= mol->NumAtoms(); ++i) charges.push_back(mol->GetAtom(i)->GetPartialCharge());

    
    mol->BeginModify();
    for (int i = 0; i < mcc_matches.size(); ++i)
    {
      int idx = mcc_matches[i];
      OpenBabel::OBAtom* a = mol->GetAtom(idx);
      Multipoles *m = static_cast<Multipoles *>(a->GetData("Multipoles"));
      if ( !m ) continue; // no multipoles on this atom, so nothing to do
      if ( !( m->hasDipole() || m->hasQuadrupole() ) ) continue;
      
      double monopole = m->monopole().value();
      Vector3 dipole; if ( m->hasDipole() && (!iDipole) ) dipole = m->dipole().value();
      Matrix3 quadrupole; if ( m->hasQuadrupole() && (!iQuadrupole) ) quadrupole = m->quadrupole().value();
      
      if (dipole.length() < eps && fabs(quadrupole(2,2)) < eps) continue;
      
      Matrix3 orient = mol->matrix(idx);
      Eigen::Matrix<double, 7, 7> matrixQ;
      matrixQ <<
      //  charge positions
      //  0       +x             -x             +y             -y             +z             -z
      1,       1,             1,             1,             1,             1,             1,
      0,      rx,           -rx,             0,             0,             0,             0,
      0,       0,             0,            ry,           -ry,             0,             0,
      0,       0,             0,             0,             0,            rz,           -rz,
      0,    rx * rx,       rx * rx,    -.5 * ry * ry, -.5 * ry * ry, -.5 * rz * rz, -.5 * rz * rz,
      0, -.5 * rx * rx, -.5 * rx * rx,    ry * ry,       ry * ry,    -.5 * rz * rz, -.5 * rz * rz,
      0,       0,             0,             0,             0,             0,             0;
      Eigen::Matrix<double, 7, 1> vectorB;
      vectorB << monopole, dipole.x(), dipole.y(), dipole.z(), quadrupole(0,0), quadrupole(1,1), 0.0 ;
        
      if ( dipole.length() < eps ) // no dipole
      {
        // zeroing charge with smallest component
        double q[3] = {fabs(quadrupole(0,0)), fabs(quadrupole(1,1)), fabs(quadrupole(2,2))};
        if ( q[0] <= q[1] && q[0] <= q[2] ) matrixQ(6,1) = 1; else
          if ( q[1] <= q[0] && q[1] <= q[2] ) matrixQ(6,3) = 1; else
            if ( q[2] <= q[1] && q[2] <= q[0] ) matrixQ(6,5) = 1;
      } else if ( fabs(quadrupole(2,2)) < eps ) // no quadrupole
      {
        // zeroing charge with zero component
        double d[3] = {fabs(dipole(0)), fabs(dipole(1)), fabs(dipole(2))};
        if ( d[0] < eps ) matrixQ(6,1) = 1; else
          if ( d[1] < eps ) matrixQ(6,3) = 1; else
            if ( d[2] < eps ) matrixQ(6,5) = 1; else
              matrixQ(6,0) = 1; // dipole has all components > 0, so zeroing central charge
      } else // both dipole and quadrupole present so dipole controls everything
      {
        // zeroing charge with smallest component
        double d[3] = {fabs(dipole(0)), fabs(dipole(1)), fabs(dipole(2))};
        double q[3] = {fabs(quadrupole(0,0)), fabs(quadrupole(1,1)), fabs(quadrupole(2,2))};
        // check if dipole has all components
        if ( d[0] > eps && d[1] > eps && d[2] > eps ) matrixQ(6,0) = 1; else
          // check if dipole has a single component
          if ( d[0] < eps && d[1] < eps ) // dipole along Z
            if ( q[0] < q[1] ) matrixQ(6,1) = 1; else matrixQ(6,3) = 1; else
              if ( d[0] < eps && d[2] < eps ) // dipole along Y
                if ( q[0] < q[2] ) matrixQ(6,1) = 1; else matrixQ(6,5) = 1; else
                  if ( d[1] < eps && d[2] < eps ) // dipole along X
                    if ( q[1] < q[2] ) matrixQ(6,3) = 1; else matrixQ(6,5) = 1; else
                      // diple has a single zero component
                      if ( d[0] < eps ) matrixQ(6,1) = 1; else
                        if ( d[1] < eps ) matrixQ(6,3) = 1; else
                          if ( d[2] < eps ) matrixQ(6,5) = 1; else
                            matrixQ(6,0) = 1; // we should never get here, but just in case ...
      }
      Eigen::Matrix<double, 7, 1> q = matrixQ.fullPivLu().solve(vectorB);
       
      OpenBabel::OBAtom* central = mol->GetAtom(idx);
      //central->SetPartialCharge(q(0));
      charges[idx] = q(0);
      int vecs[7] = {0, 0, 0, 1, 1, 2, 2};
      double rv[7] = {0, rx*bohr, -rx*bohr, ry*bohr, -ry*bohr, rz*bohr, -rz*bohr};
      for (int j = 1; j < 7; ++j)
      {
        if (fabs(q(j)) > eps)
        {
          OpenBabel::OBAtom *ep = mol->NewAtom();
          OpenBabel::vector3 v =
            OpenBabel::vector3(orient(0,vecs[j]), orient(1,vecs[j]), orient(2,vecs[j])) * rv[j];
          ep->SetVector(a->GetVector() + v);
          ep->SetAtomicNum(0);
          //ep->SetPartialCharge(q(i));
          charges.push_back(q(j));
          mol->AddBond(idx, ep->GetIdx(), 1);
        }
      }
    }
    mol->EndModify();
    
    for (int i = 1; i <= mol->NumAtoms(); ++i) mol->GetAtom(i)->SetPartialCharge(charges[i]);
    mol->SetPartialChargesPerceived();
  }
private:
  GeneralMultipoledMolecule* mol;
  double rx, ry, rz;
  bool iDipole;
  bool iQuadrupole;
  double eps;
};

#endif // MCC_CREATOR_H