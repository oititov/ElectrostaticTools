/*************************************************************************
 * mep_fitter.h - Reference MEP fitter routines.                         *
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

#ifndef MEP_FITTER_H
#define MEP_FITTER_H

#include <fstream>
#include <sstream>
#include <vector>
#include <stdexcept>

#include <openbabel/chargemodel.h>

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <Eigen/QR>

#include <boost/tuple/tuple.hpp>

#include "grid.h"
#include "multipoles.h"
#include "molecule.h"
#include "multipole_orient_rules.h"
#include "multipole_placement_rules.h"
#include "comment_skipper.h"
#include "log_adder.h"


class MEPFitter : public LogAdder
{
public:
  MEPFitter(std::ostream& l) : LogAdder(l)
  {
    numRealAtoms = 0;
    numParams = 0;
    RSMD = 0.0;
    restrictEquivalent = true;
    recalcTopology = false;
    refit = false;
    has_charges = 0;  
    log() << "MEP Fitter initialized." << std::endl;
  }
  
  void setRestrictEquivalent(bool f = true)
  {
    restrictEquivalent = f;
    log() << "Forcing tolopogical equivalency: " << f << std::endl;
  }

  void setRecalcTopology(bool f = true)
  {
    recalcTopology = f;
    log() << "Tolopogical information will be recalculated: " << f << std::endl;
  }
  
  void setRefit(bool f = true)
  {
    refit = f;
    log() << "Multipole refit requested: " << f << std::endl;
  }
  
  void setMolecule(GeneralMultipoledMolecule *m)
  {
    multMol = m;
    mol = GeneralMultipoledMolecule(*multMol); // we need a copy to safely add dummy atoms  
    numRealAtoms = multMol->NumAtoms();
    log() << "Molecule with " << numRealAtoms << " atoms loaded." << std::endl;
  }

  void setMultipolePlacementRules(const MultipolePlacementRules& r)
  {
    p_rules = r;
    log() << p_rules.numAtomRules() << " atom and " << p_rules.numGroupRules() << " group multipole positioning rules loaded." << std::endl;
  }

  void setGrid(const Grid& g)
  {
    grid = g;
    log() << grid.size() << " points of field loaded." << std::endl;
  }
  
  double rmsd() const { return RSMD; }

  void addDummyCenters()
  {
    // creating dummies for groups
    OpenBabel::OBSmartsPattern patt;
    paramFlags.resize(numRealAtoms);
    for (int i = 0; i < p_rules.numGroupRules(); ++i)
    {
      MultipolePlacementRule mult_rule = p_rules.groupRule(i);
      patt.Init(mult_rule.SMARTS);
      patt.Match(mol);
      std::vector< std::vector< int > > hits = patt.GetUMapList();
      std::vector< int > addedDummies;
      for (int j = 0; j < hits.size(); ++j)
      {
        Vector3 dummyCrd(0.0, 0.0, 0.0);
        for (int k = 0; k < hits[j].size(); ++k)
        {
          dummyCrd += Vector3(multMol->GetAtom(hits[j][k])->GetCoordinate());
        }
        dummyCrd /= hits[j].size();
        OpenBabel::OBAtom *a = multMol->GetAtom(hits[j][0]);
        OpenBabel::OBAtom *b = multMol->GetAtom(hits[j][1]);
        multMol->BeginModify();
        OpenBabel::OBAtom *c = multMol->NewAtom();
        c->SetAtomicNum(99);
        c->SetVector(dummyCrd[0], dummyCrd[1], dummyCrd[2]);
        OpenBabel::OBBond *b1 = multMol->NewBond();
        OpenBabel::OBBond *b2 = multMol->NewBond();
        b1->SetBondOrder(1); b2->SetBondOrder(1);
        b1->SetBegin(a); b1->SetEnd(c);
        b2->SetBegin(b); b2->SetEnd(c);
        a->AddBond(b1); c->AddBond(b1);
        b->AddBond(b2); c->AddBond(b2);
        multMol->EndModify();
        addedDummies.push_back(c->GetIdx());
        // need to save parameters flags
        paramFlags.push_back(mult_rule.rule);
      }
      // checking equivalent groups
      if (hits.size() > 1 && restrictEquivalent)
      {
        int hitsize = hits[0].size();
        for (int k = 0; k < hits.size(); ++k)
        {
          for (int l = k + 1; l < hits.size(); ++l)
          {
            std::vector< double > ch1, ch2;
            for (int m = 0; m < hitsize; ++m)
            {
              ch1.push_back(mol.GetAtom(hits[k][m])->GetPartialCharge());
              ch2.push_back(mol.GetAtom(hits[l][m])->GetPartialCharge());
            }
            std::sort(ch1.begin(), ch1.end());
            std::sort(ch2.begin(), ch2.end());
            if (std::equal(ch1.begin(), ch1.end(), ch2.begin()))
            {
              equal_centers.push_back(std::pair< int, int >(addedDummies[k], addedDummies[l]));
            }
          }
        }
      }
    }
    for (int i = 1; i < mol.NumAtoms() + 1; ++i) if (mol.GetAtom(i)->IsAromatic()) multMol->GetAtom(i)->SetAromatic();
    for (int i = 0; i < mol.NumBonds(); ++i) if (mol.GetBond(i)->IsAromatic()) multMol->GetBond(i)->SetAromatic();
    multMol->SetAromaticPerceived();
    log() << multMol->NumAtoms() - mol.NumAtoms() << " dummy atoms added." << std::endl;
  }

  void placeAtomParameters(const std::vector<std::string>& ignore_masks = std::vector<std::string>())
  {
    OpenBabel::OBSmartsPattern patt;
    for (unsigned i = 0; i < p_rules.numAtomRules(); ++i)
    {
      MultipolePlacementRule mult_rule = p_rules.atomRule(i);
      patt.Init(mult_rule.SMARTS);
      patt.Match(mol);
      std::vector< std::vector< int > > hits = patt.GetMapList();
      for (unsigned j = 0; j < hits.size(); ++j)
        paramFlags[hits[j][0] - 1] = mult_rule.rule;
    }

    for (unsigned i = 0; i < ignore_masks.size(); ++i)
    {
      patt.Init(ignore_masks[i]);
      patt.Match(mol);
      std::vector< std::vector< int > > hits = patt.GetMapList();
      for (unsigned j = 0; j < hits.size(); ++j)
      {
        bool monopole = paramFlags[hits[j][0] - 1].test(MONOPOLE);
        paramFlags[hits[j][0] - 1].clear();
        if (monopole) paramFlags[hits[j][0] - 1].set(MONOPOLE);
      }
    }

    if ( (!restrictEquivalent) && recalcTopology)
      log() << "Recalculation of topology requested but topology information is ignored. Skipping this step." << std::endl;

    if (restrictEquivalent && recalcTopology)
    {
      OpenBabel::OBChargeModel *gasteigerCharges = OpenBabel::OBChargeModel::FindType("gasteiger");
      if (gasteigerCharges)
      {
        log() << "Recalculating topology..." << std::endl;
        std::vector<double> chrgs; chrgs.reserve(mol.NumAtoms()+1); chrgs.push_back(0.0);
        for (int i = 1 ; i <= mol.NumAtoms(); ++i)
        {
          chrgs.push_back(mol.GetAtom(i)->GetPartialCharge());
          mol.GetAtom(i)->SetPartialCharge(0.0);
        }
        if (!gasteigerCharges->ComputeCharges(mol))
        {
          log() << "Failure recalculting topology. Using provided one." << std::endl;
          for (int i = 1 ; i <= mol.NumAtoms(); ++i) mol.GetAtom(i)->SetPartialCharge(chrgs[i]);
        }
      }
    }
    if (restrictEquivalent)
    {
      std::vector< std::pair< int, double > > charges;
      charges.reserve(mol.NumAtoms());
      for (int i = 1; i < mol.NumAtoms() + 1; ++i)
        charges.push_back(std::pair< int, double >(i, mol.GetAtom(i)->GetPartialCharge()));

      std::sort(charges.begin(), charges.end(), compare_pair_second());
      for (int i = 1; i < charges.size(); ++i)
        if (charges[i - 1].second == charges[i].second)
          equal_centers.push_back(std::pair< int, int >(charges[i - 1].first, charges[i].first));
    }
  }

  void createMatrix()
  {
    has_charges = 0;
    multMol->createOrientMap();
    multMol->calcOrientMatrices();
    //create parameteer map and adding symmetry constraints
    e_constraints.resize(0);
    q_constraints.resize(0);
    z_constraints.resize(0);
    d_constraints.resize(0);
    atomBegins.resize(multMol->NumAtoms() + 1);
    numParams = 0;
    for (unsigned i = 1; i < multMol->NumAtoms() + 1; ++i)
    {
      MultipoleFlags flags = paramFlags[i - 1];
      atomBegins[i] = ParamMarkers(-1, -1, -1);
      if (flags.test(NONE))
      {
        atomBegins[i] = ParamMarkers(-1, -1, -1);
      } else
      {
        if (flags.test(MONOPOLE)) { atomBegins[i][0] = numParams++; has_charges=1; }
        if (flags.test(DIPOLE)) 
        { 
          atomBegins[i][1] = numParams++; numParams++; numParams++; // for dx, dy, dz
          if (flags.test(DIPOLE_DIRECTED)) 
          {
            Vector3 dir = flags.dipoleDirection();
            bool nonzero_x = (fabs(dir.x()) > 1e-7); 
            bool nonzero_y = (fabs(dir.y()) > 1e-7); 
            bool nonzero_z = (fabs(dir.z()) > 1e-7); 
            if ( !nonzero_x ) z_constraints.push_back(numParams -3); // if dipole direction is zero zero out this component
            if ( !nonzero_y ) z_constraints.push_back(numParams -2);
            if ( !nonzero_z ) z_constraints.push_back(numParams -1);
            if ( nonzero_x && nonzero_y && nonzero_z )  // a little bit overcomplicated logics to avoid excessive constraint
            {                                           // (x/x0 == y/y0, y/y0=z/z0, and x/x0 = z/z0)
              d_constraints.push_back(boost::tuple<int,int,double,double>(numParams-3,numParams-2,1.0/dir.x(),-1.0/dir.y()));
              d_constraints.push_back(boost::tuple<int,int,double,double>(numParams-2,numParams-1,1.0/dir.y(),-1.0/dir.z()));
            } 
            else if ( nonzero_x && nonzero_y ) 
              d_constraints.push_back(boost::tuple<int,int,double,double>(numParams-3,numParams-2,1.0/dir.x(),-1.0/dir.y()));
            else if ( nonzero_y && nonzero_z ) 
              d_constraints.push_back(boost::tuple<int,int,double,double>(numParams-2,numParams-1,1.0/dir.y(),-1.0/dir.z()));
            else if ( nonzero_x && nonzero_z ) 
              d_constraints.push_back(boost::tuple<int,int,double,double>(numParams-3,numParams-1,1.0/dir.x(),-1.0/dir.z()));
          }
        }
        if (flags.test(QUADRUPOLE)) // as with dipole for consistency reasons 
                                    // we fit all Qxx, Qyy and Qzz with relevant constraints
        { 
          atomBegins[i][2] = numParams++; numParams++; numParams++; 
          q_constraints.push_back(boost::tuple<int,int,int>(numParams-3, numParams-2, numParams-1));
          if (flags.test(QUADRUPOLE_X)) { e_constraints.push_back(std::pair<int,int>(numParams -2, numParams -1)); } else
          if (flags.test(QUADRUPOLE_Y)) { e_constraints.push_back(std::pair<int,int>(numParams -3, numParams -1)); } else
          if (flags.test(QUADRUPOLE_Z)) { e_constraints.push_back(std::pair<int,int>(numParams -3, numParams -2)); }
        }
      }
    }
    //adding equal centers constraints
    if (restrictEquivalent)
      for (unsigned i = 0; i < equal_centers.size(); ++i)
      {
        int a = equal_centers[i].first, b = equal_centers[i].second;
        int k, l;
        if (( k = atomBegins[a][0] ) != -1 && ( l = atomBegins[b][0] ) != -1) // constraints for charges
          e_constraints.push_back(std::pair< int, int >(k, l));
        if (( k = atomBegins[a][1] ) != -1 && ( l = atomBegins[b][1] ) != -1) // constraints for dipoles
        {
          MultipoleFlags flags = paramFlags[a - 1];
          if (!flags.test(DIPOLE_DIRECTED))  // equlize everything if not directed
          {
            e_constraints.push_back(std::pair< int, int >(k, l));
            e_constraints.push_back(std::pair< int, int >(k+1, l+2));
            e_constraints.push_back(std::pair< int, int >(k+1, l+2));
          } else {                           // equlize single nonzero parameter otherwise
            Vector3 dir = flags.dipoleDirection();
            if ( fabs(dir.x()) > 1e-7 ) e_constraints.push_back(std::pair< int, int >(k, l)); else
            if ( fabs(dir.y()) > 1e-7 ) e_constraints.push_back(std::pair< int, int >(k+1, l+1)); else
            if ( fabs(dir.z()) > 1e-7 ) e_constraints.push_back(std::pair< int, int >(k+2, l+2));
          }
        }
        if (( k = atomBegins[a][2] ) != -1 && ( l = atomBegins[b][2] ) != -1) // constraints for quadrupoles
        {
          MultipoleFlags flags = paramFlags[a - 1];
          if (flags.test(QUADRUPOLE_X)) { e_constraints.push_back(std::pair<int,int>(k,l)); } else
          if (flags.test(QUADRUPOLE_Y)) { e_constraints.push_back(std::pair<int,int>(k+1,l+1)); } else
          if (flags.test(QUADRUPOLE_Z)) { e_constraints.push_back(std::pair<int,int>(k+2,l+2)); } 
          else { e_constraints.push_back(std::pair<int,int>(k,l)); e_constraints.push_back(std::pair<int,int>(k+1,l+1)); }
        }
      }
    //creating distance matrix
    A.resize(grid.size(), numParams);
    for (unsigned i = 1; i < multMol->NumAtoms() + 1; ++i)
    {
      Matrix3 Orient = multMol->matrix(i);
      int col = -1;
      Vector3 acrd(multMol->GetAtom(i)->GetCoordinate());
      acrd /= bohr;
      MultipoleFlags flags = paramFlags[i - 1];
      // filling for monopole
      if (( col = atomBegins[i][0] ) != -1)
        for (unsigned j = 0; j < grid.size(); ++j)
          A(j, col) = 1.0 / ( grid.lattice(j) - acrd ).length();
      // filling for dipole
      if (( col = atomBegins[i][1] ) != -1)
      {
        for (unsigned j = 0; j < grid.size(); ++j)
        {
          Vector3 r = ( grid.lattice(j) - acrd ) * Orient;
          double r_l = r.length();
          double r_l3 = r_l*r_l*r_l;
          A(j, col) = r[0] / r_l3;
          A(j, col+1) = r[1] / r_l3;
          A(j, col+2) = r[2] / r_l3;
        }
      }
      // filling for quarupole
      if (( col = atomBegins[i][2] ) != -1)
      {
        for (unsigned j = 0; j < grid.size(); ++j)
        {
          Vector3 r = ( grid.lattice(j) - acrd ) * Orient;
          double r_l = r.length();
          double r_l5 = r_l*r_l*r_l*r_l*r_l;
          A(j, col) =  r[0] * r[0] / r_l5;
          A(j, col+1) =  r[1] * r[1] / r_l5;
          A(j, col+2) =  r[2] * r[2] / r_l5;
        }
      }
    }
    log() << "Created distance matrix " << grid.size() << "x" << numParams << std::endl;
    //creating field matrix
    f = Eigen::VectorXd(grid.size());
    for (unsigned i = 0; i < grid.size(); ++i) f(i) = grid.field(i);
    if (refit) // removing reference multipoles from the field
    {
      for (unsigned i = 1; i < multMol->NumAtoms() + 1; ++i)
      {
        MultipoleFlags flags = paramFlags[i - 1];
        Multipoles* mult = static_cast<Multipoles*>(multMol->GetAtom(i)->GetData("Multipoles"));
        if (!mult) 
        {
          mult = new Multipoles(Vector3(multMol->GetAtom(i)->GetCoordinate()));
          mult->setMonopole(multMol->GetAtom(i)->GetPartialCharge());
        }
        Multipoles tmult(Vector3(multMol->GetAtom(i)->GetCoordinate()));

        if ( ( ! flags.test(MONOPOLE) ) && mult->hasMonopole() ) tmult.setMonopole(mult->monopole());
        if ( ( ! flags.test(DIPOLE) ) && mult->hasDipole() ) tmult.setDipole(mult->dipole()); 
        if ( ( ! flags.test(QUADRUPOLE) ) && mult->hasQuadrupole() ) tmult.setQuadrupole(mult->quadrupole());
        for (int j = 0; j < grid.size(); ++j) f(j) -= tmult.potential(grid.lattice(j)*bohr);
      }

    }
    // calculating q_tot
    q_tot = grid.charge();
    if (refit)
    {
      for (unsigned i = 1; i < multMol->NumAtoms() + 1; ++i)
        if ( atomBegins[i][0] == -1 ) q_tot -= multMol->GetAtom(i)->GetPartialCharge();
    }
  }

  virtual void makeFit() = 0;

  void getResults()
  {
    for (unsigned i = 1; i < multMol->NumAtoms() + 1; ++i)
    {
      MultipoleFlags flags = paramFlags[i - 1];
      Multipoles *mult;
      if (refit)
      {
        mult = static_cast<Multipoles*>(multMol->GetAtom(i)->GetData("Multipoles"));
        if (!mult)
        {
          mult = new Multipoles(Vector3(multMol->GetAtom(i)->GetCoordinate()));
          mult->setMonopole(multMol->GetAtom(i)->GetPartialCharge());
          multMol->GetAtom(i)->SetData(mult);
        }
      } else {
        multMol->GetAtom(i)->DeleteData("Multipoles");
        mult = new Multipoles(Vector3(multMol->GetAtom(i)->GetCoordinate()));
        mult->setMonopole(multMol->GetAtom(i)->GetPartialCharge());
        multMol->GetAtom(i)->SetData(mult);
      }
      int col = -1;
      if (( col = atomBegins[i][0] ) != -1) mult->setMonopole(Monopole(q(col)));
      if (( col = atomBegins[i][1] ) != -1) mult->setDipole(Dipole(q(col), q(col+1), q(col+2)));
      if (( col = atomBegins[i][2] ) != -1) mult->setQuadrupole(Quadrupole(q(col), q(col + 1), q(col+2)));
      multMol->GetAtom(i)->SetPartialCharge(mult->monopole().value());
    }
    multMol->SetPartialChargesPerceived();
    log() << "Fit completed." << std::endl;
  }

  virtual void printResult(std::ostream& sink) const = 0;

  void fit()
  {
    addDummyCenters();
    placeAtomParameters();
    createMatrix();
    makeFit();
    getResults();
  }
  
  virtual ~MEPFitter() { }

private:

  struct ParamMarkers
  {
    int data[3];

    ParamMarkers(int a = -1, int b = -1, int c = -1)
    {
      data[0] = a;
      data[1] = b;
      data[2] = c;
    }

    int& operator[](int i) { return data[i]; }

    int operator[](int i) const { return data[i]; }
  };

  struct compare_pair_second
  {
    bool operator()(const std::pair< int, double >& l, const std::pair< int, double >& r) { return l.second > r.second; }
  };

protected:
  Grid grid;
  Eigen::MatrixXd A;  // distance matrix
  Eigen::MatrixXd f; // field matrix
  Eigen::VectorXd q; // result vector
  GeneralMultipoledMolecule mol;
  GeneralMultipoledMolecule *multMol;
  int numRealAtoms;
  std::vector< ParamMarkers > atomBegins;
  std::vector< bool > symmFlags;
  std::vector< std::pair< int, int > > equal_centers;
  std::vector< MultipoleFlags > paramFlags;
  MultipolePlacementRules p_rules;
  MultipoleOrientRules o_rules;
  int numParams;
  bool restrictEquivalent;
  bool recalcTopology;
  bool refit;
  double q_tot;
  std::vector< std::pair< int, int > > e_constraints; // constraints to equilize parameters
  std::vector< boost::tuple< int, int, int> > q_constraints; // constraints to make valid quadrupoles
  std::vector< int > z_constraints; // constraints to zero-out parameters
  std::vector< boost::tuple< int, int, double, double> > d_constraints; // constraints to align dipoles
  double RSMD;
  int has_charges;
};

class MEPSVDFitter : public MEPFitter
{
public:
  MEPSVDFitter(std::ostream& log) : MEPFitter(log), cutoff_criterion(0.0) { }

  void setCutoff(const double c) { cutoff_criterion = c; }

  virtual void makeFit()
  {
    log() << "Using SVD fitter." << std::endl;
    int n_constraints = e_constraints.size() + z_constraints.size() + q_constraints.size() + d_constraints.size() + has_charges; 
    B = Eigen::MatrixXd::Zero(n_constraints, numParams); 
    C = Eigen::VectorXd::Zero(n_constraints);
    unsigned c_offset = 0;
    for (unsigned i = 0; i < e_constraints.size(); ++i) // constraints for equal parameters
    {
      B(c_offset+i, e_constraints[i].first) = 1.0;
      B(c_offset+i, e_constraints[i].second) = -1.0;
    }
    c_offset += e_constraints.size();
    for (unsigned i = 0; i < z_constraints.size(); ++i) // constraints for zero parameters
    {
      B(c_offset+i, z_constraints[i]) = 1.0;
    }
    c_offset += z_constraints.size();
    for (unsigned i = 0; i < d_constraints.size(); ++i) // constraints for dipoles
    {
      B(c_offset+i, d_constraints[i].get<0>()) = d_constraints[i].get<2>();
      B(c_offset+i, d_constraints[i].get<1>()) = d_constraints[i].get<3>();
    }
    c_offset += d_constraints.size();
    for (unsigned i = 0; i < q_constraints.size(); ++i) // constraints for zero-trace quadrupoles
    {
      B(c_offset+i, q_constraints[i].get<0>()) = 1.0;
      B(c_offset+i, q_constraints[i].get<1>()) = 1.0;
      B(c_offset+i, q_constraints[i].get<2>()) = 1.0;
    }
    c_offset += q_constraints.size();
    if (has_charges)
    {
      for (unsigned i = 1; i < atomBegins.size(); ++i) // q_tot constraint
      {
        int k = -1;
        if (( k = atomBegins[i][0] ) != -1) B(c_offset, k) = 1.0; 
      }
      C(c_offset) = q_tot;  // q_tot constrained
    }
    log() << "Created constraints matrix " << ( n_constraints ) << "x" << numParams << std::endl;
    log() << "Fitting" << std::endl;
    // now ready for pseudoinversion
    Eigen::JacobiSVD< Eigen::MatrixXd > Bsvd(B, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd Bsv = Bsvd.singularValues();
    Eigen::MatrixXd Bsigma = Eigen::MatrixXd::Zero(n_constraints, n_constraints);
    for (int i = 0; i < n_constraints; ++i) Bsigma(i, i) = ( Bsv(0) / Bsv(i) < cutoff_criterion ) ? 1.0 / Bsv(i) : 0.0;
    Eigen::MatrixXd Binv = Bsvd.matrixV() * Bsigma * Bsvd.matrixU().transpose();
    Eigen::MatrixXd E_minus_BB = Eigen::MatrixXd::Identity(B.cols(), B.cols()) - ( Binv * B );
    Eigen::MatrixXd Anew = A * E_minus_BB;
    Eigen::JacobiSVD< Eigen::MatrixXd > Anewsvd(Anew, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd Anewsv = Anewsvd.singularValues();
    Eigen::MatrixXd Anewsigma = Eigen::MatrixXd::Zero(Anewsv.size(), Anewsv.size());
    for (int i = 0; i < Anewsv.size(); ++i) Anewsigma(i, i) = ( Anewsv(0) / Anewsv(i) < cutoff_criterion ) ? 1.0 / Anewsv(i) : 0.0;
    Eigen::MatrixXd Anewinv = Anewsvd.matrixV() * Anewsigma * Anewsvd.matrixU().transpose();
    q = ( Binv * C ) + ( E_minus_BB * Anewinv * ( f - ( ( A * Binv ) * C ) ) );
    // now q has found parameters
    // qnew to assign it to molecule and print results;
    RSMD = sqrt(1.0 / grid.size() * ( ( A * q ) - f ).squaredNorm());
    singular_values = Anewsvd.singularValues();
    VT = Anewsvd.matrixV().transpose();
  }

  virtual void printResult(std::ostream& sink) const
  { 
    sink << *multMol << std::endl;
    sink << "/*" << std::endl;
    sink << "Fitted " << numParams << " parameters with " << B.rows() << " constraints (" 
                      << grid.size() << " reference points)" << std::endl;
    sink << "RMSD: " << ( RSMD * 629.5 ) << " kcal/mol" << std::endl;
    sink << "*/" << std::endl;

    log() << "Fitted " << numParams << " parameters with " << B.rows()<< " constraints ("
                      << grid.size() << " reference points)" << std::endl;
    log() << "RMSD: " << ( RSMD * 629.5 ) << " kcal/mol" << std::endl;
    log() << "----------------- Constraints -----------------" << std::endl;
    std::ostringstream paramstr;
    paramstr << "Parameters: ";
    for (unsigned i = 1; i < multMol->NumAtoms() + 1; ++i)
    {
      if (atomBegins[i][0] != -1) paramstr << "Z" << i << "\t";
      if (atomBegins[i][1] != -1) paramstr << "Dx" << i << "\t" << "Dy" << i << "\t" << "Dz" << i << "\t";
      if (atomBegins[i][2] != -1) paramstr << "Qxx" << i << "\t" << "Qyy" << i << "\t" << "Qzz" << i << "\t";
    }
    paramstr << std::endl;
    log() << paramstr.str();
    log() << "Constraints matrix (B): " << std::endl << B << std::endl;
    log() << "Constraints vecctor (C): " << std::endl << C << std::endl;
    log() << "----------------- SVD analysis (after constraints applied) -----------------" << std::endl;
    log() << "Condition number cutoff: " << cutoff_criterion << std::endl;
    log() << "Condition number of the system: " << (singular_values(0)/singular_values(numParams-B.rows()-1)) << std::endl;
    log() << "Singular values:" << std::endl << singular_values << std::endl;
    log() << paramstr.str();
    log() << "Right singular vectors (V**T): " << std::endl << VT << std::endl;
  }

private:
  Eigen::MatrixXd B;  // constraints pt1
  Eigen::VectorXd C;  // constraints pt2
  Eigen::VectorXd singular_values;
  Eigen::MatrixXd VT;
  double cutoff_criterion;
};

class MEPLagrangeFitter : public MEPFitter
{
public:
  MEPLagrangeFitter(std::ostream& log) : MEPFitter(log) { }

  virtual void makeFit()
  {
    // performs simple fit with Lagrange constraints
    // creating least squares matrix
    log() << "Using simple fitter." << std::endl;
    int n_constraints = e_constraints.size() + z_constraints.size() + q_constraints.size() + d_constraints.size() + has_charges;
    Eigen::MatrixXd LSMat;
    LSMat=Eigen::MatrixXd(numParams + n_constraints, numParams + n_constraints);
    LSMat.topLeftCorner(A.cols(), A.cols()) = A.transpose() * A;
    LSMat.rightCols(LSMat.cols() - A.cols()) = Eigen::MatrixXd::Zero(LSMat.rows(), LSMat.cols() - A.cols());
    LSMat.bottomLeftCorner(LSMat.rows() - A.cols(), A.cols()) = Eigen::MatrixXd::Zero(LSMat.rows() - A.cols(), A.cols());
    Eigen::VectorXd b(numParams + n_constraints);
    b.head(A.cols()) = A.transpose() * f;
    b.tail(b.rows() - A.cols()) = Eigen::VectorXd::Zero(b.rows() - A.cols());

    unsigned c_offset = numParams;
    for (unsigned i = 0; i < e_constraints.size(); ++i) // constraints for equal parameters
    {
      LSMat(c_offset + i, e_constraints[i].first) = 1.0;
      LSMat(c_offset + i, e_constraints[i].second) = -1.0;
      LSMat(e_constraints[i].first, c_offset + i) = 1.0;
      LSMat(e_constraints[i].second, c_offset + i) = -1.0;
    }
    c_offset += e_constraints.size();
    for (unsigned i = 0; i < z_constraints.size(); ++i) // constraints for zero parameters
    {
      LSMat(c_offset + i, z_constraints[i]) = 1.0;
      LSMat(z_constraints[i], c_offset + i) = 1.0;
    }
    c_offset += z_constraints.size();
    for (unsigned i = 0; i < d_constraints.size(); ++i) // constraints for dipoles
    {
      LSMat(c_offset + i, d_constraints[i].get<0>()) = d_constraints[i].get<2>();
      LSMat(c_offset + i, d_constraints[i].get<1>()) = d_constraints[i].get<3>();
      LSMat(d_constraints[i].get<0>(), c_offset + i) = d_constraints[i].get<2>();
      LSMat(d_constraints[i].get<1>(), c_offset + i) = d_constraints[i].get<3>();
    }
    c_offset += d_constraints.size();
    for (unsigned i = 0; i < q_constraints.size(); ++i) // constraints for zero-trace quadrupoles
    {
      LSMat(c_offset + i, q_constraints[i].get<0>()) = 1.0;
      LSMat(c_offset + i, q_constraints[i].get<1>()) = 1.0;
      LSMat(c_offset + i, q_constraints[i].get<2>()) = 1.0;
      LSMat(q_constraints[i].get<0>(), c_offset + i) = 1.0;
      LSMat(q_constraints[i].get<1>(), c_offset + i) = 1.0;
      LSMat(q_constraints[i].get<2>(), c_offset + i) = 1.0;
    }
    c_offset += q_constraints.size();
    if (has_charges)
    {
      for (unsigned i = 1; i < atomBegins.size(); ++i)
      {
        int k = -1;
        if (( k = atomBegins[i][0] ) != -1)
        {
          LSMat(c_offset, k) = 1.0;
          LSMat(k, c_offset) = 1.0;
        }
      }
      b(c_offset) = q_tot;
    }
    log() << "Matrices created." << std::endl;
    //performing fit
    q = solve(LSMat, b);
    Eigen::VectorXd q_rmsd = q.head(A.cols());
    RSMD = sqrt(1.0 / grid.size() * ( ( A * q_rmsd ) - f ).squaredNorm());
  }

  virtual void printResult(std::ostream& sink) const
  {
    int n_constraints = e_constraints.size() + z_constraints.size() + q_constraints.size() + d_constraints.size() + has_charges;
    sink << *multMol << std::endl;
    sink << "/*" << std::endl;
    sink << "Fitted " << numParams << " parameters with " << n_constraints 
         << " constraints  (" << grid.size() << " reference points)" << std::endl;
    sink << "RMSD : " << ( RSMD * 629.5 ) << " kcal/mol" << std::endl;
    sink << "*/" << std::endl;

    log() << "Fitted " << numParams << " parameters with " << n_constraints 
          << " constraints  (" << grid.size() << " reference points)" << std::endl;
    log() << "RMSD : " << ( RSMD * 629.5 ) << " kcal/mol" << std::endl;
  }

protected:
  virtual Eigen::VectorXd solve(const Eigen::MatrixXd& m, const Eigen::VectorXd& b) = 0;
};

class MEPLagrangeLDLTFitter : public MEPLagrangeFitter
{
public:
  MEPLagrangeLDLTFitter(std::ostream& log) : MEPLagrangeFitter(log) { }

protected:
  virtual Eigen::VectorXd solve(const Eigen::MatrixXd& m, const Eigen::VectorXd& b)
  {
    log() << "Using LDLT solver." << std::endl;
    return m.ldlt().solve(b);
  }
};

class MEPLagrangeFullPivLUFitter : public MEPLagrangeFitter
{
public:
  MEPLagrangeFullPivLUFitter(std::ostream& log) : MEPLagrangeFitter(log) { }

protected:
  virtual Eigen::VectorXd solve(const Eigen::MatrixXd& m, const Eigen::VectorXd& b) 
  {
    log() << "Using Complete Pivoting LU solver." << std::endl;
    Eigen::FullPivLU<Eigen::MatrixXd> lu(m);
    return lu.solve(b);
  }
};

class MEPLagrangePartialPivLUFitter : public MEPLagrangeFitter
{
public:
  MEPLagrangePartialPivLUFitter(std::ostream& log) : MEPLagrangeFitter(log) { }

protected:
  virtual Eigen::VectorXd solve(const Eigen::MatrixXd& m, const Eigen::VectorXd& b) 
  {
    log() << "Using Partial Pivoting LU solver." << std::endl;
    Eigen::PartialPivLU<Eigen::MatrixXd> lu(m);
    return lu.solve(b);
  }
};

class MEPLagrangeHouseholderQRFitter : public MEPLagrangeFitter
{
public:
  MEPLagrangeHouseholderQRFitter(std::ostream& log) : MEPLagrangeFitter(log) { }

protected:
  virtual Eigen::VectorXd solve(const Eigen::MatrixXd& m, const Eigen::VectorXd& b) 
  {
    log() << "Using Householder QR solver." << std::endl;
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(m);
    return qr.solve(b);
  }
};

class MEPLagrangeColPivHouseholderQRFitter : public MEPLagrangeFitter
{
public:
  MEPLagrangeColPivHouseholderQRFitter(std::ostream& log) : MEPLagrangeFitter(log) { }

protected:
  virtual Eigen::VectorXd solve(const Eigen::MatrixXd& m, const Eigen::VectorXd& b)
  {
    log() << "Using Column Pivoting Householder QR solver." << std::endl;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(m);
    return qr.solve(b);
  }
};

class MEPLagrangeFullPivHouseholderQRFitter : public MEPLagrangeFitter
{
public:
  MEPLagrangeFullPivHouseholderQRFitter(std::ostream& log) : MEPLagrangeFitter(log) { }

protected:
  virtual Eigen::VectorXd solve(const Eigen::MatrixXd& m, const Eigen::VectorXd& b)
  {
    log() << "Using Full Pivoting Householder QR solver." << std::endl;
    Eigen::FullPivHouseholderQR<Eigen::MatrixXd> qr(m);
    return qr.solve(b);
  }
};


class MEPLagrangeLLTFitter : public MEPLagrangeFitter
{
public:
  MEPLagrangeLLTFitter(std::ostream& log) : MEPLagrangeFitter(log) { }

protected:
  virtual Eigen::VectorXd solve(const Eigen::MatrixXd& m, const Eigen::VectorXd& b)
  {
    log() << "Using LLT solver." << std::endl;
    return m.llt().solve(b);
  }
};

class FitterFactory
{
public:
  FitterFactory(std::ostream& log = std::clog) : log_(log) { }

  MEPFitter *getFitter(const std::string key)
  {
    if (!key.compare("SVD")) return new MEPSVDFitter(log_);
    if (!key.compare("PartialPivLU")) return new MEPLagrangePartialPivLUFitter(log_);
    if (!key.compare("FullPivLU")) return new MEPLagrangeFullPivLUFitter(log_);
    if (!key.compare("HouseholderQR")) return new MEPLagrangeHouseholderQRFitter(log_);
    if (!key.compare("ColPivHouseholderQR")) return new MEPLagrangeColPivHouseholderQRFitter(log_);
    if (!key.compare("FullPivHouseholderQR")) return new MEPLagrangeFullPivHouseholderQRFitter(log_);
    if (!key.compare("LLT")) return new MEPLagrangeLLTFitter(log_);
    if (!key.compare("LDLT")) return new MEPLagrangeLDLTFitter(log_);
    throw std::invalid_argument("Unknown fitter type: " + key);
  }

private:
  std::ostream& log_;
};

#endif

