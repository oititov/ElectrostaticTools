#ifndef BCC_FITTER_HEADER
#define BCC_FITTER_HEADER

#include <vector>
#include <string>
#include <stdexcept>

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <Eigen/QR>


#include "molecule.h"
#include "grid.h"
#include "par_table.h"
#include "bcc_typer.h"

class BCCFitter
{
public: 
  BCCFitter(ParameterTable& table, std::vector<GeneralMultipoledMolecule>& vmols, std::vector<Grid>& vgrids,
            std::vector<std::vector<std::string> >& vtypes) 
          : partable(table), mols(vmols), grids(vgrids), types(vtypes) 
  { 
    if (mols.size() != grids.size()) throw std::domain_error("Number of molecules does not correspond to the number of grids."); 
  }
  
  void makeFit()
  {
    param_map.clear();
    createParamMap();
    activeMols.clear();
    findActiveMols();
    int gridLength = calcGridLength();
    int numParams = param_map.size();
    A = Eigen::MatrixXd::Zero(gridLength, numParams);
    C = Eigen::VectorXd::Zero(gridLength);
    fillMatrices();
    R = solve(A, C);
    fillParameters();
  }
  std::vector<double> rmsds()
  {
    std::vector<double> result;
    double global_rmsd = 0;
    double global_before = 0;
    result.reserve(mols.size()*2 + 2);
    int row_offset = 0;
    Eigen::VectorXd E = C - (A*R);
    for (int i = 0, j = 0; j < activeMols.size(); ++i, ++j)
    {
      while (activeMols[j] != i) 
      {
        result.push_back(99999);
        result.push_back(99999);
        ++i;
      }
      double rmsd = C.segment(row_offset, grids[i].size()).squaredNorm();
      global_before += rmsd;
      result.push_back(sqrt(rmsd/grids[i].size())); 
      rmsd = E.segment(row_offset, grids[i].size()).squaredNorm();
      global_rmsd += rmsd;
      result.push_back(sqrt(rmsd/grids[i].size())); 
      row_offset += grids[i].size();     
    }
    result.push_back(sqrt(global_before/E.size()));
    result.push_back(sqrt(global_rmsd/E.size()));
    return result;
  }
private:
  void createParamMap() 
  {
    int n = 0;
    for (int i = 0; i < mols.size(); ++i)
    {
      BCCAtomTyper atyper(types[i]);
      BCCBondTyper btyper(atyper);
    
      for (OpenBabel::OBBondIterator bond = mols[i].BeginBonds(), ebond = mols[i].EndBonds(); bond != ebond; ++bond)
      {
        std::string bond_type = btyper.type(*bond);
        if ( partable.exists(bond_type, "bcc") && (param_map.find(bond_type) == param_map.end()))
        {
          param_map[bond_type] = n++;
        }
      }  
    }
  }
//  bool isMolFitActive(OpenBabel::OBMol* mol)
  bool isMolFitActive(OpenBabel::OBMol* mol, int molidx)
  {
    BCCAtomTyper atyper(types[molidx]);
    BCCBondTyper btyper(atyper);
    
    for (OpenBabel::OBBondIterator bond = mol->BeginBonds(), ebond = mol->EndBonds(); bond != ebond; ++bond)
    {
      std::string bond_type = btyper.type(*bond);
      if (partable.exists(bond_type, "bcc"))
      {
        return true;
      }
    }
    return false;
  }
  void findActiveMols()
  {
//    for (int i = 0; i < mols.size(); ++i) if (isMolFitActive(&(mols[i]))) activeMols.push_back(i);
    for (int i = 0; i < mols.size(); ++i) if (isMolFitActive(&(mols[i]),i)) activeMols.push_back(i);
  }
  int calcGridLength()
  {
    int result = 0;
//    for (int i = 0; i < mols.size(); ++i) if (isMolFitActive(&(mols[i]))) result += grids[i].size();
    for (int i = 0; i < mols.size(); ++i) if (isMolFitActive(&(mols[i]),i)) result += grids[i].size();
    return result;
  }
  
  void fillMatrices()
  {
    unsigned row_offset = 0;
    for (int i = 0; i < activeMols.size(); ++i)
    {
      int j = activeMols[i];
      Grid grid = grids[j];
      // remove potential from existing charges
      
      for (int k = 0; k < grid.size(); ++k) grid.modifyField(k, -mols[j].MEP(grid.lattice(k)*bohr)); 
 
      BCCAtomTyper atyper(types[j]);
      BCCBondTyper btyper(atyper);

      for (OpenBabel::OBBondIterator bond = mols[j].BeginBonds(), ebond = mols[j].EndBonds(); bond != ebond; ++bond)
      {
        std::string bond_type = btyper.type(*bond);
        if (partable.exists(bond_type, "bcc"))
        {
          partable.mark(bond_type, "bcc");
          int col_pos = param_map[bond_type];
          OpenBabel::OBAtom* pa = (*bond)->GetBeginAtom();
          OpenBabel::OBAtom* pb = (*bond)->GetEndAtom();
          std::string atype = atyper.type(pa);
          std::string btype = atyper.type(pb);
          Multipoles* am = static_cast<Multipoles*>(pa->GetData("Multipoles"));
          Multipoles* bm = static_cast<Multipoles*>(pb->GetData("Multipoles"));
          if (!am) { am = new Multipoles(Vector3(pa->GetX(), pa->GetY(), pa->GetZ())); pa->SetData(am); }
          if (!bm) { bm = new Multipoles(Vector3(pb->GetX(), pb->GetY(), pb->GetZ())); pb->SetData(bm); }
          for (int l = 0; l < grid.size(); ++l)
          {
            if ( atype < btype )          
            {
              A(row_offset+l, col_pos) += 1.0 / ( ( am->center() / bohr ) - grid.lattice(l) ).length();
              A(row_offset+l, col_pos) -= 1.0 / ( ( bm->center() / bohr ) - grid.lattice(l) ).length();
            } else {
              A(row_offset+l, col_pos) -= 1.0 / ( ( am->center() / bohr ) - grid.lattice(l) ).length();
              A(row_offset+l, col_pos) += 1.0 / ( ( bm->center() / bohr ) - grid.lattice(l) ).length();
            } 
          }
        }
         
      }
      for (int l = 0; l < grid.size(); ++l) C(row_offset+l) = grid.field(l);
      row_offset += grid.size();
    }
  }
  void fillParameters()
  {
    for (std::map<std::string, int>::const_iterator i = param_map.begin(), ie = param_map.end(); i != ie; ++i)
      partable.set(i->first, "bcc", R(i->second));
  }

  BCCFitter(const BCCFitter& other);
  BCCFitter& operator = (const BCCFitter& other);

  ParameterTable& partable;
  std::vector<GeneralMultipoledMolecule>& mols;
  std::vector<Grid>& grids;
  std::vector<std::vector<std::string> >& types;

  std::map<std::string, int> param_map;
  std::vector<int> activeMols;

  std::vector<std::pair<std::string,std::string> > atom_type_array;
protected:
  virtual Eigen::VectorXd solve(Eigen::MatrixXd& A, Eigen::VectorXd& C)
  {
    Eigen::MatrixXd LSMat = A.transpose() * A;
    Eigen::MatrixXd BMat = A.transpose() * C;
    return LSMat.ldlt().solve(BMat);
  }
  Eigen::MatrixXd A;
  Eigen::VectorXd C;
  Eigen::VectorXd R;
};


#endif // BCC_FITTER_HEADER
