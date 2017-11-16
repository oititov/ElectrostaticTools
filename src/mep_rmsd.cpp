#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include <boost/program_options.hpp>

#include "iotools.h"
#include "elec_tools_helpers.h"
#include "elec_tools_errors.h"
#include "config.h"
#include "molecule.h"
#include "grid.h"

struct smarts
{
  smarts(std::string params)
  {
    std::istringstream iss(params);
    std::vector<std::string> tokens;
    std::copy(std::istream_iterator<std::string>(iss), 
              std::istream_iterator<std::string>(),
              std::back_inserter(tokens));

    pattern.Init(tokens[0]);
    double alpha = 0;
    if (tokens.size() > 1 ) 
      cosalpha = cos(atof(tokens[1].c_str())/180.0*3.14159265); 
    else 
      cosalpha = cos(        90.0           /180.0*3.14159265);
    if (tokens.size() > 2) rmax = atof(tokens[2].c_str())/bohr/bohr; else rmax = 9999999999999;
  }
  OpenBabel::OBSmartsPattern pattern;
  double cosalpha;
  double rmax;
};

int main(int argc, char* argv[])
{
  namespace po = boost::program_options;

  po::options_description all_opts;
  po::options_description generic_opts("MEP RMSD calculator v " + project_version + "\n\n"
      "Usage:\n"
      "mep_rmsd <options>\n");
  add_generic_options(generic_opts);
  po::options_description input("Input options");
  input.add_options()
    ("molecules,m", po::value<std::vector<std::string> >()->composing()->multitoken()->required(), "input files")
    ("grids,g", po::value<std::vector<std::string> >()->composing()->multitoken()->required(), "grid files")
    ("mask,M", po::value<std::vector<std::string> >()->composing(), "SMARTS mask of anchor atoms followed by "
                                                                    "minimum angle between neighbouring_atom-anchor_atom "
                                                                    "and anchor_atom-point vectors (90 deg default) "
                                                                    "followed by "
                                                                    "maximum distance from anchor atoms to ESP points, "
                                                                    " unlimited by default")
    ("no-coordinates-check", po::bool_switch()->default_value(false), "check the coordinates in molecules to match "
                                                                  "with the coordinates in the ESP files");
  all_opts.add(generic_opts);
  all_opts.add(input);

  po::variables_map vm;
  try 
  {
    po::store(po::parse_command_line(argc, argv, all_opts), vm);
    if (process_generic_opts(vm, all_opts, "MEP RMSD calculator")) return 0;
    po::notify(vm);
  } catch (const boost::program_options::required_option& e)
  {
    std::cout << all_opts << std::endl;
    std::cerr << "Error: " << e.what() << std::endl;
    return ELEC_TOOLS_BAD_OPTION;
  } catch (const boost::program_options::unknown_option& e)
  {
    std::cout << all_opts << std::endl;
    std::cerr << "Error: " << e.what() << std::endl;
    return ELEC_TOOLS_BAD_OPTION;
  }
 
  if (vm.count("help"))
  {
    std::cout << all_opts << std::endl;
    return 0;
  }

  std::vector<std::string> masks;
  std::vector<smarts> smartss;
  if ( vm.count("mask") ) 
  { 
    masks = vm["mask"].as<std::vector<std::string> >();
    for (int i = 0; i < masks.size(); ++i) smartss.push_back(smarts(masks[i]));
  }

  IOTools iotools;

  std::vector<std::string> gridfiles = vm["grids"].as<std::vector<std::string> >();
  std::vector<std::string> molfiles = vm["molecules"].as<std::vector<std::string> >();
  if ( gridfiles.size() != molfiles.size() )
  {
    std::cerr << "The number of grids does not match the number of molecules. Aborting." << std::endl; return -1;
  }

  std::vector<Vector3> coords;
  for (int i = 0; i < molfiles.size(); ++i)  
  {
    GeneralMultipoledMolecule mol; 
    mol.SetAtomTypesPerceived();
    Grid grid;
    int err = IOTools::read_mol(molfiles[i], "", mol);
    if (err) 
    { 
      std::cerr << "Error reading molecule: " << molfiles[i] << std::endl; 
      if (!vm["skip-unreadable"].as<bool>()) return err; else continue; 
    }
  
    err = IOTools::read_grid(grid, gridfiles[i], &coords);
    if (err)
    {
      std::cerr << "Error reading grid: " << gridfiles[i] << std::endl;
      if (!vm["skip-unreadable"].as<bool>()) return err; else continue;
    }
 
    if (! vm["no-coordinates-check"].as<bool>())
    {
      if (coords.size() != mol.NumAtoms())
      {
        std::cerr << "Error in molecule: " << molfiles[i] << " grid " << gridfiles[i] 
                  << " - The number of atoms in molecule does not match the number of atoms in grid file. Skipping." << std::endl;
        continue;
      }
      bool fail = false;
      for (int i = 0; i < coords.size(); ++i) 
      {
        OpenBabel::OBAtom* a = mol.GetAtom(i+1);
        Vector3 v(a->GetX(), a->GetY(), a->GetZ());
        if ((coords[i] - v/bohr).length_squared() > 1e-3)
        {
          std::cerr << "Error in molecule: " << molfiles[i] << " grid " << gridfiles[i] 
                  << " - The coordinates in molecule does not match the coordinates in grid file. Skipping." << std::endl;
          fail = true; break;
        }
      }
      if (fail) continue;
    }

// filter points
   
    std::vector<bool> active_grid(grid.size(), true);
    if ( smartss.size() ) 
    {
      active_grid = std::vector<bool>(grid.size(), false);
      std::vector<Vector3> centers;
      std::vector<double> calphas;
      std::vector<double> rmaxs;
      std::vector< std::vector<Vector3> > axes;
      for (int p = 0; p < smartss.size(); ++p) 
      {
        std::vector<std::vector<int> > matches;
        smartss[p].pattern.Match(mol, matches, OpenBabel::OBSmartsPattern::AllUnique);
        for (int j = 0; j < matches.size(); ++j)
        {
          OpenBabel::OBAtom* aa = mol.GetAtom(matches[j][0]);
          Vector3 va(aa->GetX(), aa->GetY(), aa->GetZ()); va*=bohr;
          centers.push_back(va);
          axes.push_back(std::vector<Vector3>());
          FOR_NBORS_OF_ATOM(ab, &*aa)
          {
            Vector3 vb(ab->GetX(), ab->GetY(), ab->GetZ()); vb*=bohr;
            axes.back().push_back( (vb-va).normalized() );
          }
          rmaxs.push_back(smartss[p].rmax);
          calphas.push_back(smartss[p].cosalpha);
        }
      }
      for (int g = 0; g < grid.size(); ++g)
      {
        Vector3 point = grid.lattice(g);
        for (int j = 0; j < centers.size(); ++j)
        {
          bool good = false;
          Vector3 r = point - centers[j];
          if ( r.length_squared() < rmaxs[j] )
          {
            r.normalize();
            std::vector<Vector3>& nbrs = axes[j];
            good = true;
            for (int l = 0; l < nbrs.size(); ++l)
              if (r*nbrs[l] > calphas[j]) // comparing cos(alpha)
              {
                good = false; break;
              }
          }
          if ( good ) 
          {
            active_grid[g] = true;
            break;
          }
        }
      }
    }

// get rmsd
    double rmsd = 0;
    int active_size=0;
    for (int k = 0; k < grid.size(); ++k)
    {
      if (!active_grid[k]) continue;
      ++active_size;
      double err = grid.field(k) - mol.MEP(grid.lattice(k)*bohr) ;
      rmsd += err*err;
    }
    std::cout << molfiles[i] << " : " << sqrt(rmsd/active_size) * 627.54 << " kcal/mol*e ; N = " <<active_size<< std::endl;
  }

  return 0;
}

