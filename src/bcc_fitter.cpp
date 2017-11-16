#include <iostream>
#include <vector>
#include <string>

#include <boost/program_options.hpp>

#include "iotools.h"
//#include "elec_tools_helpers.h"
#include "elec_tools_errors.h"
#include "config.h"
#include "molecule.h"
#include "grid.h"
#include "par_table.h"
#include "bcc_fitter.h"

int main(int argc, char* argv[])
{
  namespace po = boost::program_options;

  po::options_description all_opts;
  po::options_description generic_opts("BCC value fitter v " + project_version + "\n\n"
      "Usage:\n"
      "bcc_fitter <options>\n");
  generic_opts.add_options()("help,h", "print help message");
  all_opts.add(generic_opts);
  po::options_description input("Input options");
  input.add_options()
    ("input,i", po::value<std::string>()->required(), "input parameter file")
    ("molecules,M", po::value<std::vector<std::string> >()->composing()->multitoken()->required(), "input files")
    ("grids,G", po::value<std::vector<std::string> >()->composing()->multitoken()->required(), "grid files")
    ("types,T", po::value<std::vector<std::string> >()->composing()->multitoken()->required(), "atom types files")
    ("no-coordinates-check", po::bool_switch()->default_value(false), "check the coordinates in molecules to match "
                                                                  "with the coordinates in the ESP files")
    ("skip-unreadable", po::bool_switch()->default_value(false), "skip molecules and grids in unreadable files");
  all_opts.add(input);
  po::options_description output("Output options");
  output.add_options()  
    ("output,o", po::value<std::string>(), "output parameter file")
    ("force-output,f", po::bool_switch()->default_value(false), "owerwrite existing files");
  all_opts.add(output);

  po::variables_map vm;
  try 
  {
    po::store(po::parse_command_line(argc, argv, all_opts), vm);
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
 
 
  std::ifstream partypes_file;
  int err = IOTools::open_in_file(partypes_file, vm["input"].as<std::string>());
  if (err) return err;
  if (!partypes_file) { std::cerr << "Error: Cannot read parameter types file" << std::endl; return 1; }
  CommentSkipper skipper;
  std::istringstream partypes_strm(skipper.skip(partypes_file));
  ParameterTable partable(partypes_strm);
  partypes_file.close();

  IOTools iotools;
  
  std::string opartypes;
  if (vm.count("output")) opartypes = vm["output"].as<std::string>();
  else opartypes = IOTools::add_filename_suffix(vm["input"].as<std::string>(), "fit", "out");
  
  std::ostream& opartypesstrm = iotools.output(opartypes, vm["force-output"].as<bool>());
  if (iotools.error())
  {
    std::cerr << iotools.what() << " Aborting." << std::endl;
    return iotools.error();
  }
 
  std::vector<GeneralMultipoledMolecule> molecules;
  std::vector<Grid> grids;

  std::vector<std::string> gridfiles = vm["grids"].as<std::vector<std::string> >();
  std::vector<std::string> molfiles = vm["molecules"].as<std::vector<std::string> >();
  std::vector<std::string> typefiles = vm["types"].as<std::vector<std::string> >();
  if ( gridfiles.size() != molfiles.size() && molfiles.size() != typefiles.size() )
  {
    std::cerr << "The numbers of grids, molecules and types does not match. Aborting." << std::endl; return -1;
  }

  std::vector<std::string> mol_read;
  std::vector<Vector3> coords;
  std::vector<std::vector<std::string> > types;
  for (int i = 0; i < molfiles.size(); ++i)  
  {
    GeneralMultipoledMolecule mol; 
    Grid grid;
    std::vector<std::string> typesv;

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

    std::ifstream typestrm;
    IOTools::open_in_file(typestrm, typefiles[i]);
    while (typestrm) 
    {
      std::string type; typestrm >> type; 
      if (!type.empty()) typesv.push_back(type);
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

    molecules.push_back(mol);
    grids.push_back(grid);
    types.push_back(typesv);
    mol_read.push_back(molfiles[i]);
  }

  BCCFitter fitter(partable, molecules, grids, types);
  fitter.makeFit(); 
  
  opartypesstrm << partable << std::endl;

  std::vector<double> rmsds = fitter.rmsds();

  opartypesstrm << "/* Fit results: " << std::endl;
  for (int i = 0; i < mol_read.size(); ++i) opartypesstrm << mol_read[i] << " : " << 
                                 (rmsds[i*2]*627.54) << " -> " << (rmsds[i*2+1]*627.54) << " kcal/mol" << std::endl;
  opartypesstrm << "Overall : " << rmsds[rmsds.size()-2]*627.54 << 
                                  " -> " << rmsds[rmsds.size()-1]*627.54 << " kcal/mol" << std::endl;
  opartypesstrm << std::endl;
  partable.printUseCount(opartypesstrm);
  opartypesstrm << std::endl;
  
  opartypesstrm << "*/" << std::endl;
  
  return 0;
}

