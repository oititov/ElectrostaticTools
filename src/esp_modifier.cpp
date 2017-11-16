/*************************************************************************
 * esp_modifier.cpp - Program for conversion of mmol files to common     *
 *                file formats.                                          *
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

#include <fstream>
#include <iomanip>
#include <cstring>
#include <vector>
#include <cctype>

#include <boost/program_options.hpp>

#include "molecule.h"
#include "grid.h"
#include "comment_skipper.h"
#include "elec_tools_errors.h"
#include "elec_tools_helpers.h"
#include "config.h"
#include "iotools.h"


namespace po = boost::program_options;

int main(int argc, char const *const *argv)
{
  po::options_description all_opts;

  po::options_description generic_opts("ESP modifier " + project_version + "\n\n"
      "Usage:\n"
      "esp_modifier <options>\n"
      "esp_modiifer <in_esp> <out_esp> [options]\n\n"
      "Options available");
  add_generic_options(generic_opts);

  po::options_description input_opts("Input control");
  input_opts.add_options()
                ("input,I", po::value<std::string>()->required(), "<required> input ESP file")
                ("add,A", po::value<std::vector<std::string> >(), "molecules, which MEP should be added to the ESP")
                ("remove,R", po::value<std::vector<std::string> >(), "molecules, which MEP should be removed from the ESP")
		("coordinates,c", po::value<std::string>(), "replace atom coordinates in ESP with the specified ones");
  
  po::options_description output_opts("Output control");
  output_opts.add_options()
                 ("output,O", po::value<std::string>()->required(), "output ESP file")
                 ("force-output,f", po::bool_switch()->default_value(false), "owerwrite existing files");

  all_opts.add(generic_opts);
  all_opts.add(input_opts);
  all_opts.add(output_opts);

  po::positional_options_description positional;
  positional.add("input", 1);
  positional.add("output", 1);

  po::variables_map vm;
  try
  {
    po::store(po::command_line_parser(argc, argv).options(all_opts).positional(positional).run(), vm);
    if (process_generic_opts(vm, all_opts, "ESP modifier")) return 0;
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

  try
  {
    std::string outname = vm["output"].as<std::string>();
    std::string inname = vm["input"].as<std::string>();
    std::string cname;
    int err = 0;
    GeneralMultipoledMolecule cmol;
    
    if (vm.count("coordinates")) 
    {
      err = IOTools::read_mol(vm["coordinates"].as<std::string>(), "guess", cmol); if(err) return err;
    }

    IOTools iotools;
    
    iotools.output(outname, !vm["force-output"].as< bool >());
    
    if (iotools.error())
    {
      std::cerr << iotools.what() << " Aborting." << std::endl;
      return iotools.error();
    }
    
    std::vector< Vector3 > coords;
    Grid grid;
    err = IOTools::read_grid(grid, inname, &coords); if (err) return err;
    std::vector<std::string> add, remove;
    
    if (vm.count("add")) add = vm["add"].as< std::vector < std::string > >();
    if (vm.count("remove")) remove = vm["remove"].as<std::vector<std::string> >();

    GeneralMultipoledMolecule mol;
    for (std::vector<std::string>::const_iterator i = add.begin(); i != add.end(); i++)
    {
      err = IOTools::read_mol(*i, "guess", mol); if(err) return err;
      for (int k = 0; k < grid.size(); ++k) grid.modifyField(k, mol.MEP(grid.lattice(k)*bohr));
      grid.charge(mol.GetTotalCharge());
    }
    for (std::vector<std::string>::const_iterator i = remove.begin(); i != remove.end(); i++)
    {
      err = IOTools::read_mol(*i, "guess", mol); if(err) return err;
      for (int k = 0; k < grid.size(); ++k) grid.modifyField(k, -mol.MEP(grid.lattice(k)*bohr));
      grid.charge(-mol.GetTotalCharge());
    }

    if (vm.count("coordinates"))
    {
      coords.clear();
      coords.reserve(cmol.NumAtoms());
      for (int i = 1; i <= cmol.NumAtoms(); ++i)
      {
        OpenBabel::OBAtom* a = cmol.GetAtom(i);
        coords.push_back(Vector3(a->GetX(), a->GetY(), a->GetZ())/bohr);
      }
    }

    grid.toESP(iotools.output(outname), &coords);

  } catch (const std::exception& e)
  {
    // generally nothing should fail here but just in case
    std::cerr << "Error: " << e.what() << std::endl;
    return ELEC_TOOLS_UNKNOWN_ERROR;
  }


  return 0;
}
