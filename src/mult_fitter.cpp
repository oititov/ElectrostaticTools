/*************************************************************************
 * mult_fitter.cpp - Program for optimization of atomic multipoles.      *
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
#include <string>
#include <exception>

#include <boost/program_options.hpp>

#include <openbabel/obconversion.h>

#include "config.h"
#include "mep_fitter.h"
#include "elec_tools_errors.h"
#include "elec_tools_helpers.h"
#include "iotools.h"

namespace po = boost::program_options;

int main(int argc, char const *const *argv)
{
  po::options_description all_opts;
  po::positional_options_description positional;
         
  po::options_description generic_opts("Multipole fitter v " + project_version + "\n\nUsage:\nmult_fitter <options>\nmult_fitter <in_mol> <in_grid> <out_mol> [options]\n\nOptions available");
  add_generic_options(generic_opts);
  all_opts.add(generic_opts);
  create_mep_fitter_options(all_opts, positional, "mmol");

  po::variables_map vm;
  try
  {
    po::store(po::command_line_parser(argc, argv).options(all_opts).positional(positional).run(), vm);
    if (process_generic_opts(vm, all_opts, "Multipole fitter")) return 0;
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

  // check input
  
  try
  {
    GeneralMultipoledMolecule multMol;
    int err = IOTools::read_mol(vm["input"].as<std::string>(), vm["filetype"].as<std::string>(), multMol); if (err) return err;
    
    std::ostream* logstrm;
    std::ostream* outstrm;

    std::string logname, outname;
    if ( ! vm.count("output") )
      outname = IOTools::add_filename_suffix(vm["input"].as<std::string>(), "", "mmol");
    else outname = vm["output"].as<std::string>();

    if ( ! vm.count("log") )
      logname = IOTools::add_filename_suffix(outname, "", "log");
    else logname = vm["log"].as<std::string>();

    IOTools iotools;
    
    iotools.output(outname, vm["force-output"].as< bool >());
    iotools.output(logname, vm["force-output"].as< bool >());
    
    if (iotools.error())
    {
      std::cerr << iotools.what() << " Aborting." << std::endl;
      return iotools.error();
    }
    
    Grid grid;
    err = IOTools::read_grid(grid, vm["grid"].as< std::string >()); if (err) return err;
    
    MultipoleOrientRules o_rules;
    err = IOTools::read_rules(o_rules, vm, &iotools.output(logname)); if (err) return err;
    
    MultipolePlacementRules p_rules;
    err = IOTools::read_rules(p_rules, vm, &iotools.output(logname)); if (err) return err;
    
// done with checks!
// try to fit something finally
    multMol.setOrientRules(o_rules);

    multMol.raiseDummies();

    iotools.output(logname) << "Multipole fitter v " << project_version << std::endl;
    
    MEPFitter *fitter;
    fitter = FitterFactory(iotools.output(logname)).getFitter(vm["algorithm"].as<std::string>());
    MEPSVDFitter *f;
    if (f = dynamic_cast<MEPSVDFitter *>(fitter)) f->setCutoff(vm["cutoff"].as<double>());
    
    fitter->setRestrictEquivalent( ! vm["break-equivalency"].as<bool>() );
    fitter->setRecalcTopology( vm["recalculate-topology"].as<bool>() );
    fitter->setRefit( vm["refit"].as<bool>() );
    fitter->setMolecule(&multMol);
    fitter->setMultipolePlacementRules(p_rules);
    fitter->setGrid(grid);
    fitter->fit();
    multMol.lowerDummies();
    fitter->printResult(iotools.output(outname));
    delete fitter;
  } catch (const std::exception& e)
  {
    std::cerr << e.what() << std::endl;
    return ELEC_TOOLS_FITTING_FAILURE;
  }

  return 0;
}
