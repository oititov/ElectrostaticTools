/*************************************************************************
 * elec_tools_helpers.h - Electrostatic Tools helper functions.          *
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

#ifndef ELEC_TOOLS_HELPERS_HEADER
#define ELEC_TOOLS_HELPERS_HEADER

#include <string>
#include <exception>

#include <boost/program_options.hpp>

#include <openbabel/obconversion.h>

namespace po = boost::program_options;

#include "comment_skipper.h"
#include "multipole_placement_rules.h"
#include "config.h"
#include "elec_tools_errors.h"

inline void add_generic_options(po::options_description& generic_opts)
{
  generic_opts.add_options()
                  ("version", "print version string")
                  ("help", "show this message");
}

inline void create_mep_fitter_options(po::options_description& all_opts, po::positional_options_description& positional, std::string out_suffix)
{
  po::options_description input_opts("Input control");
  input_opts.add_options()
                ("input,I", po::value<std::string>()->required(), "<required> input molecule file")
                ("filetype,t", po::value<std::string>()->default_value("GUESS"), "input molecule file type (extension or GUESS)")
                ("grid,G", po::value<std::string>()->required(), "<required> input reference ESP grid file")
                ("placement-rules,p", po::value<std::string>(), "multipole placement rules")
                ("orient-rules,o", po::value<std::string>(), "multipole orient rules");

  po::options_description output_opts("Output control");
  output_opts.add_options()
                 ("output,O", po::value<std::string>(), ("output molecule file. Default: <input>." + out_suffix).c_str())
                 ("log,L", po::value<std::string>(), "log output file. Default: <output>.log")
                 ("force-output,f", po::bool_switch()->default_value(false), "owerwrite existing files");

  po::options_description fit_opts("MEP fit control");
  fit_opts.add_options()
              ("algorithm,a", po::value<std::string>()->default_value("SVD"), "optimization algorithm. "
                  "Valid values are: SVD, LDLT, PartialPivLU, FullPivLU, HouseholderQR, ColPivHouseholderQR, FullPivHouseholderQR")
              ("cutoff,c", po::value<double>()->default_value(1e7, "10^7"), "condition value cutoff (only for SVD optimizer)")
              ("recalculate-topology,r", po::bool_switch()->default_value(false), "ignore specified topology information")
              ("break-equivalency,b", po::bool_switch()->default_value(false), "remove topological equivalency constraint")
              ("refit", po::bool_switch()->default_value(false), "save existing multipoles instead of removing them");

  all_opts.add(input_opts);
  all_opts.add(output_opts);
  all_opts.add(fit_opts);

  positional.add("input", 1);
  positional.add("grid", 1);
  positional.add("output", 1);
}

inline bool process_generic_opts(po::variables_map& vm, po::options_description& all_opts, std::string progname)
{
  bool processed = false;
  if (vm.count("help"))
  {
    std::cout << all_opts << std::endl;
    processed = true;
  }
  if (vm.count("version"))
  {
    std::cout << progname <<" v " << project_version << std::endl;
    processed = true;
  }
  return processed;
}


#endif // ELEC_TOOLS_HELPERS_HEADER
