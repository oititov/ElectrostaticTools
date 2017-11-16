/*************************************************************************
 * ep_fitter.cpp - Program for optimization of extra-point positions.    *
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
#include <string>
#include <cstring>
#include <cctype>

#include <boost/bind.hpp>
#include <boost/program_options.hpp>

#include <openbabel/obconversion.h>

#include "mep_fitter.h"
#include "optimizer.h"
#include "smarts_matcher.h"
#include "comment_skipper.h"
#include "elec_tools_errors.h"
#include "config.h"
#include "elec_tools_helpers.h"
#include "iotools.h"

void setPos(OpenBabel::OBAtom *ep, OpenBabel::vector3 hal_pos, OpenBabel::vector3 axis, double pos, std::ostream* log)
{
  *log << "  Setting EP position to: " << pos << std::endl;
  ep->SetVector(hal_pos + axis * pos);
}

double test(MEPFitter *fitter, std::ostream* log)
{
  fitter->createMatrix();
  fitter->makeFit();
  fitter->getResults();
  *log << "  RMSD: " << (fitter->rmsd()*629.5) << " kcal/mol" << std::endl;
  return fitter->rmsd();
}

namespace po = boost::program_options;

int main(int argc, char const *const *argv)
{
  po::options_description all_opts;
  po::positional_options_description positional;

  po::options_description generic_opts("EP position fitter v " + project_version + "\n\n"
      "Usage:\n"
      "ep_fitter <options>\n"
      "ep_fitter <in_mol> <grid> <out_mol> [options]\n\n"
      "Options available");
  add_generic_options(generic_opts);
  all_opts.add(generic_opts);
  create_mep_fitter_options(all_opts, positional, "ep.mmol");

  po::options_description ep_fit_opts("EP position fit control");
  ep_fit_opts.add_options()
                 ("host-mask,M", po::value< std::vector< std::string > >(),
                  "SMARTS mask of EP hosts (multiple keys are posible). Note: EP scan is performed along local Z-axis. Note: experimental feature for non-halogens.")
                 ("fixed-position,x", po::value< double >(), "fixed EP distance from halogen atom in Angstroms. No optimization will be applied.")
                 ("init-position,d", po::value< double >()->default_value(1.5, "1.5"), "initial EP distance from halogen atom in Angstroms")
                 ("init-step,s", po::value< double >()->default_value(0.1, "0.1"), "initial position step for search in Angstroms")
                 ("precision,e", po::value< double >()->default_value(1e-3, "0.001"), "required EP position precision (epsilon)")
                 ("max-steps,m", po::value< int >()->default_value(100, "100"), "maximum number of position optimization algorithm steps")
                 ("keep-multipoles,k", po::bool_switch()->default_value(false), "keep multipoles on EP host atoms");

  all_opts.add(ep_fit_opts);

  po::variables_map vm;
  try
  {
    po::store(po::command_line_parser(argc, argv).options(all_opts).positional(positional).run(), vm);
    if (process_generic_opts(vm, all_opts, "EP position fitter")) return 0;
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
    GeneralMultipoledMolecule multMol;
    int err = IOTools::read_mol(vm["input"].as<std::string>(), vm["filetype"].as<std::string>(), multMol);
    if (err) return err;
    multMol.raiseDummies();

    std::string logname, outname;

    if (!vm.count("output"))
      outname = IOTools::add_filename_suffix(vm["input"].as< std::string >(), "ep", "mmol");
    else outname = vm["output"].as< std::string >();

    if (!vm.count("log"))
      logname = IOTools::add_filename_suffix(outname, "", "log");
    else logname = vm["log"].as< std::string >();

    IOTools iotools;
    
    iotools.output(outname, !vm["force-output"].as< bool >());
    iotools.output(logname, !vm["force-output"].as< bool >());
  
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

    multMol.setOrientRules(o_rules);

    // ready to process molecule

    iotools.output(logname) << "EP position fitter v " << project_version << std::endl;

    multMol.createOrientMap();
    multMol.calcOrientMatrices();

    std::vector< double > charges;
    for (int i = 1; i <= multMol.NumAtoms(); ++i) charges.push_back(multMol.GetAtom(i)->GetPartialCharge());

    std::vector< OpenBabel::vector3 > axis, positions;
    std::vector< OpenBabel::OBAtom * > eps;

    // first find all EP hosts

    bool keep_multipoles = vm["keep-multipoles"].as< bool >();
    std::vector< std::string > host_patterns;
    if (vm.count("host-mask"))
      host_patterns = vm["host-mask"].as< std::vector< std::string > >();
    SMARTSMatcher matcher;
    matcher.addPatterns(host_patterns);
    std::vector< int > host_indices = matcher.getMatches(multMol);

    if (host_indices.empty()) // no need for extra-points
    {
      iotools.output(outname) << multMol;
      return 0;
    }

    // add EPs

    double init_position = vm["init-position"].as<double>();
    if (vm.count("fixed-position")) init_position = vm["fixed-position"].as<double>(); 
    
    multMol.BeginModify();
    for (int i = 0; i < host_indices.size(); ++i)
    {
      OpenBabel::OBAtom *a = multMol.GetAtom(host_indices[i]);
      Matrix3 orient = multMol.matrix(host_indices[i]);
      // check if vector is aligned outside of molecule
      OpenBabel::vector3 v = OpenBabel::vector3(orient(0, 2), orient(1, 2), orient(2, 2));
      v.normalize();
      OpenBabel::vector3 v0 = a->GetVector();
      OpenBabel::vector3 vi;
      using OpenBabel::OBAtomBondIter;
      FOR_NBORS_OF_ATOM(nbr, *a)
      {
        vi = nbr->GetVector();
        break;
      }
      if (( v0 + v - vi ).length_2() < ( v0 - v - vi ).length_2()) v = -v;
      positions.push_back(a->GetVector());
      axis.push_back(v);
      eps.push_back(multMol.NewAtom());
      eps.back()->SetAtomicNum(99);
      eps.back()->SetVector(positions.back() + axis.back() * init_position);
      multMol.AddBond(a->GetIdx(), eps.back()->GetIdx(), 1);
    }
    multMol.EndModify();

    // restore partial atomic charges
    for (int i = 1; i <= charges.size(); ++i)
    {
      multMol.GetAtom(i)->SetPartialCharge(charges[i - 1]);
    }
    multMol.SetPartialChargesPerceived();

    // initialize MEP fitter

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
    fitter->addDummyCenters();
    if (vm["keep-multipoles"].as<bool>())
      fitter->placeAtomParameters();
    else
      fitter->placeAtomParameters(host_patterns);
    
    if (vm.count("fixed-position"))
    {
      // no optimization needed
      iotools.output(logname) << "Using fixed EP positions. No optimization required." << std::endl;
      test(fitter, &iotools.output(logname));
      fitter->printResult(iotools.output(outname));
      iotools.output(logname) << "EP Positions : ";
      for (int i = 0; i < eps.size(); ++i) iotools.output(logname) << init_position << " ";
      iotools.output(logname) << std::endl;
    } else {
      // set up the optimizer
      Optimizers::Optimizer *optimizer = Optimizers::OptimizerFactory::getOptimizer("Nelder-Mead");
      optimizer->setEpsilon(vm["precision"].as<double>());
      optimizer->setMaxIter(vm["max-steps"].as<int>());
      optimizer->setTargetFunction(boost::bind(test, fitter, &iotools.output(logname)));
      for (int i = 0; i < eps.size(); ++i)
        optimizer->addParameter
                   (
                     boost::bind(setPos, eps[i], positions[i], axis[i], _1, &iotools.output(logname)),
                       vm["init-position"].as<double>(), vm["init-step"].as<double>()
                   );
      // finally optimize
      bool converged = optimizer->optimize();
      Optimizers::Optimizer::DataPoint point = optimizer->result();
      optimizer->setAll(point.x);
      test(fitter, &iotools.output(logname));
      fitter->printResult(iotools.output(outname));
      if (converged)
        iotools.output(logname) << "Optimization converged." << std::endl;
      else
        iotools.output(logname) << "Optimization NOT converged." << std::endl;
      iotools.output(logname) << "Final point:" << std::endl
               << "RMSD: " << point.f * 629.5 << "kcal/mol" << std::endl
               << "EP Positions : ";
      for (int i = 0; i < point.x.size(); ++i) iotools.output(logname) << point.x[i] << " ";
      iotools.output(logname) << std::endl;
    }
   
  } catch (const std::exception& e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    return ELEC_TOOLS_UNKNOWN_ERROR;
  }

  return 0;
}
