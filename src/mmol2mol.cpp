/*************************************************************************
 * mmol2mol.cpp - Program for conversion of mmol files to common         *
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

#include <openbabel/obconversion.h>

//#include "mep_fitter.h"
#include "molecule.h"
#include "comment_skipper.h"
#include "elec_tools_errors.h"
#include "elec_tools_helpers.h"
#include "mcc_creator.h"
#include "multipole_setter.h"
#include "config.h"
#include "iotools.h"
#include "par_table.h"
#include "bcc_typer.h"

namespace po = boost::program_options;

int main(int argc, char const *const *argv)
{
  po::options_description all_opts;

  po::options_description generic_opts("MMol to common format converter v " + project_version + "\n\n"
      "Usage:\n"
      "mmol2mol <options>\n"
      "mmol2mol <in_mol> <out_mol> [options]\n\n"
      "Options available");
  add_generic_options(generic_opts);

  po::options_description input_opts("Input control");
  input_opts.add_options()
                ("input,I", po::value<std::string>()->required(), "<required> input molecule file")
                ("input-type", po::value<std::string>()->default_value("GUESS"), "input molecule file type")
                ("orient-rules,o", po::value<std::string>(), "multipole orient rules");


  po::options_description output_opts("Output control");
  output_opts.add_options()
                 ("output,O", po::value<std::string>(), "output molecule file. Default: <input>.mol2")
                 ("output-type,t", po::value<std::string>()->default_value("GUESS"), "output molecule file type extension")
                 ("force-output,f", po::bool_switch()->default_value(false), "owerwrite existing files");

  po::options_description conv_opts("Conversion control");
  conv_opts.add_options()
               ("radius,r", po::value<double>()->default_value(0.1, "0.1"), "charge cluster radius in Angstroms")
               ("ignore-dipole,d", po::bool_switch()->default_value(false), "do not include dipole in charge cluster")
               ("ignore-quadrupole,q", po::bool_switch()->default_value(false), "do not include quadrupole in charge cluster")
               ("mcc-mask,M", po::value<std::vector<std::string> >(), "SMARTS mask of atoms which multipoles "
                                                                      "are to be converted to MCC (multiple keys are posible)");
               
  po::options_description mod_opts("Multipole modification control");
  mod_opts.add_options()
               ("copy-multipoles", po::value<std::vector<std::string> >(), "molecule files from which the multipoles will be"
                                                                        "copied. No atom matching is performed.") 
               ("set-charge,C", po::value<std::vector<std::string> >(), "a pair of SMARTS mask and charge value, eg \"[Cl-] -1\"."
                                                                        " Beware of breaking molecule's total charge!")
               ("set-dipole,D", po::value<std::vector<std::string> >(), "a pair of SMARTS mask and dipole value, "
                                                                        "eg \"[Br] 0.0 0.0 1.0\"")
               ("set-quadrupole,Q", po::value<std::vector<std::string> >(), "a pair of SMARTS mask and quadrupole value, "
                                                                        "eg \"[I]a 1.5\"")
               ("add-charge", po::value<std::vector<std::string> >(), "a pair of SMARTS mask and charge value to add on atom,"
                                                                        "eg \"[Cl-] -1\"."
                                                                        " Beware of breaking molecule's total charge!")
               ("add-dipole", po::value<std::vector<std::string> >(), "a pair of SMARTS mask and dipole value to add on atom, "
                                                                        "eg \"[Br] 0.0 0.0 1.0\"")
               ("add-quadrupole", po::value<std::vector<std::string> >(), "a pair of SMARTS mask and quadrupole value "
                                                                        "to add on atom, "
                                                                        "eg \"[I]a 1.5\"");

  po::options_description bcc_opts("Bond Charge Correction control");
  bcc_opts.add_options()
               ("bcc", po::value<std::string> (), "BCC rule file.")
               ("bcc-types", po::value<std::string>(), "A file with bcc atom types. Required for BCC.");
              
  all_opts.add(generic_opts);
  all_opts.add(input_opts);
  all_opts.add(output_opts);
  all_opts.add(mod_opts);
  all_opts.add(bcc_opts);
  all_opts.add(conv_opts);

  po::positional_options_description positional;
  positional.add("input", 1);
  positional.add("output", 1);

  po::variables_map vm;
  try
  {
    po::store(po::command_line_parser(argc, argv).options(all_opts).positional(positional).run(), vm);
    if (process_generic_opts(vm, all_opts, "MMol to common format converter")) return 0;
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
    std::string ofiletype = vm["output-type"].as<std::string>();
    std::transform(ofiletype.begin(), ofiletype.end(), ofiletype.begin(), tolower);
    if (ofiletype == "guess")
    {
      if ( ! vm.count("output") )
      {
        ofiletype = "mol2";
      } else {
        std::string t = vm["output"].as<std::string>();
        ofiletype = t.substr(t.find_last_of(".")+1);
      }
    }

    std::string outname;
    if ( ! vm.count("output") )
      outname = IOTools::add_filename_suffix(vm["input"].as<std::string>(), "mcc", ofiletype);
    else outname = vm["output"].as<std::string>();
    
    IOTools iotools;
    iotools.output(outname, vm["force-output"].as< bool >());
    if (iotools.error()) 
    {
      std::cerr << iotools.what() << " Aborting." << std::endl;
      return iotools.error();
    }

    double radius = vm["radius"].as< double >();
    if (radius > 1.0)
      std::cout << "Your MCC radius is larger than 1.0 Angstroms. This is strange since it is expected to be as small as possible." << std::endl;

    GeneralMultipoledMolecule mol;
    std::string type_read;
    int err = IOTools::read_mol( vm["input"].as< std::string >(), vm["input-type"].as< std::string >(), mol, &type_read); if(err) return err;
    
    if (type_read != "mmol")
    {
      MultipoleOrientRules o_rules;
      err = IOTools::read_rules(o_rules, vm, &std::clog); if(err) return err;
      mol.setOrientRules(o_rules);
    }
    
    MultipoleSetter multSetter;
    multSetter.setMolecule(&mol);

    if (vm.count("copy-multipoles"))
    {
      std::vector<Multipoles*> mults;
      std::vector<std::string> molfiles = vm["copy-multipoles"].as<std::vector<std::string> > ();
      MultipoleOrientRules refrules = mol.getOrientRules();
      for (int i = 0; i < molfiles.size(); ++i)
      {
        GeneralMultipoledMolecule refmol;
        std::string type_read;
        int err = IOTools::read_mol( molfiles[i], "GUESS", refmol, &type_read); 
          if(err) return err;
        if (type_read == "mmol")
          refrules += refmol.getOrientRules();
       
        FOR_ATOMS_OF_MOL(atom, refmol)
          mults.push_back(new Multipoles(static_cast<Multipoles*>(atom->GetData("Multipoles"))));
      }
      for (int i = 1; i <= mol.NumAtoms(); ++i)
      {
        if (i > mults.size()+1) 
        {
          std::cout << "Warning! Copying multipoles: not enougth multipoles!" << std::endl;
          break;
        } 
        OpenBabel::OBAtom* a = mol.GetAtom(i);
        Multipoles* m = new Multipoles(mults[i-1]);
        m->moveTo(Vector3(a->GetX(), a->GetY(), a-> GetZ()));
        if ( a->HasData("Multipoles") ) a->DeleteData("Multipoles");
        a->SetData(m);
      }
      mol.setOrientRules(refrules);
      for (int i = 0; i < mults.size(); ++i) delete mults[i];
    }    

    if (vm.count("set-charge"))
    {
      std::vector<std::string> vs = vm["set-charge"].as<std::vector<std::string> > ();
      for (std::vector<std::string>::iterator i = vs.begin(), ie = vs.end(); i != ie; ++i)
      {
        multSetter.placeMonopole(*i);
      }
    }
    
    if (vm.count("set-dipole"))
    {
      std::vector<std::string> vs = vm["set-dipole"].as<std::vector<std::string> > ();
      for (std::vector<std::string>::iterator i = vs.begin(), ie = vs.end(); i != ie; ++i)
      {
        multSetter.placeDipole(*i);
      }
    }
   
    if (vm.count("set-quadrupole"))
    {
      std::vector<std::string> vs = vm["set-quadrupole"].as<std::vector<std::string> > ();
      for (std::vector<std::string>::iterator i = vs.begin(), ie = vs.end(); i != ie; ++i)
      {
        multSetter.placeQuadrupole(*i);
      }
    }
   
    if (vm.count("add-charge"))
    {
      std::vector<std::string> vs = vm["add-charge"].as<std::vector<std::string> > ();
      for (std::vector<std::string>::iterator i = vs.begin(), ie = vs.end(); i != ie; ++i)
      {
        multSetter.placeMonopole(*i, true);
      }
    }

    if (vm.count("add-dipole"))
    {
      std::vector<std::string> vs = vm["add-dipole"].as<std::vector<std::string> > ();
      for (std::vector<std::string>::iterator i = vs.begin(), ie = vs.end(); i != ie; ++i)
      {
        multSetter.placeQuadrupole(*i, true);
      }
    }

    if (vm.count("add-quadrupole"))
    {
      std::vector<std::string> vs = vm["add-quadrupole"].as<std::vector<std::string> > ();
      for (std::vector<std::string>::iterator i = vs.begin(), ie = vs.end(); i != ie; ++i)
      {
        multSetter.placeQuadrupole(*i, true);
      }
    }
  
    if (vm.count("bcc"))
    {
      if (!vm.count("bcc-types")) 
      {
        std::cerr << "Error: BCC requested but no BCC types file specified. Aborting." << std::endl;
        return -1;
      } 
      
      std::ifstream partypes_file;
      int err = IOTools::open_in_file(partypes_file, vm["bcc"].as<std::string>());
      if (err) return err;
      if (!partypes_file) { std::cerr << "Error: Cannot read parameter types file" << std::endl; return 1; }
      CommentSkipper skipper;
      std::istringstream partypes_strm(skipper.skip(partypes_file));
      ParameterTable partable(partypes_strm);
      partypes_file.close();

      std::vector<std::string> typesv;
      std::ifstream typestrm;
      IOTools::open_in_file(typestrm, vm["bcc-types"].as<std::string>());
      while (typestrm)
      {
        std::string type; typestrm >> type;
        if (!type.empty()) typesv.push_back(type);
      }
         
      BCCAtomTyper atyper(typesv);
      BCCBondTyper btyper(atyper);
      FOR_BONDS_OF_MOL(bond, mol)
      {
        std::string bond_type = btyper.type(&(*bond));
        if (partable.exists(bond_type, "bcc"))
        {
          double bcc = partable.get(bond_type, "bcc");
          OpenBabel::OBAtom* pa = bond->GetBeginAtom();
          OpenBabel::OBAtom* pb = bond->GetEndAtom();
          std::string atype = atyper.type(pa);
          std::string btype = atyper.type(pb);
          Multipoles* am = static_cast<Multipoles*>(pa->GetData("Multipoles"));
          Multipoles* bm = static_cast<Multipoles*>(pb->GetData("Multipoles"));
          if (!am) { am = new Multipoles(Vector3(pa->GetX(), pa->GetY(), pa->GetZ())); pa->SetData(am); }
          if (!bm) { bm = new Multipoles(Vector3(pb->GetX(), pb->GetY(), pb->GetZ())); pb->SetData(bm); }
          if ( atype < btype )
          {
            am->setMonopole(am->monopole().value() + bcc);
            bm->setMonopole(bm->monopole().value() - bcc);
          } else {
            am->setMonopole(am->monopole().value() - bcc);
            bm->setMonopole(bm->monopole().value() + bcc);
          }
        }
      }
      mol.syncronizeCharges(); 
    }
 
    if (vm.count("mcc-mask"))
    {  
      std::vector< std::string > mcc_patterns = vm["mcc-mask"].as< std::vector< std::string > >();
      MCCCreator mccCreator;
      mccCreator.setMolecule(&mol);
      mccCreator.ignoreDipole(vm["ignore-dipole"].as< bool >());
      mccCreator.ignoreQuadrupole(vm["ignore-quadrupole"].as< bool >());
      mccCreator.setRadius(radius);
      mccCreator.setEps(1e-4);
      mccCreator.installMCCs(mcc_patterns);
    }
    
    iotools.write_mol( outname, ofiletype, mol );
    if (iotools.error())
    {
      std::cerr << iotools.what() << " Aborting." << std::endl;
      return iotools.error();
    }

  } catch (const std::exception& e)
  {
    // generally nothing should fail here but just in case
    std::cerr << "Error: " << e.what() << std::endl;
    return ELEC_TOOLS_UNKNOWN_ERROR;
  }


  return 0;
}
