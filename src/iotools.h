/*************************************************************************
 * iotools.h - Tools for fast input/output of molecule files             *
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

#ifndef IOTOOLS_H
#define IOTOOLS_H

#include <string>
#include <fstream>
#include <map>

#include <boost/program_options.hpp>

#include <openbabel/obconversion.h>

#include "molecule.h"
#include "grid.h"
#include "comment_skipper.h"
#include "multipole_placement_rules.h"
#include "elec_tools_errors.h"
#include "config.h"

class IOTools
{
public:
  IOTools() : err(0), human_error("") { }
  ~IOTools() 
  {
    for (map_type::iterator i = ostreams.begin(), ie = ostreams.end(); i != ie; ++i) 
    {
      i->second->close();
      delete i->second;
    }
  }
  
  int error() { return err; }
  const std::string& what() { return human_error; }
  void clear_errors() { err = 0; human_error = ""; }
  
  std::ostream& output(const std::string& filename, bool force = false)
  {
    map_type::iterator handle = ostreams.find(filename);
    if ( handle != ostreams.end() ) return *(handle->second);
    
    if ( !( force || ! check_file_exists(filename) ) )
    {
      human_error = "Error: output file exists." ;
      err = ELEC_TOOLS_OUTPUT_EXISTS;
      return std::cerr;
    }
    std::ofstream* out_strm = new std::ofstream(filename.c_str());
    if ( !out_strm->good() ) {
      human_error = std::string("Error opening output file: \"") + filename + "\"";
      err = ELEC_TOOLS_UNWRITABLE_OUTPUT;
      return std::cerr;
    }
    ostreams[filename] = out_strm;
    return *out_strm;
  }  
  
  void close(const std::string& filename)
  {
    map_type::iterator handle = ostreams.find(filename);
    if ( handle != ostreams.end() ) 
    {
      handle->second->close();
      delete handle->second;
      ostreams.erase(handle);
    }
  }
  
  static int read_grid(Grid& grid, const std::string& filename, std::vector<Vector3> *coords = 0)
  {
    std::ifstream igstream; 
    int err = open_in_file(igstream, filename); if (err) return err;
    grid.fromESP(igstream, coords);
    igstream.close();
    return 0;
  }
  
  static int read_mol( const std::string& filename, const std::string& type, GeneralMultipoledMolecule& mol, std::string* read_type = 0)
  {
    int lastindex = filename.find_last_of(".");
    std::string ext = filename.substr(lastindex+1, filename.length());
    std::transform(ext.begin(), ext.end(), ext.begin(), tolower);
    std::string type_l = type;
    std::transform(type.begin(), type.end(), type_l.begin(), tolower);
        
    std::ifstream molfilestrm;
    int err = open_in_file(molfilestrm, filename); if(err) return err;
    
    if ( ext == "mmol" || type_l == "mmol")
    {
      int err = read_mmol(molfilestrm, mol); if (err) return err; 
      if (read_type) *read_type = "mmol";
    } else {
      type_l = type;
      int err = read_obmol(molfilestrm, mol, filename, type_l); if (err) return err;  
      if (read_type) *read_type = type_l;
      mol.reverseSyncronizeCharges();
    }    
    molfilestrm.close();
    
    return 0;
  }
  
  void write_mol( const std::string& filename, const std::string& type, GeneralMultipoledMolecule& mol )
  {
    int lastindex = filename.find_last_of(".");
    std::string ext = filename.substr(lastindex+1, filename.length());
    std::transform(ext.begin(), ext.end(), ext.begin(), tolower);
    std::string type_l = type;
    std::transform(type.begin(), type.end(), type_l.begin(), tolower);
        
    if ( ext == "mmol" || type_l == "mmol")
    {
      write_mmol(output(filename), mol); 
    } else {
      write_obmol(output(filename), mol, filename, type); 
    }    
  }
  
  static bool check_file_exists(const std::string& filename)
  {
    std::ifstream om_test(filename.c_str());
    if (om_test.good())
    {
      om_test.close();
      return true;
    } else return false;
  }
  static std::string add_filename_suffix(std::string name, std::string suffix, std::string filetype = "")
  {
    std::transform(suffix.begin(), suffix.end(), suffix.begin(), tolower);
    std::transform(filetype.begin(), filetype.end(), filetype.begin(), tolower);
    int lastindex = name.find_last_of(".");
    std::string result = name.substr(0, lastindex);
    if (!suffix.empty()) result += "." + suffix;
    if (!filetype.empty()) result += "." + filetype;
    return result;
  }
  static std::string trim_path(const std::string& filename)
  {
    int lastindex = filename.find_last_of("/");
    std::string result = filename.substr(lastindex+1, filename.length());
    lastindex = result.find_last_of("\\");
    return result.substr(lastindex+1, result.length());
  } 
  // some template magic
  
  
  
  template <class RulesT>
  static int read_rules(RulesT& rules, boost::program_options::variables_map& vm, std::ostream* logstrm)
  {
    RuleHelper<RulesT> helper;
    std::string paramname = helper.param_name();
    std::string filename;
    if (vm.count(paramname))
      filename = vm[paramname].as< std::string >();
    else filename = config_path + helper.default_filename();
    std::istream* fstrm;
    std::ifstream *test_strm = new std::ifstream(filename.c_str());
    if (!test_strm->good())
    {
      *logstrm << "Cannot read " << helper.human_name() << ": \"" << filename << "\". Continuing with defaults." << std::endl;
      delete test_strm;
      if (! vm.count(paramname))
      {
        std::cerr << "Default " << helper.human_name() << " file not found. Aborting." << std::endl;
        return ELEC_TOOLS_BAD_INSTALLATION;
      }
      fstrm = new std::ifstream((config_path+helper.default_filename()).c_str());
    } else fstrm = test_strm;
    
    CommentSkipper skipper;
    std::istringstream strm(skipper.skip(*fstrm));
    
    std::string result = rules.readMe(strm);
    if( ! result.empty() )
    {
      std::cerr << "Reading " << helper.human_name() << " failed. Found unexpected: \"" << result << "\"." << std::endl;
      return ELEC_TOOLS_RULE_READ_FAILURE;
    }
    delete fstrm;
    return 0;
  }

  static int open_in_file(std::ifstream& strm, const std::string& filename)
  {
    strm.open(filename.c_str());
    if (!strm.good()) 
    {
      std::cerr << "Error opening input file: \"" << filename << "\"." << std::endl;
      return ELEC_TOOLS_NO_INPUT_FILE; 
    }
    return 0;
  }

private:
  typedef std::map<std::string, std::ofstream*> map_type;
  map_type ostreams; 
  int err;
  std::string human_error;
  
  void close_file(std::ostream* strm) 
  {
    if ( strm )
    {
      std::ofstream * t = dynamic_cast<std::ofstream *> (strm);
      if (t) t->close();
    }
  }
  
  static int read_mmol(std::ifstream& strm, GeneralMultipoledMolecule& mol)
  {
    CommentSkipper skipper;
    std::string istring = skipper.skip(strm);
    std::istringstream ifile(istring);
    if (!mol.readMe(ifile).empty())
    {
      std::cerr << "Failure reading input molecule." << std::endl;
      return ELEC_TOOLS_MOL_READ_FAILURE;
    }
    return 0;
  }
  
  static int read_obmol(std::ifstream& strm, OpenBabel::OBMol& mol, const std::string& filename, std::string& type)
  {
    OpenBabel::OBConversion conv(&strm);
    if (type.empty()) type = std::string("guess");
    std::transform(type.begin(), type.end(), type.begin(), tolower);
    if (type != "guess") {
      if(!conv.SetInFormat(type.c_str()))
      {
        std::cerr <<"Unsuppoted input molecule file format: \"" << type << "\"." << std::endl;
        return ELEC_TOOLS_BAD_MOL_FORMAT; 
      }
    } else {
      OpenBabel::OBFormat* fmt = OpenBabel::OBConversion::FormatFromExt(filename.c_str());
      if (!conv.SetInFormat(fmt))
      {
        std::cerr << "Input molecule format guess failed." << std::endl;
        return ELEC_TOOLS_BAD_MOL_FORMAT; 
      }
      type = std::string(fmt->TypeID());
    }
    
    if(!conv.Read(&mol)){
      std::cerr << "Failure reading input molecule." << std::endl;
      return ELEC_TOOLS_MOL_READ_FAILURE;
    }
    return 0;
  }
  
  void write_mmol(std::ostream& strm, GeneralMultipoledMolecule& mol)
  {    
    strm << mol << std::endl;
  }
  
  void write_obmol(std::ostream& strm, OpenBabel::OBMol& mol, const std::string& filename, std::string type = std::string("guess"))
  {
    OpenBabel::OBConversion conv;
    conv.SetOutStream(&strm);    
    std::transform(type.begin(), type.end(), type.begin(), tolower);
    if (type != "guess") {
      if(!conv.SetOutFormat(type.c_str()))
      {
        human_error = "Unsuppoted input molecule file format: \"" + type + "\".";
        err = ELEC_TOOLS_BAD_MOL_FORMAT; 
      }
    } else {
      if (!conv.SetOutFormat(OpenBabel::OBConversion::FormatFromExt(filename.c_str())))
      {
        human_error = "Input molecule format guess failed.";
        err = ELEC_TOOLS_BAD_MOL_FORMAT; 
      }
    }
    
    if (!conv.Write(&mol))
    {
      human_error = "Error writing final molecule.";
      err = ELEC_TOOLS_OUTPUT_FAILURE; 
    } 
  }
  
  template <class T>
  struct RuleHelper
  {
    std::string default_filename() const;
    std::string param_name() const;
    std::string human_name() const;
  };
  
};

template <>
struct IOTools::RuleHelper<MultipoleOrientRules>
{
  std::string default_filename() const { return "orient.rules"; }
  std::string param_name() const { return "orient-rules"; }
  std::string human_name() const { return "multipole orient rules"; }
};

template <>
struct IOTools::RuleHelper<MultipolePlacementRules>
{
  std::string default_filename() const { return "placement.rules"; }
  std::string param_name() const { return "placement-rules"; }
  std::string human_name() const { return "multipole placement rules"; }
};

#endif // IOTOOLS_H
