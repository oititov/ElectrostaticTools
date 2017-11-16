/*************************************************************************
 * elec_tools.i - script bindings interface file.                        *
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

%include <stl.i>
%include <carrays.i>
%include <cpointer.i>

%define extend_toString(class_name)
%extend class_name{
    const std::string toString() {
      std::ostringstream out;
      $self->printMe(out);
      return out.str();
    }
}
%enddef // extend_toString

#ifdef SWIGPYTHON
%define __STR__(class_name)
// this feature does not work without -builtin flag, 
// however OpenBabel has problems with -builtin
//%feature("python:slot", "tp_str", functype="reprfunc") class_name::toString();

extend_toString(class_name)
%enddef
#endif // SWIGPYTHON

#ifdef SWIGPERL
%define __STR__(class_name)
extend_toString(class_name)
%enddef
#endif // SWIGPERL

%module ElectrostaticTools
%{
typedef double doubleArray; // this is not generated by SWIG 3 for some reason.
namespace OpenBabel 
{ 
// declarations to make bindings compilable

  class OBReaction; 
  class OBAtomClassData;
  class OBStereoBase;
  class OBTetraPlanarStereo;
  class OBCisTransStereo;
  class OBTetraNonPlanarStereo;
  class OBTetrahedralStereo;
  class OBSquarePlanarStereo;
  class OBRotamerList;
  class AliasData;
  class OBClassData;
  class OBGridData;
};

#include "grid.h"
#include "molecule.h"
%}

#ifdef SWIGPYTHON
%import "openbabel-python.i"
#endif // SWIGPYTHON
#ifdef  SWIGPERL
%import "openbabel-perl.i"
#endif // SWIGPERL

%include "vector3.h"
%include "matrix3.h"
%include "grid.h"
%include "multipoles.h"
%include "molecule.h"
%include "multipole_orient_rules.h"

