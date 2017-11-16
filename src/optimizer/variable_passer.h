/*************************************************************************
 * variable_passer.h - Convenience class to pass varibles to optimizers. *
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

#ifndef VARIABLE_PASSER_HEADER
#define VARIABLE_PASSER_HEADER

namespace Optimizers {

    template< typename T >
    class VariablePasserT
    {
    public:
      VariablePasserT(T& val) : v(val) { }

      VariablePasserT(const VariablePasserT& other) : v(other.v) { }

      VariablePasserT& operator=(const VariablePasserT& other) { v = other.v; return *this; }

      void operator()(T x) { v = x; }

      T& operator()() { return v; }

    private:
      T& v;
    };

} // namespace Optimizers

#define DECLARE_VARIABLE_PASSER \
  typedef VariablePasserT<real_type> VariablePasser;

#endif // VARIABLE_PASSER_HEADER
