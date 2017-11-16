/*************************************************************************
 * numerical_derivative.h - Class for calculation of                     *
 *                          numerical derivatives.                       *
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

#ifndef NUMERICAL_DERIVATIVE_HEADER
#define NUMERICAL_DERIVATIVE_HEADER

#include <boost/function.hpp>

#include "optimizer_defines.h"

namespace Optimizers {

    template< typename T >
    class NumericalDerivativeT
    {
    public:
      DECLARE_REAL_TYPE
      DECLARE_FUNC_TYPES
    public:
      NumericalDerivativeT(const target_function& f, const set_function& s, const get_function& g, real_type dx = 1.0e-4)
      {
        target = f;
        get = g;
        set = s;
        step = dx;
      }

      real_type operator()() const
      {
        real_type x = get();
        real_type result = 0.0;
        set(x + step);
        result += target();
        set(x - step);
        result -= target();
        result /= 2.0 * step;
        set(x);
        return result;
      }

    private:

      target_function target;
      get_function get;
      set_function set;
      real_type step;
    };

    typedef NumericalDerivativeT< double > NumericalDerivative;

} // namespace Optimizers

#define DECLARE_NUMERICAL_DERIVATIVE \
  typedef NumericalDerivativeT<real_type> NumericalDerivative;

#endif // NUMERICAL_DERIVATIVE_HEADER

