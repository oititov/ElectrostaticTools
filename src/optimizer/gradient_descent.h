/*************************************************************************
 * gradient_descent.h - Gradient descent optimizer.                      *
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

#ifndef GRADIENT_DESCENT_HEADER
#define GRADIENT_DESCENT_HEADER

#include "abstract_optimizers.h"

namespace Optimizers {

    template< typename T >
    class GradientDescentOptimizerT : public GradientOptimizerBaseT< T >
    {
    public:
      DECLARE_REAL_TYPE
      DECLARE_DATA_VECTOR
      DECLARE_DATA_POINT

      GradientDescentOptimizerT() { setLambda(); }

      void setLambda(real_type l = 1.0) { lambda = l; }

      virtual bool optimize()
      {
        real_type error = HeuristicBase::eps * 10;
        DataPoint x;
        x.x = HeuristicBase::init;
        this->setAll(x.x);
        real_type f_prev;
        x.f = HeuristicBase::target();
        int counter = 0;
        do {
          f_prev = x.f;
          DataVector step(lambda * GradientBase::gradient());
          real_type step_length = step.squaredNorm();
          if (step_length >= GradientBase::max_step * GradientBase::max_step) step *= GradientBase::max_step / sqrt(step_length);
          x.x -= step;
          this->setAll(x.x);
          x.f = HeuristicBase::target();
        } while (( ++counter ) < HeuristicBase::max_iter && fabs(f_prev - x.f) > HeuristicBase::eps);
        HeuristicBase::best_point = x;
        if (counter >= HeuristicBase::max_iter) return false; else return true;
      }

    private:
      DECLARE_HEURISTIC_BASE
      DECLARE_GRADIENT_BASE
      real_type lambda;
    };

    typedef GradientDescentOptimizerT< double > GradientDescentOptimizer;

} // namespace Optimizers

#endif // GRADIENT_DESCENT_HEADER
