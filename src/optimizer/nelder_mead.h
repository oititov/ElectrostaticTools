/*************************************************************************
 * nelder_mead.h - Nelder-Mead optimizer.                                *
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

#ifndef NELDER_MEAD_OPTIMIZER_HEADER
#define NELDER_MEAD_OPTIMIZER_HEADER

#include "abstract_optimizers.h"

namespace Optimizers {

    template< typename T >
    class NelderMeadOptimizerT : public HeuristicOptimizerT< T >
    {
    public:
      DECLARE_REAL_TYPE
      DECLARE_DATA_VECTOR
      DECLARE_DATA_POINT

      NelderMeadOptimizerT()
      {
        setReflection();
        setExpansion();
        setContraction();
        setShrink();
      }
  
      virtual bool optimize()
      {
        // preparation
        DataPoint init_point;
        init_point.x = HeuristicBase::init;
        int n = HeuristicBase::init.size();
        std::vector< DataPoint > simplex(n + 1, init_point);
        for (int i = 1; i < simplex.size(); ++i) simplex[i].x[i - 1] += HeuristicBase::steps[i - 1];
        for (int i = 0; i < simplex.size(); ++i) {
          this->setAll(simplex[i].x);
          simplex[i].f = HeuristicBase::target();
        }
        real_type radius = 9999999;
        int counter = 0;
        while (counter++ < HeuristicBase::max_iter && radius > HeuristicBase::eps) {
          // sorting
          int h = 0, g = 0, l = 0;
          for (int i = 0; i < simplex.size(); ++i) {
            if (simplex[i].f >= simplex[h].f) {
              g = h;
              h = i;
            } else if (simplex[i].f >= simplex[g].f) g = i; else if (simplex[i].f <= simplex[l].f) l = i;
          }

          // simplex center
          DataVector center = DataVector::Zero(n);
          for (int i = 0; i < simplex.size(); ++i) if (i != h) center += 1.0 / n * simplex[i].x;

          // reflection
          DataPoint reflected;
          reflected.x = ( 1.0 + reflection ) * center - reflection * simplex[h].x;
          this->setAll(reflected.x);
          reflected.f = HeuristicBase::target();
          bool skip_to_end = false;

          // expansion
          if (reflected.f < simplex[l].f) {
            DataPoint expanded;
            expanded.x = ( 1.0 - expansion ) * center + expansion * reflected.x;
            this->setAll(expanded.x);
            expanded.f = HeuristicBase::target();
            if (expanded.f < simplex[l].f) simplex[l] = expanded;
            else simplex[l] = reflected;
            skip_to_end = true;
          } else if (reflected.f < simplex[g].f) {
            simplex[h] = reflected;
            skip_to_end = true;
          } else if (reflected.f < simplex[h].f) {
            DataPoint tmp(reflected);
            reflected = simplex[h];
            simplex[h] = tmp;
          }

          // contraction
          if (!skip_to_end) {
            DataPoint contracted;
            contracted.x = contraction * simplex[h].x + ( 1.0 - contraction ) * center;
            this->setAll(contracted.x);
            contracted.f = HeuristicBase::target();
            if (contracted.f < simplex[h].f) {
              simplex[h] = contracted;
            } else {
              // shrink
              for (int i = 0; i < simplex.size(); ++i)
                if (i != l) {
                  simplex[i].x = simplex[l].x + shrink * ( simplex[i].x - simplex[l].x );
                  this->setAll(simplex[i].x);
                  simplex[i].f = HeuristicBase::target();
                }
            }
          }
          // convergence test
          real_type dist;
          center = DataVector::Zero(n);
          radius = 0.0;
          for (int i = 0; i < simplex.size(); ++i) center += 1.0 / ( n + 1 ) * simplex[i].x;
          for (int i = 0; i < simplex.size(); ++i) if (( dist = ( center - simplex[i].x ).norm() ) > radius) radius = dist;
        }

        // find minimum
        int l = 0;
        for (int i = 0; i < simplex.size(); ++i) {
          if (simplex[i].f <= simplex[l].f) l = i;
        }
        HeuristicBase::best_point = simplex[l];

        if (counter >= HeuristicBase::max_iter) return false; else return true;
      }

      void setReflection(real_type c = 1.0) { reflection = c; }

      void setExpansion(real_type c = 2.0) { expansion = c; }

      void setContraction(real_type c = 0.5) { contraction = c; }

      void setShrink(real_type c = 0.5) { shrink = c; }

    private:
      DECLARE_HEURISTIC_BASE

      real_type reflection;
      real_type expansion;
      real_type contraction;
      real_type shrink;
    };

    typedef NelderMeadOptimizerT< double > NelderMeadOptimizer;

} // namespace Optimizers

#endif // NELDER_MEAD_OPTIMIZER_HEADER
