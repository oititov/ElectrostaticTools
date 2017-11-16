/*************************************************************************
 * abstract_optimizers.h - Optimizer interfaces.                         *
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

#ifndef ABSTRACT_OPTIMIZERS_HEADER
#define ABSTRACT_OPTIMIZERS_HEADER

#include <cmath>
#include <vector>

#include "optimizer_defines.h"
#include "data_point.h"
#include "variable_passer.h"
#include "numerical_derivative.h"

#define COMMON_DECLARATIONS \
DECLARE_REAL_TYPE \
DECLARE_DATA_VECTOR \
DECLARE_DATA_POINT \
DECLARE_FUNC_TYPES

namespace Optimizers {

    template< typename T >
    class OptimizerT
    {
    public:
      COMMON_DECLARATIONS

      virtual void setTargetFunction(const target_function& f) = 0;

      virtual void addParameter(const set_function& s, real_type start, real_type step) = 0;

      virtual void addParameter(real_type& s, real_type step) = 0;

      virtual void setEpsilon(real_type e) = 0;

      virtual void setMaxIter(real_type i) = 0;

      virtual bool optimize() = 0;

      virtual DataPoint result() const = 0;

      virtual void setAll(const DataVector& x) const = 0;

      virtual ~OptimizerT() { }
    };

    typedef OptimizerT< double > Optimizer;

    template< typename T >
    class HeuristicOptimizerT : public OptimizerT< T >
    {
    public:
      COMMON_DECLARATIONS
    private:
      DECLARE_VARIABLE_PASSER
    public:

      HeuristicOptimizerT()
      {
        HeuristicOptimizerT< T >::setEpsilon();
        HeuristicOptimizerT< T >::setMaxIter();
      }

      virtual void setTargetFunction(const target_function& f) { target = f; }

      virtual void addParameter(const set_function& s, real_type start, real_type step)
      {
        set.push_back(s);
        init.conservativeResize(init.size() + 1);
        init[init.size() - 1] = start;
        steps.conservativeResize(steps.size() + 1);
        steps[steps.size() - 1] = step;
      }

      virtual void addParameter(real_type& s, real_type step)
      {
        vars.push_back(VariablePasser(s));
        addParameter(vars.back(), s, step);
      }

      virtual void setEpsilon(real_type e = 1.0e-6) { eps = e; }

      virtual void setMaxIter(real_type i = 100) { max_iter = i; }

      virtual void setAll(const DataVector& v) const
      {
        for (int i = 0; i < set.size(); ++i) set[i](v[i]);
      }

      virtual DataPoint result() const
      {
        return best_point;
      }

    protected:
      target_function target;
      DataVector init;
      DataVector steps;
      std::vector< set_function > set;
      std::vector< VariablePasser > vars;
      real_type eps;
      int max_iter;
      DataPoint best_point;
    };

#define DECLARE_HEURISTIC_BASE \
  typedef HeuristicOptimizerT<real_type> HeuristicBase;

    template< typename T >
    class GradientOptimizerT : public HeuristicOptimizerT< T >
    {
    public:
      DECLARE_REAL_TYPE
      DECLARE_FUNC_TYPES

      virtual void addParameter(const set_function& s, const derivative_function& d, real_type start) = 0;

      virtual void addParameter(real_type& s, const derivative_function& d) = 0;

      virtual void addParameter(const set_function& s, const get_function& g) = 0;

      virtual void addParameter(real_type& var) = 0;

      virtual void setMaxStep(real_type l) = 0;
    };

    typedef GradientOptimizerT< double > GradientOptimizer;

    template< typename T >
    class GradientOptimizerBaseT : public GradientOptimizerT< T >
    {
    public:
      COMMON_DECLARATIONS

    private:
      DECLARE_VARIABLE_PASSER
      DECLARE_NUMERICAL_DERIVATIVE
      DECLARE_HEURISTIC_BASE

    public:
      GradientOptimizerBaseT()
      {
        GradientOptimizerBaseT< T >::setEpsilon();
        GradientOptimizerBaseT< T >::setMaxIter();
        GradientOptimizerBaseT< T >::setMaxStep();
      }

      virtual void addParameter(const set_function& s, const derivative_function& d, real_type start)
      {
        HeuristicBase::set.push_back(s);
        derivative.push_back(d);
        HeuristicBase::init.conservativeResize(HeuristicBase::init.size() + 1);
        HeuristicBase::init[HeuristicBase::init.size() - 1] = start;
      }

      virtual void addParameter(real_type& s, const derivative_function& d)
      {
        HeuristicBase::vars.push_back(VariablePasser(s));
        addParameter(HeuristicBase::vars.back(), d, s);
      }

      virtual void addParameter(const set_function& s, const get_function& g)
      {
        ders.push_back(NumericalDerivative(HeuristicBase::target, s, g));
        addParameter(s, ders.back(), g());
      }

      virtual void addParameter(real_type& var)
      {
        HeuristicBase::vars.push_back(VariablePasser(var));
        addParameter(HeuristicBase::vars.back(), HeuristicBase::vars.back());
      }

      virtual void setMaxStep(real_type l = 1.0) { max_step = l; }

    protected:
      DataVector gradient() const
      {
        DataVector result(derivative.size());
        for (int i = 0; i < derivative.size(); ++i) result[i] = derivative[i]();
        return result;
      }

      std::vector< derivative_function > derivative;
      std::vector< NumericalDerivative > ders;
      real_type max_step;
    };

#define DECLARE_GRADIENT_BASE \
  typedef GradientOptimizerBaseT<real_type> GradientBase;

    template< typename T >
    class GlobalOptimizerT
    {
    public:
      DECLARE_REAL_TYPE
      DECLARE_DATA_POINT

      virtual std::vector< DataPoint > results() const = 0;
    };

    typedef GlobalOptimizerT< double > GlobalOptimizer;

    template< typename T >
    class GlobalOptimizerBaseT : public GlobalOptimizerT< T >
    {
    public:
      DECLARE_REAL_TYPE
      DECLARE_DATA_POINT

      virtual std::vector< DataPoint > results() const
      {
        return points;
      }

    protected:
      std::vector< DataPoint > points;
    };

#define DECLARE_GLOBAL_BASE \
  typedef GlobalOptimizerBaseT<real_type> GlobalBase;

} // namespace Optimizers

#endif // ABSTRACT_OPTIMIZERS_HEADER
