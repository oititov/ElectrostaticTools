/*************************************************************************
 * genetic.h - Genetic optimizer definitions.                            *
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

#ifndef GENETIC_HEADER
#define GENETIC_HEADER

#include <boost/random.hpp>

#include "abstract_optimizers.h"

namespace Optimizers {

    template< typename T, typename D >
    class DistributionGeneratorT
    {
    public:
      D operator()(T a, T b) const { }
    };

    template< typename T >
    class DistributionGeneratorT< T, boost::uniform_real< T > >
    {
    public:
      boost::uniform_real< T > operator()(T a, T b) const
      {
        return boost::uniform_real< T >(a - b, a + b);
      }
    };

    template< typename T >
    class DistributionGeneratorT< T, boost::normal_distribution< T > >
    {
    public:
      boost::normal_distribution< T > operator()(T a, T b) const
      {
        return boost::normal_distribution< T >(a, b / 3.0);
      }
    };

    template< typename T, typename D = boost::uniform_real< T >, typename G = boost::mt19937 >
    class GeneticOptimizerT : public GlobalOptimizerBaseT< T >, public HeuristicOptimizerT< T >
    {
    public:
      DECLARE_REAL_TYPE
      DECLARE_DATA_VECTOR
      DECLARE_DATA_POINT

      GeneticOptimizerT()
      {
        setPopulation();
        setSelectionPartition();
        setMutationP();
        setSeed(time(0));
        setPreserveSelected();
        setSelection();
        setOnlySelectedGenes();
      }
  
      enum Selection
      {
        FitnessProportionate, Best
      };

      virtual bool optimize()
      {
        // initialization
        std::vector< DataPoint > population(numPopulation);
        int numDim = HeuristicBase::init.size();
        DistributionGeneratorT< real_type, distribution_type > dist_gen;
        std::vector< distribution_type > distributions(numDim);
        for (int i = 0; i < numDim; ++i) distributions[i] = dist_gen(HeuristicBase::init[i], HeuristicBase::steps[i]);
        for (int i = 0; i < numPopulation; ++i) {
          population[i].x = DataVector(numDim);
          for (int j = 0; j < numDim; ++j) population[i].x[j] = distributions[j](generator);
          this->setAll(population[i].x);
          population[i].f = HeuristicBase::target();
        }

        uniform_distribution_type zero_one(0, 1);

        int counter = 0;
        while (( counter++ ) < HeuristicBase::max_iter) {
          // selection
          std::vector< bool > selected;
          switch (selection_algo) {
            case FitnessProportionate:
              selectionFitnessProportionate(population, selected);
              break;
            case Best:
              selectionBest(population, selected);
              break;
          }
          // crossover
          std::vector< DataPoint > new_generation;
          std::vector< bool > modified;
          crossoverUniform(population, new_generation, selected, modified);
          // mutation

          for (int i = 0; i < numPopulation; ++i)
            if (preserve_selected && selected[i]) continue; else
              for (int j = 0; j < numDim; ++j) {
                real_type dice = zero_one(generator);
                if (dice < mutationP) {
                  new_generation[i].x[j] = distributions[j](generator);
                  modified[i] = true;
                }
              }
          // finalization
          for (int i = 0; i < numPopulation; ++i)
            if (modified[i]) {
              this->setAll(new_generation[i].x);
              new_generation[i].f = HeuristicBase::target();
            }
          population = new_generation;
        }
        std::sort(population.begin(), population.end());
        GlobalBase::points = population;
        HeuristicBase::best_point = population.front();

        if (counter >= HeuristicBase::max_iter) return false; else return true;
      }

      void setSeed(int s)
      {
        seed = s;
        generator = generator_type(seed);
      }

      void setPopulation(int num = 100) { numPopulation = num; }

      void setSelectionPartition(real_type R = 0.5) { selectionPartition = R; }

      void setMutationP(real_type p = 0.01) { mutationP = p; }

      void setPreserveSelected(bool flag = true) { preserve_selected = flag; }

      void setOnlySelectedGenes(bool flag = false) { only_selected_genes = flag; }

      void setSelection(Selection s = Best) { selection_algo = s; }

    private:
      DECLARE_HEURISTIC_BASE
      DECLARE_GLOBAL_BASE
      typedef D distribution_type;
      typedef G generator_type;
      typedef boost::uniform_int< int > int_distribution_type;
      typedef boost::uniform_real< real_type > uniform_distribution_type;

      void selectionFitnessProportionate(std::vector< DataPoint >& population, std::vector< bool >& selected)
      {
        real_type f_max = population[0].f;
        for (int i = 0; i < numPopulation; ++i) if (population[i].f > f_max) f_max = population[i].f;
        real_type fSum = f_max * numPopulation;
        for (int i = 0; i < numPopulation; ++i) fSum -= population[i].f;
        std::vector< real_type > accum_fitness(numPopulation);
        accum_fitness[0] = ( f_max - population[0].f ) / fSum;
        for (int i = 1; i < numPopulation; ++i) accum_fitness[i] = accum_fitness[i - 1] + ( f_max - population[i].f ) / fSum;
    
        int numSelected = static_cast<int>(numPopulation * selectionPartition);
        int selectedCounter = 0;
        selected = std::vector< bool >(numPopulation, false);
        uniform_distribution_type zero_one(0, 1);
        while (selectedCounter < numSelected) {
          real_type R = zero_one(generator);
          for (int i = 0; i < numPopulation; ++i) {
            if (accum_fitness[i] > R) {
              if (!selected[i]) {
                selected[i] = true;
                selectedCounter++;
              }
              break;
            }
          }
        }
      }

      void selectionBest(std::vector< DataPoint >& population, std::vector< bool >& selected)
      {
        int numSelected = static_cast<int>(numPopulation * selectionPartition);
        selected = std::vector< bool >(numPopulation, false);
        std::partial_sort(population.begin(), population.begin() + numSelected, population.end());
        for (int i = 0; i < numSelected; ++i) selected[i] = true;
      }

      void crossoverUniform(const std::vector< DataPoint >& oldPopulation, std::vector< DataPoint >& newPopulation, const std::vector< bool >& selected, std::vector< bool >& modified)
      {
        newPopulation = std::vector< DataPoint >();
        newPopulation.reserve(numPopulation);
        modified = std::vector< bool >(numPopulation, false);
        uniform_distribution_type zero_one(0, 1);
        int_distribution_type population_distr(0, numPopulation - 1);
        int numDim = oldPopulation.front().x.size();
        for (int i = 0; i < numPopulation; ++i) {
          if (selected[i]) {
            newPopulation.push_back(oldPopulation[i]);
            modified[i] = false;
          }
          else {
            int parents[2];
            if (!only_selected_genes) {
              parents[0] = population_distr(generator);
              while (( parents[1] = population_distr(generator) ) == parents[0]);
            } else {
              do { parents[0] = population_distr(generator); } while (!selected[parents[0]]);
              do { parents[1] = population_distr(generator); } while (!selected[parents[1]] && parents[0] == parents[1]);
            }
            DataPoint child;
            child.x = DataVector(numDim);
            for (int j = 0; j < numDim; ++j)
              child.x[j] = ( zero_one(generator) > 0.5 ) ? oldPopulation[parents[0]].x[j] : oldPopulation[parents[1]].x[j];
            newPopulation.push_back(child);
            modified[i] = true;
          }
        }
      }

      int numPopulation;
      int seed;
      generator_type generator;
      real_type selectionPartition;
      real_type mutationP;
      bool preserve_selected;
      bool only_selected_genes;
      Selection selection_algo;
    };

    typedef GeneticOptimizerT< double > GeneticOptimizer;

} // namespace Optimizers

#endif // GENETIC_HEADER
