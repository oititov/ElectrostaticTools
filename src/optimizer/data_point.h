/*************************************************************************
 * data_point.h - Convenience functions for optimizers.                  *
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

#ifndef DATA_POINT_HEADER
#define DATA_POINT_HEADER

#include <Eigen/Dense>

#include "optimizer_defines.h"

namespace Optimizers {

#define DECLARE_DATA_VECTOR \
  typedef Eigen::Matrix<real_type, Eigen::Dynamic, 1> DataVector;

    template< typename T >
    struct DataPointT
    {
      DECLARE_REAL_TYPE
      DECLARE_DATA_VECTOR
  
      DataVector x;
      real_type f;

      bool operator<(const DataPointT& other) const { return f < other.f; }

      bool operator>(const DataPointT& other) const { return f > other.f; }

      bool operator<=(const DataPointT& other) const { return f <= other.f; }

      bool operator>=(const DataPointT& other) const { return f >= other.f; }
    };

#define DECLARE_DATA_POINT \
  typedef DataPointT<real_type> DataPoint;

} // namespace Optimizers

#endif // DATA_POINT_HEADER
