/*************************************************************************
 * optimizer_defines.h - convenience defines for optimizers.             *
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

#ifndef OPTIMIZER_DEFINES_OPTIMIZER_HEADER
#define OPTIMIZER_DEFINES_OPTIMIZER_HEADER

#include <boost/function.hpp>

#define DECLARE_REAL_TYPE \
  typedef T real_type;

#define DECLARE_FUNC_TYPES \
  typedef boost::function< real_type (void) > target_function; \
  typedef boost::function< real_type (void) > get_function; \
  typedef boost::function< real_type (void) > derivative_function; \
  typedef boost::function< void (real_type) > set_function;

#endif // OPTIMIZER_DEFINES_OPTIMIZER_HEADER
