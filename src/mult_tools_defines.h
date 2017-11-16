/*************************************************************************
 * mult_tools_defines.h - Convenience macro                              *
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

#ifndef COMMON_DEFINES_MULT_TOOLS_HEADER
#define COMMON_DEFINES_MULT_TOOLS_HEADER

#define DECLARE_IO_STREAM_OPERATORS(t) \
inline std::ostream& operator<<(std::ostream& strm, const t & x) \
{\
  x.printMe(strm); \
  return strm; \
} \
\
inline std::istream& operator>>(std::istream& strm, t & x) \
{ \
  x.readMe(strm); \
  return strm; \
}

#endif // COMMON_DEFINES_MULT_TOOLS_HEADER
