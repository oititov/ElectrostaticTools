/*************************************************************************
 * elec_tools_errors.h - Electrostatic Tools error codes.                *
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

#ifndef ELEC_TOOLS_ERRORS_H
#define ELEC_TOOLS_ERRORS_H

#define ELEC_TOOLS_UNKNOWN_ERROR -1
#define ELEC_TOOLS_BAD_OPTION 1
#define ELEC_TOOLS_NO_INPUT_FILE 2
#define ELEC_TOOLS_NO_GRID_FILE 3
#define ELEC_TOOLS_OUTPUT_EXISTS 4
#define ELEC_TOOLS_UNWRITABLE_OUTPUT 5
#define ELEC_TOOLS_BAD_MOL_FORMAT 6
#define ELEC_TOOLS_MOL_READ_FAILURE 7
#define ELEC_TOOLS_RULE_READ_FAILURE 8
#define ELEC_TOOLS_FITTING_FAILURE 9
#define ELEC_TOOLS_OUTPUT_FAILURE 10
#define ELEC_TOOLS_BAD_INSTALLATION 11

#endif // ELEC_TOOLS_ERRORS