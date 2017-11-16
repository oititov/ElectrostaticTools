/*************************************************************************
 * grid.h - Handler of reference electrostatic grids.                    *
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

#ifndef GRID_HEADER
#define GRID_HEADER

#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>

#include "comment_skipper.h"
#include "mult_tools_defines.h"
#include "vector3.h"

class GridNode
{
public:
  typedef double real_type;
  typedef Vector3 point_type;

  GridNode(const point_type p = point_type(), const real_type v = real_type()) : point(p), value(v) { }
  
  void printMe(std::ostream& strm) const
  {
    strm << "Point: " << point << " Field: " << value;
  }
  
  std::string readMe(std::istream& strm, std::string header = "")
  {
    point_type p;
    real_type f;
    if (header.empty()) strm >> header;
    if (header != "Point:") return header;
    strm >> p;
    header = "";
    strm >> header;
    if (header != "Field:") throw std::domain_error("Reading grid point: expected \"Field:\" found " + header + "\"");
    strm >> f;
    point = p;
    value = f;
    return std::string();
  }
  
  point_type point;
  real_type value;
};

#ifndef SWIG
DECLARE_IO_STREAM_OPERATORS(GridNode)
#endif // SWIG

#ifdef SWIG
__STR__(GridNode);
#endif // SWIG

class Grid
{
public:
  const Vector3& lattice(const int i) const { return grid_[i].point; }

  const double field(const int i) const { return grid_[i].value; }

  void addPoint(const Vector3& position, const double value) { grid_.push_back(GridNode(position, value)); }

  void modifyField(const int i, const double value) { grid_[i].value += value; }

  const GridNode& grid(const int i) const { return grid_[i]; }

  double charge() const { return q_tot; }

  void charge(const double q) { q_tot += q; }

  int size() const { return grid_.size(); }
  
  void toESP(std::ostream& strm, const std::vector< Vector3 > *coords = 0) const
  // note that without coordinates this will be incomplete ESP file
  {
    strm.setf(std::ios::right);
    strm.setf(std::ios::scientific);
    strm << std::setw(4) << (coords ? coords->size() : 0) << 
            std::setw(6) << grid_.size() << std::setw(5) << static_cast<int>(charge()) << std::endl;
    if (coords)
      for (int i = 0; i < coords->size(); ++i)
      {
        Vector3 p = ( *coords )[i];
        strm << "                " <<
            std::setw(16) << std::setprecision(7) << p.x() <<
            std::setw(16) << std::setprecision(7) << p.y() <<
            std::setw(16) << std::setprecision(7) << p.z() << std::endl;
      }
    for (int i = 0; i < grid_.size(); ++i)
    {
      Vector3 p = lattice(i);
      strm << std::setw(16) << std::setprecision(7) << field(i) <<
          std::setw(16) << std::setprecision(7) << p.x() <<
          std::setw(16) << std::setprecision(7) << p.y() <<
          std::setw(16) << std::setprecision(7) << p.z() << std::endl;
    }
    strm.unsetf(std::ios::right);
    strm.unsetf(std::ios::scientific);
  }
  
  bool fromESP(const std::string& filename )
  {
    std::ifstream input(filename.c_str());
    if (!input.good()) return false;
    CommentSkipper skipper;
    std::istringstream ifile(skipper.skip(input));
    fromESP(ifile);
    return true;
  }

  void fromESP(std::istream& strm, std::vector< Vector3 > *coords = 0)
  {
    // reads grid in esp format
    // we're reading positions since spaces between 
    // numbers are not guaranteed by esp format
    int natoms, npoints;
    double x, y, z, v;
    char line[256];
    char buff[17];
    buff[16] = '\0';
    strm.getline(line, 256);
    memcpy(buff, line, 4);
    buff[4] = '\0';
    natoms = atoi(buff);
    memcpy(buff, line + 4, 6);
    buff[6] = '\0';
    npoints = atoi(buff);
    memcpy(buff, line + 10, 5);
    buff[5] = '\0';
    q_tot = atoi(buff);
    grid_.reserve(npoints);
    
    if (!coords)
      for (int i = 0; i < natoms; ++i) strm.getline(line, 256);
    else
    {
      coords->clear();
      coords->reserve(natoms);
      for (int i = 0; i < natoms; ++i)
      {
        strm.getline(line, 256);
        memcpy(buff, line + 16, 16);
        x = atof(buff);
        memcpy(buff, line + 32, 16);
        y = atof(buff);
        memcpy(buff, line + 48, 16);
        z = atof(buff);
        coords->push_back(Vector3(x, y, z));
      }
    }
    for (int i = 0; i < npoints; ++i)
    {
      strm.getline(line, 256);
      memcpy(buff, line + 1, 16);
      v = atof(buff);
      memcpy(buff, line + 16, 16);
      x = atof(buff);
      memcpy(buff, line + 32, 16);
      y = atof(buff);
      memcpy(buff, line + 48, 16);
      z = atof(buff);
      grid_.push_back(GridNode(Vector3(x, y, z), v));
    }
  }
  
  void printMe(std::ostream& strm) const
  {
    strm << "Grid:" << std::endl;
    strm << "Charge: " << q_tot << std::endl;
    for (int i = 0; i < grid_.size(); ++i) strm << grid_[i] << std::endl;
  }
  
  bool readMe(const std::string& filename )
  {
    std::ifstream input(filename.c_str());
    if (!input.good()) return false;
    CommentSkipper skipper;
    std::istringstream ifile(skipper.skip(input));
    if (!readMe(ifile).empty()) return false;
  }

  std::string readMe(std::istream& strm, std::string header = "")
  {
    if (header.empty()) strm >> header;
    if (header != "Grid:") return header;
    if (strm.good())
    {
      double q;
      std::string buff;
      strm >> buff;
      if (buff.compare("Charge:")) throw std::domain_error("Reading grid: expected \"Charge:\" but found \"" + buff + "\"");
      strm >> q;
      q_tot = q;
    }
    grid_.clear();
    std::string buff;
    while (strm.good())
    {
      GridNode n;
      buff = n.readMe(strm);
      if (!buff.empty()) return buff;
      grid_.push_back(n);
    }
    return std::string();
  }

private:
  std::vector< GridNode > grid_;
  double q_tot;
};

#ifndef SWIG
DECLARE_IO_STREAM_OPERATORS(Grid)
#endif // SWIG

#ifdef SWIG
__STR__(Grid);

%extend Grid{
  const std::string toEspString() {
    std::ostringstream out;
    $self->toESP(out);
    return out.str();
  }
}
#endif // SWIG

#endif // GRID_HEADER
