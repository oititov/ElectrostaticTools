/*************************************************************************
 * multipoles.h - Multipole classes.                                     *
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

#ifndef MULTIPOLES_H
#define MULTIPOLES_H

#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <stdexcept>

#include <openbabel/generic.h>

#include "mult_tools_defines.h"
#include "vector3.h"
#include "matrix3.h"

enum MULTIPOLES
{
  NONE = 0, 
  MONOPOLE = 1 << 0, 
  DIPOLE = 1 << 1, 
  QUADRUPOLE = 1 << 2,
  DIPOLE_DIRECTED = 1 << 3,
  QUADRUPOLE_X = 1 << 4,
  QUADRUPOLE_Y = 1 << 5,
  QUADRUPOLE_Z = 1 << 6
};

const double bohr = 0.52917720859;

class Multipole
{
public:
  virtual void rotate(const Matrix3& rot) { }

  virtual void printMe(std::ostream& strm = std::cout) const = 0;

  virtual std::string readMe(std::istream& strm, std::string header = "") = 0 ;

};

#ifndef SWIG
DECLARE_IO_STREAM_OPERATORS(Multipole)
#endif // SWIG

// all energy() methods accepts r as vector TO multipole, which method was called

class Dipole;

class Quadrupole;

class Monopole : public Multipole
{
public:
  Monopole(const double charge = 0.0) : _q(charge)
  {
  }

  Monopole(const Monopole& other) : _q(other._q)
  {
  }

  double value() const { return _q; }

  void setValue(const double v) { _q = v; }

  double energy(const Monopole& other, const double dist) const
  {
    return value() * other.value() / dist;
  }

  double energy(const Dipole& other, const Vector3& r, const double dist_3) const;

  double energy(const Quadrupole& other, const Vector3& r, const double dist_5) const;

  double potential(const double dist) const
  {
    return energy(Monopole(1.0), dist);
  }

  Vector3 field(const Vector3& r, const double dist_3) const
  {
    return r * _q / dist_3;
  }

  Vector3 field(const Vector3& r) const
  {
    double l2 = r.length_squared();
    return r * _q / l2 / sqrt(l2);
  }

  bool operator==(const Monopole& other) const
  {
    return _q == other._q;
  }

  bool operator!=(const Monopole& other) const
  {
    return _q != other._q;
  }

  void printMe(std::ostream& strm = std::cout) const
  {
    strm << "Monopole: " << std::setprecision(5) << _q;
  }

  std::string readMe(std::istream& strm, std::string header = "")
  {
    double t;
    if (header.empty())
    {
      strm >> header;
    }
    if (header != "Monopole:") return header;
    strm >> t;
    _q = t;
    return "";
  }

private:
  double _q;
};

#ifdef SWIG
__STR__(Monopole);
#endif // SWIG

class Dipole : public Multipole
{
public:
  Dipole(const Vector3& dip = Vector3()) : d(dip)
  {
  }

  Dipole(const Dipole& other) : d(other.d)
  {
  }

  Dipole(const double x, const double y, const double z) : d(Vector3(x, y, z)) { }

  Dipole(const double z) : d(Vector3(0.0, 0.0, z)) { }
  
  const Vector3& value() const { return d; }

  void setValue(const double x, const double y, const double z)
  {
    d = Vector3(x, y, z);
  }

  void setValue(const double z)
  {
    setValue(0.0, 0.0, z);
  }

  double energy(const Dipole& other, const Vector3& r, const double dist_2, const double dist_5) const
  {
    return ( dist_2 * d * ( other.value() ) - 3.0 * d.dot(r) * other.value().dot(r) ) / dist_5;
  }

  double energy(const Monopole& other, const Vector3& r, const double dist_3) const
  {
    return -other.value() * ( d * r ) / dist_3;
  }

  double energy(const Quadrupole& other, const Vector3& r, const double dist_5, const double dist_7) const;

  double potential(const Vector3& r, const double dist_3) const
  {
    return energy(Monopole(1.0), r, dist_3);
  }

  Vector3 field(const Vector3& r, const double dist_3, const double dist_5) const
  {
    return 3.0 * ( d * r ) * r / dist_5 - d / dist_3;
  }

  Vector3 field(const Vector3& r) const
  {
    double l2 = r.length_squared();
    double l3 = l2 * sqrt(l2);
    return field(r, l3, l3 * l2);
  }

  void rotate(const Matrix3& rot) { d = rot * d; };

  bool operator==(const Dipole& other) const
  {
    return d == other.d;
  }

  bool operator!=(const Dipole& other) const
  {
    return d != other.d;
  }

  void printMe(std::ostream& strm = std::cout) const
  {
    strm << "Dipole: " << d;
  }

  std::string readMe(std::istream& strm, std::string header = "")
  {
    Vector3 t;
    if (header.empty())
    {
      strm >> header;
    }
    if (header != "Dipole:") return header;
    strm >> t;
    d = t;
    return header;
  }

private:
  Vector3 d;
};

#ifdef SWIG
__STR__(Dipole);
#endif // SWIG

class Quadrupole : public Multipole
{
public:
  Quadrupole(const Matrix3& quad = Matrix3()) : q(quad)
  {
  }

  Quadrupole(const Quadrupole& other) : q(other.q)
  {
  }

  Quadrupole(const double xx, const double yy, const double zz, const double xy = 0.0, const double xz = 0.0, const double yz = 0.0)
      : q(Matrix3(xx, xy, xz, xy, yy, yz, xz, yz, zz)) { }

  Quadrupole(const double zz) : q(Matrix3(-zz/2.0, 0.0, 0.0, 0.0, -zz/2.0, 0.0, 0.0, 0.0, zz)) { }
  
  const Matrix3& value() const { return q; }

  void setComponent(const int i, const int j, const double v) { q[i][j] = v; }

  void setValue(const double xx, const double yy, const double zz, const double xy = 0.0, const double xz = 0.0, const double yz = 0.0)
  {
    q = Matrix3(xx, xy, xz, xy, yy, yz, xz, yz, zz);
  }

  double energy(const Quadrupole& other, const Vector3& r, const double dist_5, const double dist_7, const double dist_9) const
  {
    return 
    (2.0 / 3.0) * ( q * other.value() ).trace() / dist_5
    //- (20.0 / 3.0) * ( Matrix3(q * r, other.value() * r) ).trace() / dist_7
    - (20.0 / 3.0) * ( r * other.value() * q * r ) / dist_7
    //+ (35.0 / 3.0) * (r * other.value() * r ) / dist_9 * ( Matrix3(r, r) * q ).trace();
    + (35.0 / 3.0) * (r * other.value() * r) * (r * q * r)  / dist_9;
  }

  double energy(const Dipole& other, const Vector3& r, const double dist_5, const double dist_7) const
  {
    return -2.0 * other.value() * q * r / dist_5
        + 5.0 * ( other.value() * r ) * ( r * q * r ) / dist_7;
  }

  double energy(const Monopole& other, const Vector3& r, const double dist_5) const
  {
    return other.value() * ( r * q * r ) / dist_5;
  }

  double potential(const Vector3& r, const double dist_5) const
  {
    return energy(Monopole(1.0), r, dist_5);
  }

  Vector3 field(const Vector3& r, const double r_5, const double r_7) const
  {
    return ( -2.0 * ( q * r ) / r_5 ) + ( 5.0 * r * ( r * q * r ) / r_7 );
  }

  Vector3 field(const Vector3& r) const
  {
    double l2 = r.length_squared();
    double l5 = l2 * l2 * sqrt(l2);
    return field(r, l5, l5 * l2);
  }

  void rotate(const Matrix3& rot) { q = rot.transform(q); }

  bool operator==(const Quadrupole& other) const
  {
    return q == other.q;
  }

  bool operator!=(const Quadrupole& other) const
  {
    return q != other.q;
  }

  void printMe(std::ostream& strm = std::cout) const
  {
    strm << "Quadrupole: " << std::setprecision(5) << q;
  }

  std::string readMe(std::istream& strm, std::string header = "")
  {
    Matrix3 t;
    if (header.empty())
    {
      strm >> header;
    }
    if (header != "Quadrupole:") return header;
    strm >> t;
    q = t;
    return "";
  }

private:
  Matrix3 q;
};

#ifdef SWIG
__STR__(Quadrupole);
#endif // SWIG

inline double Monopole::energy(const Dipole& other, const Vector3& r, const double dist_3) const
{
  return other.energy(*this, -r, dist_3);
}

inline double Monopole::energy(const Quadrupole& other, const Vector3& r, const double dist_5) const
{
  return other.energy(*this, -r, dist_5);
}

inline double Dipole::energy(const Quadrupole& other, const Vector3& r, const double dist_5, const double dist_7) const
{
  return other.energy(*this, -r, dist_5, dist_7);
}


class Multipoles : public OpenBabel::OBGenericData
{
public:
  Multipoles(const Vector3& center = Vector3()) : cent(center)
  {
    _type = OpenBabel::OBGenericDataType::CustomData3;
    SetOrigin(OpenBabel::userInput);
    SetAttribute("Multipoles");
    mult[0] = 0;
    mult[1] = 0;
    mult[2] = 0;
  }

  Multipoles(const Multipoles& other) : cent(other.cent)
  {
    _type = OpenBabel::OBGenericDataType::CustomData3;
    SetOrigin(OpenBabel::userInput);
    SetAttribute("Multipoles");
    
    if (other.mult[0]) mult[0] = new Monopole(*( static_cast<Monopole *>(other.mult[0]) )); else mult[0] = 0;
    if (other.mult[1]) mult[1] = new Dipole(*( static_cast<Dipole *>(other.mult[1]) )); else mult[1] = 0;
    if (other.mult[2]) mult[2] = new Quadrupole(*( static_cast<Quadrupole *>(other.mult[2]) )); else mult[2] = 0;
  }

  Multipoles(const OpenBabel::OBGenericData* data)
  {
    _type = OpenBabel::OBGenericDataType::CustomData3;
    SetOrigin(OpenBabel::userInput);
    SetAttribute("Multipoles");
    const Multipoles* other = static_cast<const Multipoles*>(data);
    
    if (other)
    {
      if (other->mult[0]) mult[0] = new Monopole(*( static_cast<Monopole *>(other->mult[0]) )); else mult[0] = 0;
      if (other->mult[1]) mult[1] = new Dipole(*( static_cast<Dipole *>(other->mult[1]) )); else mult[1] = 0;
      if (other->mult[2]) mult[2] = new Quadrupole(*( static_cast<Quadrupole *>(other->mult[2]) )); else mult[2] = 0;
      cent = other->cent;
    } else {
      mult[0] = mult[1] = mult[2] = 0;
      cent = Vector3();
    }
      
  }

  ~Multipoles()
  {
    if (mult[0]) delete mult[0];
    if (mult[1]) delete mult[1];
    if (mult[2]) delete mult[2];
  }

  void setMonopole(const Monopole& mon)
  {
    if (mult[0]) delete mult[0];
    mult[0] = new Monopole(mon);
  }

  void setDipole(const Dipole& dip)
  {
    if (mult[1]) delete mult[1];
    mult[1] = new Dipole(dip);
  }

  void setQuadrupole(const Quadrupole& quad)
  {
    if (mult[2]) delete mult[2];
    mult[2] = new Quadrupole(quad);
  }

  bool hasMonopole() const { return mult[0]; }

  bool hasDipole() const { return mult[1]; }

  bool hasQuadrupole() const { return mult[2]; }

  void deleteMonopole() { delete mult[0]; mult[0] = 0; }

  void deleteDipole() { delete mult[1]; mult[1] = 0; }

  void deleteQuadrupole() { delete mult[2]; mult[2] = 0; }

  int order() const { if (hasQuadrupole()) return 2; else if (hasDipole()) return 1; else return 0; }

  const Monopole& monopole() const { if (mult[0]) return *( static_cast<Monopole *>(mult[0]) ); else throw std::logic_error("Request for monopole, but no monopole exists"); }

  const Dipole& dipole() const { if (mult[1]) return *( static_cast<Dipole *>(mult[1]) ); else throw std::logic_error("Request for dipole, but no dipole exists"); }

  const Quadrupole& quadrupole() const { if (mult[2]) return *( static_cast<Quadrupole *>(mult[2]) ); else throw std::logic_error("Request for quadrupole, but no quadrupole exists"); }

  Monopole& monopole() { if (mult[0]) return *( static_cast<Monopole *>(mult[0]) ); else throw std::logic_error("Request for monopole, but no monopole exists"); }

  Dipole& dipole() { if (mult[1]) return *( static_cast<Dipole *>(mult[1]) ); else throw std::logic_error("Request for dipole, but no dipole exists"); }

  Quadrupole& quadrupole() { if (mult[2]) return *( static_cast<Quadrupole *>(mult[2]) ); else throw std::logic_error("Request for monopole, but no monopole exists"); }

  Multipoles& operator=(const Multipoles& other)
  {
    if (this != &other)
    {
      if (other.hasMonopole()) setMonopole(other.monopole());
      if (other.hasDipole()) setDipole(other.dipole());
      if (other.hasQuadrupole()) setQuadrupole(other.quadrupole());
      cent = other.cent;
    }
    return *this;
  }

  const Multipole& operator[](const int i) const
  {
    return *( mult[i] );
  }

  Multipole& operator[](const int i)
  {
    return *( mult[i] );
  }

  Multipole *multipole(const MULTIPOLES i)
  {
    switch (i)
    {
      case MONOPOLE:
        return mult[0];
      case DIPOLE:
        return mult[1];
      case QUADRUPOLE:
        return mult[2];
    }
  }

  virtual OpenBabel::OBGenericData *Clone(OpenBabel::OBBase *) const { return new Multipoles(*this); }

  const Vector3& center() const { return cent; }

  void move(const Vector3& dr) { cent += dr; }

  void moveTo(const Vector3& r) { cent = r; }

  void rotate(const Matrix3& rot)
  {
    if (mult[1]) mult[1]->rotate(rot);
    if (mult[2]) mult[2]->rotate(rot);
  }

  double potential(const Vector3& r, const double dist, const double dist_2, const double dist_3, const double dist_5) const
  {
    double result = 0.0;
    if (hasMonopole()) result += monopole().potential(dist);
    if (hasDipole()) result += dipole().potential(r, dist_3);
    if (hasQuadrupole()) result += quadrupole().potential(r, dist_5);
    return result;
  }

  double potential(const Vector3& point) const
  {
    Vector3 r = ( point - cent ) / bohr;
    double dist_2 = r.length_squared();
    double dist = sqrt(dist_2);
    double dist_3 = dist * dist_2;
    double dist_5 = dist_3 * dist_2;
    return potential(-r, dist, dist_2, dist_3, dist_5);
  }

  double energy(const Monopole& m, const Vector3& r, const double dist, const double dist_2, const double dist_3, const double dist_5) const
  {
    return potential(r, dist, dist_2, dist_3, dist_5) * m.value();
  }

  double energy(const Monopole& m, const Vector3& point) const
  {
    Vector3 r = ( cent - point ) / bohr;
    double dist_2 = r.length_squared();
    double dist = sqrt(dist_2);
    double dist_3 = dist * dist_2;
    double dist_5 = dist_3 * dist_2;
    return energy(m, r, dist, dist_2, dist_3, dist_5);
  }

  double energy(const Dipole& d, const Vector3& r, const double dist_2, const double dist_3, const double dist_5, const double dist_7) const
  {
    double result = 0.0;
    if (hasMonopole()) result += monopole().energy(d, r, dist_3);
    if (hasDipole()) result += dipole().energy(d, r, dist_2, dist_5);
    if (hasQuadrupole()) result += quadrupole().energy(d, r, dist_5, dist_7);
    return result;
  }

  double energy(const Dipole& d, const Vector3& point) const
  {
    Vector3 r = ( cent - point ) / bohr;
    double dist_2 = r.length_squared();
    double dist = sqrt(dist_2);
    double dist_3 = dist * dist_2;
    double dist_5 = dist_3 * dist_2;
    double dist_7 = dist_5 * dist_2;
    return energy(d, r, dist_2, dist_3, dist_5, dist_7);
  }

  double energy(const Quadrupole& q, const Vector3& r, const double dist_5, const double dist_7, const double dist_9) const
  {
    double result = 0.0;
    if (hasMonopole()) result += monopole().energy(q, r, dist_5);
    if (hasDipole()) result += dipole().energy(q, r, dist_5, dist_7);
    if (hasQuadrupole()) result += quadrupole().energy(q, r, dist_5, dist_7, dist_9);
    return result;
  }

  double energy(const Quadrupole& q, const Vector3& point) const
  {
    Vector3 r = ( cent - point ) / bohr;
    double dist_2 = r.length_squared();
    double dist = sqrt(dist_2);
    double dist_3 = dist * dist_2;
    double dist_5 = dist_3 * dist_2;
    double dist_7 = dist_5 * dist_2;
    double dist_9 = dist_7 * dist_2;
    return energy(q, r, dist_5, dist_7, dist_9);
  }

  double energy(const Multipoles& other) const
  {
    Vector3 r = ( cent - other.center() ) / bohr;
    int order_ = order() + other.order();
    double dist, dist_2, dist_3, dist_5, dist_7, dist_9;
    switch (order_)
    {
      case 0:
        dist = r.length();
        break;
      case 1:
        dist_2 = r.length_squared();
        dist = sqrt(dist_2);
        dist_3 = dist_2 * dist;
        break;
      case 2:
        dist_2 = r.length_squared();
        dist = sqrt(dist_2);
        dist_3 = dist_2 * dist;
        dist_5 = dist_3 * dist_2;
        break;
      case 3:
        dist_2 = r.length_squared();
        dist = sqrt(dist_2);
        dist_3 = dist_2 * dist;
        dist_5 = dist_3 * dist_2;
        dist_7 = dist_5 * dist_2;
        break;
      case 4:
        dist_2 = r.length_squared();
        dist = sqrt(dist_2);
        dist_3 = dist_2 * dist;
        dist_5 = dist_3 * dist_2;
        dist_7 = dist_5 * dist_2;
        dist_9 = dist_7 * dist_2;
        break;
      default:
        break;
    }
    double result = 0.0;
    if (other.hasMonopole()) result += energy(other.monopole(), r, dist, dist_2, dist_3, dist_5);
    if (other.hasDipole()) result += energy(other.dipole(), r, dist_2, dist_3, dist_5, dist_7);
    if (other.hasQuadrupole()) result += energy(other.quadrupole(), r, dist_5, dist_7, dist_9);
    return result;
  }

  Vector3 field(const Vector3 r, const double dist_3, const double dist_5, const double dist_7) const
  {
    Vector3 result(0.0, 0.0, 0.0);
    if (hasMonopole()) result += monopole().field(r, dist_3);
    if (hasDipole()) result += dipole().field(r, dist_3, dist_5);
    if (hasQuadrupole()) result += quadrupole().field(r, dist_5, dist_7);
    return result;
  }

  Vector3 field(const Vector3& point) const
  {
    Vector3 r = ( point - cent ) / bohr;
    double dist_2 = r.length_squared();
    double dist = sqrt(dist_2);
    double dist_3 = dist * dist_2;
    double dist_5 = dist_3 * dist_2;
    return field(r, dist_3, dist_5, dist_5 * dist_2);
  }

  bool operator==(const Multipoles& other) const
  {
    return cent == other.cent
        && *static_cast<Monopole *>(mult[0]) == *static_cast<Monopole *>(other.mult[0])
        && *static_cast<Dipole *>(mult[1]) == *static_cast<Dipole *>(other.mult[1])
        && *static_cast<Quadrupole *>(mult[2]) == *static_cast<Quadrupole *>(other.mult[2]);
  }

  bool operator!=(const Multipoles& other) const
  {
    return !( *this == other );
  }

  void printMe(std::ostream& strm = std::cout) const
  {
    strm << "Multipoles: " << cent << std::endl;
    if (hasMonopole()) strm << "  " << monopole() << std::endl;
    if (hasDipole()) strm << "  " << dipole() << std::endl;
    if (hasQuadrupole()) strm << "  " << quadrupole() << std::endl;
  }

  std::string readMe(std::istream& strm, std::string header = "")
  {
    if (header.empty())
    {
      strm >> header;
    }
    if (header != "Multipoles:") return header;
    Vector3 cent_;
    strm >> cent_;
    cent = cent_;
    Monopole m;
    Dipole d;
    Quadrupole q;
    std::string buff;
    for (int i = 0; i < 3; ++i)
    {
      buff = "";
      strm >> buff;
      if (buff == "Monopole:")
      {
        m.readMe(strm, buff);
        setMonopole(m);
      } else if (buff == "Dipole:")
      {
        d.readMe(strm, buff);
        setDipole(d);
      } else if (buff == "Quadrupole:")
      {
        q.readMe(strm, buff);
        setQuadrupole(q);
      }
      else return buff;
    }
    strm >> buff;
    return buff;
  }

private:
  Vector3 cent;
  Multipole *mult[3];
};

#ifndef SWIG
DECLARE_IO_STREAM_OPERATORS(Multipoles)
#endif // SWIG

#ifdef SWIG
__STR__(Multipoles);
#endif // SWIG

#endif
