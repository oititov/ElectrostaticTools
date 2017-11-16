/*************************************************************************
 * vector3.h - Hanles 3-dimensional vectors.                             *
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

#ifndef VECTOR3_HEADER
#define VECTOR3_HEADER

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <stdexcept>

class Vector3
{
public:
  typedef double float_type;

  Vector3(float_type a, float_type b, float_type c)
  {
    x_[0] = a;
    x_[1] = b;
    x_[2] = c;
  }

  Vector3()
  {
    x_[0] = 0;
    x_[1] = 0;
    x_[2] = 0;
  }

  Vector3(const float_type *v)
  {
    x_[0] = v[0];
    x_[1] = v[1];
    x_[2] = v[2];
  }

  Vector3(const Vector3& other)
  {
    x_[0] = other.x_[0];
    x_[1] = other.x_[1];
    x_[2] = other.x_[2];
  }

  float_type length() const { return std::sqrt(x_[0] * x_[0] + x_[1] * x_[1] + x_[2] * x_[2]); }

  float_type length_squared() const { return ( x_[0] * x_[0] + x_[1] * x_[1] + x_[2] * x_[2] ); }

  Vector3& normalize()
  {
    float_type n = length();
    n = float_type(1.0) / ( n ? n : 1.0 );
    this->x_[0] *= n;
    this->x_[1] *= n;
    this->x_[2] *= n;
    return *this;
  }

  Vector3 normalized() const
  {
    Vector3 v(*this);
    v.normalize();
    return v;
  }

  const float_type& x() const { return x_[0]; }

  const float_type& y() const { return x_[1]; }

  const float_type& z() const { return x_[2]; }

#ifndef SWIG   
  float_type& x() { return x_[0]; }
  
  float_type& y() { return x_[1]; }
  
  float_type& z() { return x_[2]; }
#endif // SWIG
  
  void x(const double v) { x_[0] = v; }

  void y(const double v) { x_[1] = v; }

  void z(const double v) { x_[2] = v; }
  
  const float_type& operator[](unsigned int i) const { return x_[i]; }

  float_type& operator[](unsigned int i) { return x_[i]; }

  const float_type& operator()(unsigned int i) const { return ( *this )[i]; }

#ifndef SWIG  
  float_type& operator()(unsigned int i) { return ( *this )[i]; }
#endif // SWIG
  
  const float_type *vertex() const { return x_; }

  float_type dot(const Vector3& v) const { return x_[0] * v.x_[0] + x_[1] * v.x_[1] + x_[2] * v.x_[2]; }

  float_type operator*(const Vector3& v) const { return v.dot(*this); }

  Vector3 vec(const Vector3& v) const
  {
    return Vector3(
        x_[1] * v.x_[2] - v.x_[1] * x_[2], x_[2] * v.x_[0] - v.x_[2] * x_[0], x_[0] * v.x_[1] - v.x_[0] * x_[1]
                  );
  }

  bool operator==(const Vector3& other) const
  {
    return x_[0] == other.x_[0] && x_[1] == other.x_[1] && x_[2] == other.x_[2];
  }

  bool operator!=(const Vector3& other) const
  {
    return x_[0] != other.x_[0] || x_[1] != other.x_[1] || x_[2] != other.x_[2];
  }

  Vector3 operator+(const Vector3& other) const
  {
    return Vector3(x_[0] + other.x_[0], x_[1] + other.x_[1], x_[2] + other.x_[2]);
  }

  Vector3 operator-() const
  {
    return Vector3(-x_[0], -x_[1], -x_[2]);
  }

  Vector3 operator-(const Vector3& other) const
  {
    return *this + ( -other );
  }

  Vector3& operator+=(const Vector3& other)
  {
    return *this = *this + other;
  }

  Vector3& operator-=(const Vector3& other)
  {
    return *this = *this - other;
  }

  Vector3 operator*(const float_type& a) const
  {
    return Vector3(a * x_[0], a * x_[1], a * x_[2]);
  }

  Vector3 operator/(const float_type& a) const
  {
    return Vector3(x_[0] / a, x_[1] / a, x_[2] / a);
  }

  Vector3& operator*=(const float_type& other)
  {
    return *this = *this * other;
  }

  Vector3& operator/=(const float_type& other)
  {
    return *this = *this / other;
  }
  
  void printMe(std::ostream& strm = std::cout) const
  {
    strm << "( " << std::setprecision(5) << x_[0] << ", "
        << std::setprecision(5) << x_[1] << ", "
        << std::setprecision(5) << x_[2] << " )";
  }

  std::string readMe(std::istream& strm)
  {
    char c;
    strm >> c;
    if (c != '(') throw std::domain_error("Reading vector: \"(\" not found");
    strm >> x_[0] >> c;
    if (c != ',') throw std::domain_error("Reading vector: \",\" not found");
    strm >> x_[1] >> c;
    if (c != ',') throw std::domain_error("Reading vector: \",\" not found");
    strm >> x_[2] >> c;
    if (c != ')') throw std::domain_error("Reading vector: \")\" not found");
    return std::string();
  }

private:
  float_type x_[3];
};

#ifdef SWIG
%rename(__mult__) operator*(const double& a, const Vector3& vec);
#endif // SWIG

inline Vector3 operator*(const double& a, const Vector3& vec)
{
  return vec * a;
}

#ifndef SWIG
inline std::ostream& operator<<(std::ostream& strm, const Vector3& vec)
{
  vec.printMe(strm);
  return strm;
}

inline std::istream& operator>>(std::istream& strm, Vector3& vec)
{
  vec.readMe(strm);
  return strm;
}
#endif // SWIG

#ifdef SWIG
__STR__(Vector3);
#endif // SWIG

#endif // VECTOR3_HEADER
