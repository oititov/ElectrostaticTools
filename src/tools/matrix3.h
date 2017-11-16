/*************************************************************************
 * matrix3.h - Hanles 3x3 matrices.                                      *
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

#ifndef MATRIX3_HEADER
#define MATRIX3_HEADER

#include <cstring>
#include <iostream>
#include <iomanip>
#include <stdexcept>

#include "vector3.h"

class Matrix3
{
public:
  typedef double float_type;

  Matrix3() { for (int i = 0; i < 9; ++i) data_[i] = 0.0; }

  Matrix3(const float_type *data)
  {
    for (int i = 0; i < 9; ++i) data_[i] = data[i];
  }

  Matrix3(const float_type **data)
  {
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) data_[i * 3 + j] = data[i][j];
  }

  Matrix3(const float_type data[3][3])
  {
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) data_[i * 3 + j] = data[i][j];
  }

  Matrix3(const Vector3& f, const Vector3& s)
  {
    data_[0] = f[0] * s[0];
    data_[1] = f[0] * s[1];
    data_[2] = f[0] * s[2];
    data_[3] = f[1] * s[0];
    data_[4] = f[1] * s[1];
    data_[5] = f[1] * s[2];
    data_[6] = f[2] * s[0];
    data_[7] = f[2] * s[1];
    data_[8] = f[2] * s[2];
  }

  Matrix3(const Vector3& x, const Vector3& y, const Vector3& z)
  {
    data_[0] = x[0];
    data_[1] = y[0];
    data_[2] = z[0];
    data_[3] = x[1];
    data_[4] = y[1];
    data_[5] = z[1];
    data_[6] = x[2];
    data_[7] = y[2];
    data_[8] = z[2];
  }

  Matrix3(const double v1, const double v2, const double v3,
          const double v4, const double v5, const double v6,
          const double v7, const double v8, const double v9)
  {
    data_[0] = v1;
    data_[1] = v2;
    data_[2] = v3;
    data_[3] = v4;
    data_[4] = v5;
    data_[5] = v6;
    data_[6] = v7;
    data_[7] = v8;
    data_[8] = v9;
  }

  Matrix3(const Matrix3& other)
  {
    memcpy(data_, other.data_, sizeof(float_type[9]));
  }

  Matrix3& operator=(const Matrix3& other)
  {
    memcpy(data_, other.data_, sizeof(float_type[9]));
  }

  void loadIdentity(const float_type e = 1.0) { for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) data_[i * 3 + j] = ( i == j ) ? e : 0; }

  void loadZero() { for (int i = 0; i < 9; ++i) data_[i] = 0.0; }

  const float_type *operator[](unsigned int i) const { return &( data_[i * 3] ); }

  float_type *operator[](unsigned int i) { return &( data_[i * 3] ); }

  const float_type& operator()(unsigned int i, unsigned int j) const { return data_[i * 3 + j]; }

#ifndef SWIG  
  float_type& operator()(unsigned int i, unsigned int j) { return data_[i * 3 + j]; }
#endif // SWIG
  
  const float_type *data() { return data_; }

  Vector3 operator*(const Vector3& vec) const
  {
    return Vector3(data_[0] * vec[0] + data_[1] * vec[1] + data_[2] * vec[2],
                   data_[3] * vec[0] + data_[4] * vec[1] + data_[5] * vec[2],
                   data_[6] * vec[0] + data_[7] * vec[1] + data_[8] * vec[2]);
  }

  Matrix3 operator*(const Matrix3& mat) const
  {
    Matrix3 a;
    a.loadZero();
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k)
        {
          a[i][j] += data_[i * 3 + k] * mat.data_[k * 3 + j];
        }
    return a;
  }

  const Matrix3& operator*=(const Matrix3& mat)
  {
    return *this = mat * *this;
  }
 
  const Matrix3& operator *= (const float_type f)
  {
    for (int i = 0; i < 9; ++i) data_[i] *= f;
    return *this;
  }
  
  Matrix3 operator * (const float_type f) const
  {
    Matrix3 result(*this);
    result *= f;
    return result;
  }
 
  const Matrix3 operator + (const Matrix3& other) const
  {
    return Matrix3(data_[0]+other.data_[0], data_[1]+other.data_[1], data_[2]+other.data_[2],
                   data_[3]+other.data_[3], data_[4]+other.data_[4], data_[5]+other.data_[5],
                   data_[6]+other.data_[6], data_[7]+other.data_[7], data_[8]+other.data_[8]);
  }

  const Matrix3& operator += (const Matrix3& other)
  {
    return *this = other + *this;
  }
 
  double trace() const { return data_[0] + data_[4] + data_[8]; }

  Matrix3 transposed() const
  {
    Matrix3 m;
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j) m.data_[i * 3 + j] = data_[j * 3 + i];
    }
    return m;
  }

  Matrix3 transform(const Matrix3& m) const
  {
    // return *this * m * this->transposed();
    Matrix3 t = *this * m;
    Matrix3 r;
    r.loadZero();
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k)
        {
          r[i][j] += t[i][k] * data_[j * 3 + k];
        }
    return r;
  }

  // some convenience functions follow
  void rotate(const Vector3& axis = Vector3(1.0, 0.0, 0.0), const float_type angle = 0.0)
  {
    Matrix3 rot;
    Vector3 a = axis;
    a.normalize();
    Vector3 b = Vector3(a[0] * a[0], a[1] * a[1], a[2] * a[2]);
    float_type c = cos(angle), s = sin(angle);
    rot[0][0] = b[0] + ( 1.0 - b[0] ) * c;
    rot[0][1] = a[0] * a[1] * ( 1.0 - c ) - a[2] * s;
    rot[0][2] = a[0] * a[2] * ( 1.0 - c ) + a[1] * s;
    rot[1][0] = a[0] * a[1] * ( 1.0 - c ) + a[2] * s;
    rot[1][1] = b[1] + ( 1.0 - b[1] ) * c;
    rot[1][2] = a[1] * a[2] * ( 1.0 - c ) - a[0] * s;
    rot[2][0] = a[0] * a[2] * ( 1.0 - c ) - a[1] * s;
    rot[2][1] = a[1] * a[2] * ( 1.0 - c ) + a[0] * s;
    rot[2][2] = b[2] + ( 1.0 - b[2] ) * c;
    *this *= rot;
  }

  bool operator==(const Matrix3& other) const
  {
    for (int i = 0; i < 9; ++i) if (data_[i] != other.data_[i]) { return false; }
  }

  bool operator!=(const Matrix3& other) const
  {
    return !( *this == other );
  }

  void printMe(std::ostream& strm = std::cout) const
  {
    strm << "( " << std::setprecision(5) << data_[0];
    for (int i = 1; i < 9; ++i) strm << ", " << std::setprecision(5) << data_[i];
    strm << " )";
  }

  void readMe(std::istream& strm)
  {
    char c;
    float_type a;
    strm >> c;
    if (c != '(') throw std::domain_error("Reading matrix: \"(\" not found");
    strm >> a;
    data_[0] = a;
    for (int i = 1; i < 9; ++i)
    {
      strm >> c;
      if (c != ',') throw std::domain_error("Reading matrix: \",\" not found");
      strm >> a;
      data_[i] = a;
    }
    strm >> c;
    if (c != ')') throw std::domain_error("Reading matrix: \")\" not found");
  }

private:
  float_type data_[9];
};

#ifdef SWIG
%rename(__mult__) operator*(const Vector3& vec, const Matrix3& mat);
#endif // SWIG

inline Vector3 operator*(const Vector3& vec, const Matrix3& mat)
{
  return Vector3(mat[0][0] * vec[0] + mat[1][0] * vec[1] + mat[2][0] * vec[2],
                 mat[0][1] * vec[0] + mat[1][1] * vec[1] + mat[2][1] * vec[2],
                 mat[0][2] * vec[0] + mat[1][2] * vec[1] + mat[2][2] * vec[2]);
}

#ifndef SWIG
inline std::ostream& operator<<(std::ostream& strm, const Matrix3& mat)
{
  mat.printMe(strm);
  return strm;
}

inline std::istream& operator>>(std::istream& strm, Matrix3& mat)
{
  mat.readMe(strm);
  return strm;
}
#endif // SWIG

#ifdef SWIG
__STR__(Matrix3);
#endif // SWIG

#endif // MATRIX3_HEADER
