/*************************************************************************
 * matrix4.h - Hanles 4x4 matrices.                                      *
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

#ifndef MATRIX4_HEADER
#define MATRIX4_HEADER

#include <stdexcept>

#include "vector3.h"

class Matrix4
{
public:
  typedef double float_type;

  Matrix4() { }

  Matrix4(float_type *data)
  {
    for (int i = 0; i < 16; ++i) data_[i] = data[i];
  }

  Matrix4(float_type **data)
  {
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) data_[i * 4 + j] = data[i][j];
  }

  Matrix4(const Vector3& f, const Vector3& s)
  {
    data_[0] = f[0] * s[0];
    data_[1] = f[0] * s[1];
    data_[2] = f[0] * s[2];
    data_[3] = 0.0;
    data_[4] = f[1] * s[0];
    data_[5] = f[1] * s[1];
    data_[6] = f[1] * s[2];
    data_[7] = 0.0;
    data_[8] = f[2] * s[0];
    data_[9] = f[2] * s[1];
    data_[10] = f[2] * s[2];
    data_[11] = 0.0;
    data_[12] = 0.0;
    data_[13] = 0.0;
    data_[14] = 0.0;
    data_[15] = 0.0;
  }

  void loadIdentity() { for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) data_[i * 4 + j] = ( i == j ) ? 1 : 0; }

  void loadZero() { for (int i = 0; i < 16; ++i) data_[i] = 0.0; }

  const float_type *operator[](int i) const { return &( data_[i * 4] ); }

  float_type *operator[](int i) { return &( data_[i * 4] ); }
  
  const float_type& operator()(unsigned int i, unsigned int j) const { return ( *this )[i][j]; }

#ifndef SWIG  
  float_type& operator()(unsigned int i, unsigned int j) { return ( *this )[i][j]; }
#endif // SWIG
  
  Vector3 operator*(const Vector3& vec) const
  {
    return Vector3(data_[0] * vec[0] + data_[1] * vec[1] + data_[2] * vec[2] + data_[3],
                   data_[4] * vec[0] + data_[5] * vec[1] + data_[6] * vec[2] + data_[7],
                   data_[8] * vec[0] + data_[9] * vec[1] + data_[10] * vec[2] + data_[11]);
  }

  Matrix4 operator*(const Matrix4& mat) const
  {
    Matrix4 a;
    a.loadZero();
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        for (int k = 0; k < 4; ++k)
        {
          a[i][j] += data_[i * 4 + k] * mat.data_[k * 4 + j];
        }
    return a;
  }

  const Matrix4& operator*=(const Matrix4& mat)
  {
    return *this = mat * *this;
  }

  // some convenience functions follow
  void rotate(const Vector3& axis = Vector3(1.0, 0.0, 0.0), const float_type angle = 0.0)
  {
    Matrix4 rot;
    Vector3 a = axis;
    a.normalize();
    Vector3 b = Vector3(a[0] * a[0], a[1] * a[1], a[2] * a[2]);
    float_type c = cos(angle), s = sin(angle);
    rot[0][0] = b[0] + ( 1.0 - b[0] ) * c;
    rot[0][1] = a[0] * a[1] * ( 1.0 - c ) - a[2] * s;
    rot[0][2] = a[0] * a[2] * ( 1.0 - c ) + a[1] * s;
    rot[0][3] = 0.0;
    rot[1][0] = a[0] * a[1] * ( 1.0 - c ) + a[2] * s;
    rot[1][1] = b[1] + ( 1.0 - b[1] ) * c;
    rot[1][2] = a[1] * a[2] * ( 1.0 - c ) - a[0] * s;
    rot[1][3] = 0.0;
    rot[2][0] = a[0] * a[2] * ( 1.0 - c ) - a[1] * s;
    rot[2][1] = a[1] * a[2] * ( 1.0 - c ) + a[0] * s;
    rot[2][2] = b[2] + ( 1.0 - b[2] ) * c;
    rot[2][3] = 0.0;
    rot[3][0] = 0.0;
    rot[3][1] = 0.0;
    rot[3][2] = 0.0;
    rot[3][3] = 1.0;
    *this *= rot;
  }

  void scale(const Vector3& vec)
  {
    Matrix4 sc;
    sc[0][0] = vec[0];
    sc[0][1] = 0;
    sc[0][2] = 0;
    sc[0][3] = 0;
    sc[1][0] = 0;
    sc[1][1] = vec[1];
    sc[1][2] = 0;
    sc[1][3] = 0;
    sc[2][0] = 0;
    sc[2][1] = 0;
    sc[2][2] = vec[2];
    sc[2][3] = 0;
    sc[3][0] = 0;
    sc[3][1] = 0;
    sc[3][2] = 0;
    sc[3][3] = 1;
    *this *= sc;
  }

  void translate(const Vector3& vec)
  {
    Matrix4 mat;
    mat.loadIdentity();
    mat[0][3] = vec[0];
    mat[1][3] = vec[1];
    mat[2][3] = vec[2];
    *this *= mat;
  }

  void printMe(std::ostream& strm = std::cout) const
  {
    strm << "(";
    for (int i = 0; i < 15; ++i) strm << " " << std::setprecision(5) << data_[i];
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
    for (int i = 1; i < 16; ++i)
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
  float_type data_[16];
};

#ifdef SWIG
%rename(__mult__) operator*(const Vector3& vec, const Matrix4& mat);
#endif // SWIG

inline Vector3 operator*(const Vector3& vec, const Matrix4& mat)
{
  return Vector3(mat[0][0] * vec[0] + mat[0][1] * vec[1] + mat[0][2] * vec[2] + mat[0][3],
                 mat[1][0] * vec[0] + mat[1][1] * vec[1] + mat[1][2] * vec[2] + mat[1][3],
                 mat[2][0] * vec[0] + mat[2][1] * vec[1] + mat[2][2] * vec[2] + mat[2][3]);
}

#ifndef SWIG
inline std::ostream& operator<<(std::ostream& strm, const Matrix4& mat)
{
  mat.printMe(strm);
  return strm;
}

inline std::istream& operator>>(std::istream& strm, Matrix4& mat)
{
  mat.readMe(strm);
  return strm;
}
#endif // SWIG

#ifdef SWIG
__STR__(Matrix4);
#endif // SWIG

#endif // MATRIX4_HEADER
