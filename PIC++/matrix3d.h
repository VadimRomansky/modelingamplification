#ifndef MATRIX3D_H
#define MATRIX3D_H

#include "vector3d.h"

class Matrix3d{
public:
  double** matrix;
  Matrix3d();
  Matrix3d(double m11,double m12,double m13,double m21,double m22,double m23,double m31,double m32,double m33);
  ~Matrix3d();
  Matrix3d* Inverse();
  double determinant();
  double getvalue(int i,int j);
  void setvalue(int i,int j, double value);
  Matrix3d& operator=(const Matrix3d& matr);
  vector3d operator*(const vector3d& v);
  static Matrix3d* createBasisByOneVector(const vector3d& v);
};



#endif