#include "stdafx.h"
#include "matrix3d.h"
#include "vector3d.h"

Matrix3d::Matrix3d(){
  int i,j;
  matrix=new double* [3];
  for(i=0;i<3;++i){
    matrix[i]=new double [3];
	for (j=0;j<3;++j){
      if (i==j){
	    matrix[i][j]=1;
	  }else{
	    matrix[i][j]=0;
	  }
	}
  }
}

Matrix3d::  Matrix3d(double m11,double m12,double m13,double m21,double m22,double m23,double m31,double m32,double m33){
	matrix=new double* [3];
	for(int i=0;i<3;++i){
		matrix[i]=new double [3];
	}
	matrix[0][0] = m11;
	matrix[0][1] = m12;
	matrix[0][2] = m13;
	matrix[1][0] = m21;
	matrix[1][1] = m22;
	matrix[1][2] = m23;
	matrix[2][0] = m31;
	matrix[2][1] = m32;
	matrix[2][2] = m33;
}

Matrix3d::~Matrix3d(){
  int i;
  for(i=0;i<3;++i){
	delete[] matrix[i];
  }
  delete[] matrix;
}
Matrix3d* Matrix3d::Inverse(){
  Matrix3d* inv = new Matrix3d();
  double det;
  det=matrix[0][0]*(matrix[1][1]*matrix[2][2]-matrix[1][2]*matrix[2][1])-matrix[0][1]*(matrix[1][0]*matrix[2][2]-matrix[1][2]*matrix[2][0])+matrix[0][2]*(matrix[1][0]*matrix[2][1]-matrix[1][1]*matrix[2][0]);
  if (det!=0){
	inv->matrix[0][0]=(1/det)*(matrix[1][1]*matrix[2][2]-matrix[1][2]*matrix[2][1]);
	inv->matrix[0][1]=-(1/det)*(matrix[0][1]*matrix[2][2]-matrix[0][2]*matrix[2][1]);
	inv->matrix[0][2]=(1/det)*(matrix[0][1]*matrix[1][2]-matrix[0][2]*matrix[1][1]);
	inv->matrix[1][0]=-(1/det)*(matrix[1][0]*matrix[2][2]-matrix[1][2]*matrix[2][0]);
	inv->matrix[1][1]=(1/det)*(matrix[0][0]*matrix[2][2]-matrix[0][2]*matrix[2][0]);
	inv->matrix[1][2]=-(1/det)*(matrix[0][0]*matrix[1][2]-matrix[1][0]*matrix[0][2]);
	inv->matrix[2][0]=(1/det)*(matrix[1][0]*matrix[2][1]-matrix[1][1]*matrix[2][0]);
	inv->matrix[2][1]=-(1/det)*(matrix[0][0]*matrix[2][1]-matrix[0][1]*matrix[2][0]);
	inv->matrix[2][2]=(1/det)*(matrix[0][0]*matrix[1][1]-matrix[1][0]*matrix[0][1]);
  }
  return inv;
}
double Matrix3d::determinant(){
  double temp;
  temp=matrix[0][0]*(matrix[1][1]*matrix[2][2]-matrix[1][2]*matrix[2][1])-matrix[0][1]*(matrix[1][0]*matrix[2][2]-matrix[1][2]*matrix[2][0])+matrix[0][2]*(matrix[0][1]*matrix[2][1]-matrix[1][1]*matrix[0][2]);
  return temp;
}
double Matrix3d::getvalue(int i,int j){
  return matrix[i][j];
}
void Matrix3d::setvalue(int i,int j, double value){
  matrix[i][j]=value;
}
Matrix3d& Matrix3d::operator=(const Matrix3d& matr){
  
  int i;
  int j;
  for(i=0;i<3;++i){
	for(j=0;j<3;++j){
		matrix[i][j]=matr.matrix[i][j];
	}
  }
  return *this;
}

vector3d Matrix3d::multiply(const vector3d& v){
	double x = matrix[0][0]*v.x + matrix[0][1]*v.y + matrix[0][2]*v.z;
	double y = matrix[1][0]*v.x + matrix[1][1]*v.y + matrix[1][2]*v.z;
	double z = matrix[2][0]*v.x + matrix[2][1]*v.y + matrix[2][2]*v.z;
	return vector3d(x,y,z);
}

Matrix3d* Matrix3d::createBasisByOneVector(const vector3d& v){
	double theta = acos(v.z/sqrt(v.z*v.z + v.y*v.y + v.x*v.x));
	double phi = atan2(v.y,v.x);
	return new Matrix3d(sin(phi),cos(theta)*cos(phi),sin(theta)*cos(phi),-cos(phi),cos(theta)*sin(phi),sin(theta)*sin(phi),0,-sin(theta),cos(theta));
}