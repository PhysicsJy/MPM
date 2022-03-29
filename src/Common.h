//////////////////////////////////////////////////////////////////////////
////Dartmouth Physical Computing Starter Code
////http://www.dartmouth.edu/~boolzhu/cosc89.18.html
//////////////////////////////////////////////////////////////////////////

#ifndef __Common_h__
#define __Common_h__
#include <vector>
#include <list>
#include <queue>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <iostream>
#include "Eigen/Dense"
#include "Eigen/Sparse"

#define USE_FLOAT true
//////////////////////////////////////////////////////////////////////////
//// Eigen vector and matrix types

//// To define a vector or matrix type, you may use the following aliased type names:
//// Vector<double,2>, Vector<float,3>, Vector<int,4>, Matrix<float,3>, etc.
//// or Vector2d, Vector3f, Vector4i, Matrix3f, etc.

template<class T,int d> using Vector=Eigen::Matrix<T,d,1>;
template<class T,int d> using Matrix=Eigen::Matrix<T,d,d>;

#define Declare_Eigen_Vector_Types(type,t)		\
using Vector1##t=Eigen::Matrix<type,1,1>;       \
using Vector2##t=Eigen::Vector2##t;             \
using Vector3##t=Eigen::Vector3##t;             \
using Vector4##t=Eigen::Vector4##t;             \
using VectorX##t=Eigen::VectorX##t;

#define Declare_Eigen_Matrix_Types(type,t)		\
using Matrix1##t=Eigen::Matrix<type,1,1>;       \
using Matrix2##t=Eigen::Matrix2##t;             \
using Matrix3##t=Eigen::Matrix3##t;             \
using Matrix4##t=Eigen::Matrix4##t;             \
using MatrixX##t=Eigen::MatrixX##t;  

Declare_Eigen_Vector_Types(int,i)
Declare_Eigen_Vector_Types(float,f)
Declare_Eigen_Vector_Types(double,d)
Declare_Eigen_Matrix_Types(int,i)
Declare_Eigen_Matrix_Types(float,f)
Declare_Eigen_Matrix_Types(double,d)

#define Declare_Eigen_Types(type,t)     \
using real=type;						\
using Vector1=Eigen::Matrix<real,1,1>;	\
using Vector2=Eigen::Vector2##t;        \
using Vector3=Eigen::Vector3##t;        \
using Vector4=Eigen::Vector4##t;        \
using VectorX=Eigen::VectorX##t;        \
using Matrix2=Eigen::Matrix2##t;        \
using Matrix3=Eigen::Matrix3##t;        \
using Matrix4=Eigen::Matrix4##t;        \
using MatrixX=Eigen::MatrixX##t;        \

#ifdef USE_FLOAT
Declare_Eigen_Types(float,f)
#else
Declare_Eigen_Types(double,d)
#endif

template<class T> using Array=std::vector<T>;
template<class T> using ArrayPtr=std::shared_ptr<Array<T> >;
using size_type=Array<int>::size_type;

#endif
