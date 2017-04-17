#ifndef __EIGEN_MATRIX_HPP__
#define __EIGEN_MATRIX_HPP__

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "EDType.hpp"

/** Dense complex matrix. */
typedef Eigen::Matrix<ComplexType, Eigen::Dynamic, Eigen::Dynamic, Eigen::AutoAlign|Eigen::RowMajor> ComplexMatrixType;
/** Dense real matrix. */
typedef Eigen::Matrix<RealType, Eigen::Dynamic, Eigen::Dynamic, Eigen::AutoAlign|Eigen::RowMajor> RealMatrixType;

/** Dense complex vector. */
typedef Eigen::Matrix<ComplexType, Eigen::Dynamic, 1, Eigen::AutoAlign> ComplexVectorType;
/** Dense real vector. */
typedef Eigen::Matrix<RealType, Eigen::Dynamic, 1, Eigen::AutoAlign> RealVectorType;
/** Dense vector of integers. */
typedef Eigen::Matrix<int, Eigen::Dynamic, 1, Eigen::AutoAlign> IntVectorType;

/** Solvers */
// typedef Eigen::Map<MatrixType> MapMatrix;
// typedef Eigen::SelfAdjointEigenSolver<MatrixType> Dia;
// typedef Eigen::JacobiSVD<MatrixType> SVD;

/** Map from C/C++ array*/
typedef Eigen::Map<RealVectorType> dMapVector;
typedef Eigen::Map<ComplexVectorType> zMapVector;
typedef Eigen::Map<RealMatrixType, Eigen::RowMajor> dMapMatrix;
typedef Eigen::Map<ComplexMatrixType, Eigen::RowMajor> zMapMatrix;

/** Sparse complex matrix. */
typedef Eigen::SparseMatrix<ComplexType, Eigen::AutoAlign|Eigen::RowMajor> ComplexSparseMatrixType;
/** Sparse real matrix. */
typedef Eigen::SparseMatrix<RealType, Eigen::AutoAlign|Eigen::RowMajor> RealSparseMatrixType;

/** Use to fill sparse matrix*/
typedef Eigen::Triplet<RealType> RealTriplet;
typedef Eigen::Triplet<ComplexType> ComplexTriplet;
#endif
