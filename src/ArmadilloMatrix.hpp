#ifndef __ARMODILLO_MATRIX_HPP__
#define __ARMODILLO_MATRIX_HPP__

#include <armadillo>
#include "EDType.hpp"

/** Dense complex matrix. */
typedef arma::Mat<ComplexType> ComplexMatrixType;
/** Dense real matrix. */
typedef arma::Mat<RealType> RealMatrixType;
/** Dense unsigned long matrix. */
typedef arma::umat ULongMatrixType;

/** Dense complex vector. */
typedef arma::Col<ComplexType> ComplexVectorType;
/** Dense real vector. */
typedef arma::Col<RealType> RealVectorType;
/** Dense vector of integers. */
typedef arma::Col<int> IntVectorType;

/** Sparse complex matrix. */
typedef arma::SpMat<ComplexType> ComplexSparseMatrixType;
/** Sparse real matrix. */
typedef arma::SpMat<RealType> RealSparseMatrixType;

// inline RealSparseMatrixType BuildSparseMatrix( sizt_t& dim, ULongMatrixType& Locations, RealVectorType& Values ){
//   RealSparseMatrixType out(true, Locations, Values, dim, dim);//, sort_locations = true, check_for_zeros = true);
//   return out;
// };
// inline ComplexSparseMatrixType BuildSparseMatrix( sizt_t& dim, ULongMatrixType& Locations, ComplexVectorType& Values ){
//   ComplexSparseMatrixType out(true, Locations, Values, dim, dim);//, sort_locations = true, check_for_zeros = true);
//   return out;
// };
#endif	/* end of include guard: __ARMODILLO_MATRIX_HPP__ */
