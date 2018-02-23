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

template<typename Tnum>
inline arma::SpMat<Tnum> BuildSparseHamiltonian( const size_t& dim, const std::vector<std::tuple<int, int, Tnum> >& MatElemts ){
  ULongMatrixType Locations(2, MatElemts.size());
  arma::Col<Tnum> Values( MatElemts.size() );
  typename std::vector<std::tuple<int, int, Tnum> >::const_iterator it = MatElemts.begin();
  size_t cnt = 0;
  for (; it != MatElemts.end(); ++it ){
    int row, col;
    Tnum val;
    std::tie(row, col, val) = *it;
    Locations(0,cnt) = row;
    Locations(1,cnt) = col;
    Values(cnt) = val;
    cnt++;
  }
  // First true allows repeated matrix elements
  arma::SpMat<Tnum> out(true, Locations, Values, dim, dim);//, sort_locations = true, check_for_zeros = true);
  return out;
};

#endif	/* end of include guard: __ARMODILLO_MATRIX_HPP__ */
