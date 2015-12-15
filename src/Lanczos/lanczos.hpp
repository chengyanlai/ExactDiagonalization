#ifndef __LANCZOS_HPP__
#define __LANCZOS_HPP__
#include "src/EDType.hpp"
#include "src/EigenMatrix.hpp"

bool LanczosEV(const RealSparseMatrixType A,
  RealVectorType &Vec, std::vector<RealType> &Val,
  size_t &max_iter, double err_tol = 1.0E-10);

bool LanczosEV(const ComplexSparseMatrixType A,
  ComplexVectorType &Vec, std::vector<RealType> &Val,
  size_t &max_iter, double err_tol = 1.0E-10);

#endif/* end of include guard: __LANCZOS_HPP__ */
