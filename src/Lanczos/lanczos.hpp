#ifndef __LANCZOS_H__
#define __LANCZOS_H__
#include "src/EDType.hpp"
#include "src/EigenMatrix.hpp"

bool LanczosEV(const size_t N, const RealSparseMatrixType A,
  RealVectorType &Vec, RealType &Val,
  size_t &max_iter, double err_tol = 1.0E-10);

bool LanczosEV(const size_t N, const ComplexSparseMatrixType A,
  ComplexVectorType &Vec, RealType &Val,
  size_t &max_iter, double err_tol = 1.0E-10);

#endif//__LANCZOS_H__
