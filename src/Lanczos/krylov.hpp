#ifndef __KRYLOV_HPP__
#define __KRYLOV_HPP__
#include "src/EDType.hpp"
#include "src/EigenMatrix.hpp"

void krylov(const ComplexSparseMatrixType &A, ComplexVectorType &Vec,
  const ComplexType Prefactor, const size_t Kmax = 10);

#endif /* end of include guard: __KRYLOV_HPP__ */
