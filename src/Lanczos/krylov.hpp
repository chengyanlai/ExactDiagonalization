#ifndef __KRYLOV_HPP__
#define __KRYLOV_HPP__
#include "src/EDType.hpp"
#include "src/ArmadilloMatrix.hpp"

void krylovEXP(const ComplexSparseMatrixType &A, ComplexVectorType &Vec, const ComplexType Prefactor, const size_t Kmax = 10, const double threshNorm = 1.0e-12);
void krylov(const ComplexSparseMatrixType &A, const ComplexVectorType &Vec, ComplexMatrixType& OTrans, RealVectorType& Dvec, const size_t Kmax = 10, const double threshNorm = 1.0e-12);

#endif /* end of include guard: __KRYLOV_HPP__ */
