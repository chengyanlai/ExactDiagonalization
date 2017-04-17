#ifndef LAPACK_H
#define LAPACK_H
#include <cstdint>
#include <limits.h>
#include <assert.h>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <complex>

void matrixMul(double* A, double* B, int M, int N, int K, double* C);
void vectorAdd(double* Y, double* X, size_t N);// Y = Y + X
void vectorScal(double a, double* X, size_t N);	// X = a * X
double vectorSum(double* X, size_t N, int inc);
double vectorNorm(double* X, size_t N, int inc);
void vectorExp(double a, double* X, size_t N);
void diagMM(double* diag, double* mat, size_t M, size_t N);

void syDiag(double* Kij, int N, double* Eig, double* EigVec);
void syDiag(std::complex<double>* Kij, int N, double* Eig, std::complex<double>* EigVec);

// void orthoRandomize(double* elem, int M, int N);
// void matrixSVD(double* Mij_ori, int M, int N, double* U, double* S, double* vT);
#endif /* LAPACK_H */
