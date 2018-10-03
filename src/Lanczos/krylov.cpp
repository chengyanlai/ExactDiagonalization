#include <cmath>
#include <cassert>
#ifdef MKL
    #include "mkl.h"
#else
    #include "src/numeric/lapack_wrapper.h"
#endif
#include "src/Lanczos/krylov.hpp"

#ifndef DEBUG
#define DEBUG 4
#endif

void krylovEXP(const ComplexSparseMatrixType &A, ComplexVectorType &Vec, const ComplexType Prefactor, const size_t Kmax, const double threshNorm){
  #if defined(MKL)
    mkl_set_num_threads(NumCores);
  #endif
  if (DEBUG) assert( Kmax > 2 );
  RealType alpha;
  RealType beta = 1.0;
  ComplexMatrixType Vm(Vec.size(), Kmax, arma::fill::zeros);
  std::vector<RealType> Alphas;
  std::vector<RealType> Betas;
  int cntK = 0;
  //* Normalized Vec
  Vm.col(cntK) = arma::normalise(Vec);
  while ( cntK < Kmax ) {
    ComplexVectorType work = A * Vm.col(cntK);
    if( cntK > 0 ) work -= beta * Vm.col(cntK-1);
    ComplexType alpha_c = arma::cdot( work,  Vm.col(cntK) );
    alpha = alpha_c.real();
    work -= alpha * Vm.col(cntK);
    beta = arma::norm(work);
    Alphas.push_back(alpha);
    if( DEBUG > 4 ){
      std::cout <<"@ " << cntK << " alpha is " << alpha << " beta is " << beta << std::endl;
    }
    if( beta > threshNorm ){
      // work.normalize();
      if( cntK+1 < Kmax ) {
        Vm.col(cntK+1) = arma::normalise(work);
        Betas.push_back(beta);
      }
      cntK++;
    }
    else{
      cntK++;
      break;
    }
  }
  if ( DEBUG ) {
    assert( Alphas.size() == cntK );
    assert( Betas.size() == cntK - 1 );
    if ( DEBUG > 5 ) {
      ComplexMatrixType tmpVm = Vm;
      // tmpVm.adjointInPlace();
      tmpVm.t();
      std::cout << "Vm^H * Vm " << std::endl << tmpVm * Vm << std::endl;
    }
  }
  int Kused = cntK;
  if ( Kused == 1 ){
    //* NOTE: This is a special case that input vector is eigenvector
    Vec = exp(Prefactor * Alphas.at(0)) * Vec;
  } else{
    //* NOTE: The floowing codes use MKL Lapack to solve the tri-diagonal matrix
    RealType* d = &Alphas[0];
    RealType* e = &Betas[0];
    RealType* z = (RealType*)malloc(Kused * Kused * sizeof(RealType));
    RealType* work = (RealType*)malloc(4 * Kused * sizeof(RealType));
    int info;
    //* dstev - LAPACK
    dstev((char*)"V", &Kused, d, e, z, &Kused, work, &info);
    RealMatrixType Kmat(z, Kused, Kused);
    Kmat.t();
    ComplexMatrixType Dmat(Kused, Kused, arma::fill::zeros);
    for (size_t cnt = 0; cnt < Kused; cnt++) {
      Dmat(cnt,cnt) = exp( Prefactor * d[cnt] );
    }
    if(info != 0) {
      std::cout << "Lapack INFO = " << info << std::endl;
      RUNTIME_ERROR("Error in Lapack function 'dstev'");
    }
    //* NOTE: After Solving tri-diagonal matrix, we need Kmat and Dmat to proceed further.
    if ( DEBUG > 4 ) {
      RealMatrixType tmpKmat = Kmat;
      tmpKmat.t();
      std::cout << Kmat * Dmat * tmpKmat << std::endl;
    }
    Vm.reshape( Vec.size(), Kused );
    ComplexMatrixType Otmp = Vm * Kmat;
    Vec = ( Otmp * Dmat ) * ( Otmp.t() * Vec );
  }
}

template<typename Tnum>
void krylov(const arma::SpMat<Tnum> &A, const arma::Col<Tnum> &Vec, arma::Mat<Tnum> &OTrans, RealVectorType &Dvec, const size_t Kmax, const double threshNorm){
  #if defined(MKL)
    mkl_set_num_threads(NumCores);
  #endif
  if (DEBUG) assert( Kmax > 2 );
  RealType alpha;
  RealType beta = 1.0;
  arma::Mat<Tnum> Vm(Vec.size(), Kmax, arma::fill::zeros);
  std::vector<RealType> Alphas;
  std::vector<RealType> Betas;
  int cntK = 0;
  //* NOTE: normalized Vec
  Vm.col(cntK) = arma::normalise(Vec);
  while ( cntK < Kmax ) {
    arma::Col<Tnum> work = A * Vm.col(cntK);
    if( cntK > 0 ) work -= beta * Vm.col(cntK-1);
    Tnum alpha_c = arma::cdot( work,  Vm.col(cntK) );
    alpha = RealPart(alpha_c);
    work -= alpha * Vm.col(cntK);
    beta = arma::norm(work);
    Alphas.push_back(alpha);
    if( DEBUG > 4 ){
      std::cout <<"@ " << cntK << " alpha is " << alpha << " beta is " << beta << std::endl;
    }
    if( beta > threshNorm ){
      if( cntK+1 < Kmax ) {
        Vm.col(cntK+1) = arma::normalise(work);
        Betas.push_back(beta);
      }
      cntK++;
    }
    else{
      cntK++;
      break;
    }
  }
  if ( DEBUG ) {
    assert( Alphas.size() == cntK );
    assert( Betas.size() == cntK - 1 );
    if ( DEBUG > 5 ) {
      arma::Mat<Tnum> tmpVm = Vm;
      tmpVm.t();
      std::cout << "Vm^H * Vm " << std::endl << tmpVm * Vm << std::endl;
    }
  }
  int Kused = cntK;
  if ( Kused == 1 ){
    Vm.reshape( Vec.n_rows, Kused );
    OTrans = Vm;
    Dvec = RealVectorType(Kused, arma::fill::ones);
  }else{
    //* NOTE: The floowing codes use MKL Lapack to solve the tri-diagonal matrix
    RealType* d = &Alphas[0];
    RealType* e = &Betas[0];
    RealType* z = (RealType*)malloc(Kused * Kused * sizeof(RealType));
    RealType* work = (RealType*)malloc(4 * Kused * sizeof(RealType));
    int info;
    //* dstev - LAPACK
    dstev((char*)"V", &Kused, d, e, z, &Kused, work, &info);
    RealMatrixType Kmat(z, Kused, Kused);
    Kmat.t();
    Dvec = RealVectorType(Kused, arma::fill::zeros);
    for (size_t cnt = 0; cnt < Kused; cnt++) {
      Dvec(cnt) = d[cnt];
    }
    if(info != 0) {
      std::cout <<"Lapack INFO = " << info << std::endl;
      RUNTIME_ERROR("Error in Lapack function 'dstev'");
    }
    Vm.reshape( Vec.n_rows, Kused );
    OTrans = Vm * Kmat;
  }
}
template void krylov(const arma::SpMat<RealType> &A, const arma::Col<RealType> &Vec, arma::Mat<RealType> &OTrans, RealVectorType &Dvec, const size_t Kmax, const double threshNorm);
template void krylov(const arma::SpMat<ComplexType> &A, const arma::Col<ComplexType> &Vec, arma::Mat<ComplexType> &OTrans, RealVectorType &Dvec, const size_t Kmax, const double threshNorm);
