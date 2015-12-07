#include <cmath>
#ifdef MKL
    #include "mkl.h"
#else
    #include "src/numeric/lapack_wrapper.h"
#endif
#include "src/Lanczos/krylov.hpp"

#ifndef DEBUG
#define DEBUG 4
#endif

void krylov(const ComplexSparseMatrixType &A, ComplexVectorType &Vec,
  const ComplexType Prefactor, const size_t Kmax)
{
  const RealType beta_err = 1.0E-12;
  if (DEBUG) assert( Kmax > 2 );
  RealType alpha;
  RealType beta = 1.0;
  //NOTE: Setup the memory
  ComplexMatrixType Vm = ComplexMatrixType::Zero(Vec.size(), Kmax);
  std::vector<RealType> Alphas;
  std::vector<RealType> Betas;
  int cntK = 0;
  //NOTE: normalized Vec
  Vec.normalize();
  Vm.col(cntK) = Vec;
  while ( cntK < Kmax ) {
    ComplexVectorType work = A * Vm.col(cntK);
    if( cntK > 0 ) work -= beta * Vm.col(cntK-1);
    ComplexType alpha_c = work.dot( Vm.col(cntK) );
    alpha = alpha_c.real();
    work -= alpha * Vm.col(cntK);
    beta = work.norm();
    Alphas.push_back(alpha);
    if( DEBUG > 4 ){
      INFO("@ " << cntK << " alpha is " << alpha << " beta is " << beta);
    }
    if( beta > beta_err ){
      work.normalize();
      if( cntK+1 < Kmax ) {
        Vm.col(cntK+1) = work;
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
    assert( Alphas.size() == Betas.size() + 1 );
    if ( DEBUG > 5 ) {
      ComplexMatrixType tmpVm = Vm.adjoint();
      INFO( "Vm^H * Vm " << std::endl << tmpVm * Vm);
    }
  }
  int Kused = cntK;
  if ( Kused == 1 ){
    /* NOTE: This is a special case that input vector is eigenvector */
    Vec = exp(Prefactor) * Vec;
  } else{
    RealType* d = &Alphas[0];
    RealType* e = &Betas[0];
    RealType* z = (RealType*)malloc(Kused * Kused * sizeof(RealType));
    RealType* worktmp = (RealType*)malloc(4 * Kused * sizeof(RealType));
    int info;
    //dstev - LAPACK
    dstev((char*)"V", &Kused, d, e, z, &Kused, worktmp, &info);
    if(info != 0){
      INFO("Lapack INFO = " << info);
      RUNTIME_ERROR("Error in Lapack function 'dstev'");
    }
    ComplexMatrixType Dmat = expD(Prefactor, Kused, d);
    Eigen::Map<RealMatrixType> Kmat(z, Kused, Kused);
    RealMatrixType tmpKmat = Kmat.adjoint();
    if ( DEBUG > 5 ) {
      INFO( "K^H * K " << std::endl << tmpKmat * Kmat);
    }
    ComplexMatrixType tmp = Vm * Kmat;
    Vec = tmp * Dmat * tmp.adjoint() * Vec;
    free(worktmp);
    free(z);
    /* NOTE: Can not free e and d because they are pointing to vector. */
    // free(e);
    // free(d);
  }
}

ComplexMatrixType expD( const ComplexType Prefactor, const size_t dim,
  const RealType* d )
{
  ComplexMatrixType Dmat = ComplexMatrixType::Zero(dim, dim);
  for (size_t cnt = 0; cnt < dim; cnt++) {
    Dmat(cnt, cnt) = exp(Prefactor) * d[cnt];
  }
  return Dmat;
}
