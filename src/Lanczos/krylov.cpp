#include <cmath>
// #include <Eigen/Eigenvalues>
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
  // INFO(A);
  const RealType beta_err = 1.0E-12;
  if (DEBUG) assert( Kmax > 2 );
  RealType alpha;
  RealType beta = 1.0;
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
    assert( Betas.size() == cntK - 1 );
    if ( DEBUG > 5 ) {
      ComplexMatrixType tmpVm = Vm;
      tmpVm.adjointInPlace();
      INFO( "Vm^H * Vm " << std::endl << tmpVm * Vm);
    }
  }
  int Kused = cntK;
  if ( Kused == 1 ){
    /* NOTE: This is a special case that input vector is eigenvector */
    Vec = exp(Prefactor * Alphas.at(0)) * Vec;
  } else{
    /* NOTE: The floowing codes use Eigen to solve the tri-diagonal matrix.
             If we use this, we need the Eigen header in this file.
    */
    // RealMatrixType TriDiag = RealMatrixType::Zero(Kused, Kused);
    // for (size_t cnt = 0; cnt < Kused; cnt++) {
    //   TriDiag(cnt, cnt) = Alphas.at(cnt);
    //   if (cnt > 0) {
    //     TriDiag(cnt, cnt-1) = Betas.at(cnt - 1);
    //     TriDiag(cnt-1, cnt) = Betas.at(cnt - 1);
    //   }
    // }
    // Eigen::SelfAdjointEigenSolver<RealMatrixType> es;
    // es.compute(TriDiag);
    // RealVectorType Dvec = es.eigenvalues();
    // ComplexMatrixType Dmat = ComplexMatrixType::Zero(Kused, Kused);
    // for (size_t cnt = 0; cnt < Kused; cnt++) {
    //   Dmat(cnt,cnt) = exp( Prefactor * Dvec(cnt) );
    // }
    // RealMatrixType Kmat = es.eigenvectors();
    /* NOTE: The floowing codes use MKL Lapack to solve the tri-diagonal matrix */
    RealType* d = &Alphas[0];
    RealType* e = &Betas[0];
    RealType* z = (RealType*)malloc(Kused * Kused * sizeof(RealType));
    RealType* work = (RealType*)malloc(4 * Kused * sizeof(RealType));
    int info;
    //dstev - LAPACK
    dstev((char*)"V", &Kused, d, e, z, &Kused, work, &info);
    Eigen::Map<RealMatrixType> Kmat(z, Kused, Kused);
    Kmat.transposeInPlace();
    ComplexMatrixType Dmat = ComplexMatrixType::Zero(Kused, Kused);
    for (size_t cnt = 0; cnt < Kused; cnt++) {
      Dmat(cnt,cnt) = exp( Prefactor * d[cnt] );
    }
    if(info != 0){
      INFO("Lapack INFO = " << info);
      RUNTIME_ERROR("Error in Lapack function 'dstev'");
    }
    /* NOTE: After Solving tri-diagonal matrix, we need Kmat and Dmat to proceed further. */
    if ( DEBUG > 4 ) {
      RealMatrixType tmpKmat = Kmat;
      tmpKmat.transposeInPlace();
      INFO( Kmat * Dmat * tmpKmat );
    }
    ComplexMatrixType Otmp = Vm.block(0, 0, Vec.size(), Kused) * Kmat;
    Vec = ( Otmp * Dmat ) * ( Otmp.adjoint() * Vec );
  }
}
