#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include "src/EDType.hpp"
#include "src/bitwise.h"
#include "src/Node/Node.hpp"
#include "src/Lattice/preset.hpp"
#include "src/Basis/Basis.hpp"
#include "src/Hamiltonian/Hamiltonian.hpp"
#include "src/hdf5io/hdf5io.hpp"

#ifdef MKL
  #include "mkl.h"
#endif

#ifndef DEBUG
#define DEBUG 1
#endif

#ifndef NumCores
#define NumCores 20
#endif

#define DT RealType
#define DTV RealVectorType
#define DTM RealMatrixType

// std::vector< DTV > Sz( const std::vector<Basis> &Bases, const DTV &Vec,
//   Hamiltonian<DT> &ham ){}
// DTM SzSz( const std::vector<Basis> &Bases, const DTV &Vec,
//   Hamiltonian<DT> &ham ){}

int main(int argc, char const *argv[]) {
  Eigen::setNbThreads(NumCores);
#ifdef MKL
  mkl_set_num_threads(NumCores);
#endif
  INFO("Eigen3 uses " << Eigen::nbThreads() << " threads.");
  const int L = 10;
  int OBC = 0;
  RealType hz = 0.0;
  const RealType Temp = -2.0e0;
  // const RealType Temp = -100.0e0;
  // HDF5IO *file = new HDF5IO("TIsing.h5");
  INFO("Build Lattice - ");
  std::vector<DT> Jxx;
  if ( OBC ){
    Jxx = std::vector<DT>(L - 1, 1.0e0);
  } else{
    Jxx = std::vector<DT>(L, 1.0e0);
  }
  for ( auto &val : Jxx ){
    INFO_NONEWLINE(val << " ");
  }
  INFO("");
  const std::vector< Node<DT>* > lattice = NN_1D_Chain(L, Jxx, OBC);
  // file->saveStdVector("1DChain", "Jxx", Jxx);
  for ( auto &lt : lattice ){
    if ( !(lt->VerifySite()) ) RUNTIME_ERROR("Wrong lattice setup!");
  }
  INFO("DONE!");
  INFO("Build Basis - ");
  // int N1 = (L+1)/2;
  Basis S12(L, 0, true);
  S12.TIsing();
  std::vector<int> st1 = S12.getFStates();
  std::vector<size_t> tg1 = S12.getFTags();
  INFO("DONE!");
  INFO_NONEWLINE("Build Hamiltonian - ");
  std::vector<Basis> Bases;
  Bases.push_back(S12);
  std::vector<RealType> Mz;
  for ( size_t hzi = 0; hzi <= 150; hzi++){
    hz = (RealType)hzi * 0.01;
  Hamiltonian<DT> ham( Bases );
  ham.BuildTIsingHamiltonian(hz, Bases, lattice);
  INFO(" - BuildSpinHamiltonian DONE!");
  ham.BuildTotalHamiltonian();
  INFO("DONE!");
  // std::cout << ham.getTotalHamiltonian() << std::endl;
  INFO("Diagonalize Hamiltonian - ");
  Hamiltonian<DT>::VectorType Val;
  Hamiltonian<DT>::MatrixType Vec;
  ham.diag(Val, Vec);
  Hamiltonian<DT>::MatrixType eVal = ham.expVals( Temp, Val ).asDiagonal();
  Hamiltonian<DT>::MatrixType SzOp(S12.getHilbertSpace(), S12.getHilbertSpace());
  SzOp.setZero();
  for ( size_t i = 0; i < S12.getHilbertSpace(); i++){
    SzOp(i,i) = S12.getSzTotal(i);
    // std::cout << Val(i) << " " << eVal(i,i) << " " << SzOp(i,i) << std::endl;
  }
  RealType Z0 = (Vec.transpose() * eVal * Vec).trace();
  // std::cout << Z0 << std::endl;
  std::cout << (Vec.transpose() * eVal * Vec * SzOp).trace() / ((RealType)L * Z0) << std::endl;
  Mz.push_back((Vec.transpose() * eVal * Vec * SzOp).trace() / ((RealType)L * Z0));
  INFO("DONE!");
  }
  for (auto i : Mz ){
    std::cout << i << "," << std::flush;
  }
  return 0;
}

