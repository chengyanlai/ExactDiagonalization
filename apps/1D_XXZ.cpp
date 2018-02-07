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

// #define DTYPE 0//comment out this to use complex

#ifndef DTYPE
#define DT ComplexType
#define DTV ComplexVectorType
#define DTM ComplexMatrixType
#else
#define DT RealType
#define DTV RealVectorType
#define DTM RealMatrixType
#endif

void LoadParameters( const std::string filename, int &L, RealType &Delta,
  int &OBC, int &N, int &dynamics, int &Tsteps, RealType &dt){
    HDF5IO file(filename);
    L = file.loadInt("Parameters", "L");
    Delta = file.loadReal("Parameters", "Delta");
    OBC = file.loadInt("Parameters", "OBC");
    N = file.loadInt("Parameters", "N");
    dynamics = file.loadInt("Parameters", "dynamics");
    Tsteps = file.loadInt("Parameters", "Tsteps");
    dt = file.loadReal("Parameters", "dt");
}

std::vector< DTV > Sz( const std::vector<Basis> &Bases, const DTV &Vec,
  Hamiltonian<DT> &ham ){}
DTM SzSz( const std::vector<Basis> &Bases, const DTV &Vec,
  Hamiltonian<DT> &ham ){}

int main(int argc, char const *argv[]) {
  Eigen::setNbThreads(NumCores);
#ifdef MKL
  mkl_set_num_threads(NumCores);
#endif
  INFO("Eigen3 uses " << Eigen::nbThreads() << " threads.");
  int L;
  RealType Delta;
  int OBC;
  int N;
  int dynamics, Tsteps;
  RealType dt;
  LoadParameters( "conf.h5", L, Delta, OBC, N, dynamics, Tsteps, dt);
  HDF5IO *file = new HDF5IO("XXZ.h5");
  INFO("Build Lattice - ");
  std::vector<DT> Jxy;
  if ( OBC ){
    Jxy = std::vector<DT>(L - 1, 1.0e0);
  } else{
    Jxy = std::vector<DT>(L, 1.0e0);
  }
  for ( auto &val : Jxy ){
    INFO_NONEWLINE(val << " ");
  }
  INFO("");
  const std::vector< Node<DT>* > lattice = NN_1D_Chain(L, Jxy, OBC);
  file->saveStdVector("1DChain", "Jxy", Jxy);
  for ( auto &lt : lattice ){
    if ( !(lt->VerifySite()) ) RUNTIME_ERROR("Wrong lattice setup!");
  }
  INFO("DONE!");
  INFO("Build Basis - ");
  // int N1 = (L+1)/2;
  Basis S12(L, N, true);
  S12.SpinOneHalf();
  std::vector<int> st1 = S12.GetFStates();
  std::vector<size_t> tg1 = S12.GetFTags();
  file->saveNumber("Basis", "Sup", N);
  file->saveStdVector("Basis", "SStates", st1);
  file->saveStdVector("Basis", "STags", tg1);
  INFO("DONE!");
  INFO_NONEWLINE("Build Hamiltonian - ");
  std::vector<Basis> Bases;
  Bases.push_back(S12);
  Hamiltonian<DT> ham( Bases );
  ham.BuildXXZHamiltonian(Delta, Bases, lattice);
  INFO(" - BuildSpinHamiltonian DONE!");
  ham.BuildTotalHamiltonian();
  INFO("DONE!");
  INFO_NONEWLINE("Diagonalize Hamiltonian - ");
  std::vector<RealType> Val;
  Hamiltonian<DT>::VectorType Vec;
  ham.eigh(Val, Vec);
  INFO("GS energy = " << Val.at(0));
  file->saveVector("GS", "EVec", Vec);
  file->saveStdVector("GS", "EVal", Val);
  INFO("DONE!");
  delete file;
  return 0;
}

