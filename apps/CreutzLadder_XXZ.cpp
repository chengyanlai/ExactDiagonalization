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

void LoadParameters( const std::string filename, int &L,
  RealType &JAA, RealType &JBB, RealType &JvAB, RealType &JdAB,
  RealType &JdBA, RealType &Delta, RealType &Phi,
  int &OBC, int &N, int &dynamics, int &Tsteps, RealType &dt){
    HDF5IO file(filename);
    L = file.loadInt("Parameters", "L");
    JAA = file.loadReal("Parameters", "JAA");
    JBB = file.loadReal("Parameters", "JBB");
    JvAB = file.loadReal("Parameters", "JvAB");
    JdAB = file.loadReal("Parameters", "JdAB");
    JdBA = file.loadReal("Parameters", "JdBA");
    Delta = file.loadReal("Parameters", "Delta");
    Phi = file.loadReal("Parameters", "Phi");
    OBC = file.loadInt("Parameters", "OBC");
    N = file.loadInt("Parameters", "N");
    dynamics = file.loadInt("Parameters", "dynamics");
    Tsteps = file.loadInt("Parameters", "Tsteps");
    dt = file.loadReal("Parameters", "dt");
}

ComplexVectorType Sz( const std::vector<Basis> &Bases, const ComplexVectorType &Vec,
  Hamiltonian<ComplexType> &ham ){
  ComplexVectorType tmp = ComplexVectorType::Zero(Bases.at(0).getL());
  std::vector< int > s12 = Bases.at(0).getFStates();
  int coff = 0;
  for ( auto bs : s12 ){
    for (size_t cnt = 0; cnt < Bases.at(0).getL(); cnt++) {
      if ( btest(bs, cnt) ){// spin-up
        tmp(cnt) += 0.5 * std::pow(std::abs(Vec(coff)), 2);
      }else{// spin-down
        tmp(cnt) -= 0.5 * std::pow(std::abs(Vec(coff)), 2);
      }
    }
    coff++;
  }
  return tmp;
}
// ComplexMatrixType SzSz( const std::vector<Basis> &Bases, const ComplexVectorType &Vec,
//   Hamiltonian<ComplexType> &ham ){}

int main(int argc, char const *argv[]) {
  Eigen::setNbThreads(NumCores);
#ifdef MKL
  mkl_set_num_threads(NumCores);
#endif
  INFO("Eigen3 uses " << Eigen::nbThreads() << " threads.");
  const ComplexType I = ComplexType(0.0e0, 1.0e0);
  int L;
  RealType jaa, jbb, jvab, jdab, jdba, delta, phi;
  int OBC;
  int N;
  int dynamics, Tsteps;
  RealType dt;
  LoadParameters( "conf.h5", L, jaa, jbb, jvab, jdab, jdba, delta, phi, OBC, N, dynamics, Tsteps, dt);
  assert( L % 2 == 0);
  int LUnit = L / 2;
  HDF5IO *file = new HDF5IO("CLxxz.h5");
  INFO("Build Lattice - ");
  std::vector<ComplexType> JAA, JBB, JvAB, JdAB, JdBA;
  if ( OBC ){
    JAA = std::vector<ComplexType>(LUnit - 1, jaa * exp(I * phi));
    JBB = std::vector<ComplexType>(LUnit - 1, jbb * exp(-1.0e0 * I * phi));
    JdAB = std::vector<ComplexType>(LUnit - 1, jdab);
    JdBA = std::vector<ComplexType>(LUnit - 1, jdba);
  } else{
    JAA = std::vector<ComplexType>(LUnit, jaa * exp(I * phi));
    JBB = std::vector<ComplexType>(LUnit, jbb * exp(-1.0e0 * I * phi));
    JdAB = std::vector<ComplexType>(LUnit, jdab);
    JdBA = std::vector<ComplexType>(LUnit, jdba);
  }
  // for ( auto val : JAA ) std::cout << val << std::endl;
  JvAB = std::vector<ComplexType>(LUnit, jvab);
  INFO("");
  const std::vector< Node<ComplexType>* > lattice = CreutzLadder(LUnit, JAA, JBB, JdAB, JdBA, JvAB, OBC);
  for ( auto &lt : lattice ){
    if ( !(lt->VerifySite()) ) RUNTIME_ERROR("Wrong lattice setup!");
  }
  INFO("DONE!");
  INFO("Build Basis - ");
  // int N1 = (L+1)/2;
  Basis S12(L, N, true);
  S12.SpinOneHalf();
  std::vector<int> st1 = S12.getFStates();
  std::vector<size_t> tg1 = S12.getFTags();
  file->saveNumber("Basis", "Sup", N);
  file->saveStdVector("Basis", "SStates", st1);
  file->saveStdVector("Basis", "STags", tg1);
  INFO("DONE!");
  INFO_NONEWLINE("Build Hamiltonian - ");
  std::vector<Basis> Bases;
  Bases.push_back(S12);
  Hamiltonian<ComplexType> ham( Bases );
  ham.BuildXXZHamiltonian(delta, Bases, lattice);
  INFO(" - BuildSpinHamiltonian DONE!");
  ham.BuildTotalHamiltonian();
  INFO("DONE!");
  INFO_NONEWLINE("Diagonalize Hamiltonian - ");
  std::vector<RealType> Val;
  Hamiltonian<ComplexType>::VectorType Vec;
  ham.eigh(Val, Vec);
  INFO("GS energy = " << Val.at(0));
  file->saveVector("GS", "EVec", Vec);
  file->saveStdVector("GS", "EVal", Val);
  ComplexVectorType Szi = Sz( Bases,Vec, ham );
  file->saveVector("obs", "Sz", Szi);
  INFO("DONE!");
  delete file;
  return 0;
}

