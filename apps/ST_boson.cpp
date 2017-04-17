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
#define NumCores 2
#endif

void LoadParameters( const std::string filename, int &L, int &N, 
  RealType &jaa, RealType &jab, RealType &Uloc, RealType &Vloc){
    HDF5IO file(filename);
    L = file.loadInt("Parameters", "L");
    N = file.loadInt("Parameters", "N");
    jab = file.loadReal("Parameters", "jab");
    jaa = file.loadReal("Parameters", "jaa");
    Uloc = file.loadReal("Parameters", "U");
    Vloc = file.loadReal("Parameters", "V");
}

std::vector<ComplexType> Ni( const std::vector<Basis> &Bases,
  const ComplexVectorType &Vec ){
  std::vector<ComplexType> tmp(Bases.at(0).getL(), 0.0e0);
  std::vector< std::vector<int> > b = Bases.at(0).getBStates();
  assert( b.size() == Vec.size() );
  int coff = 0;
  for ( auto &nbi : b ){
    for (size_t cnt = 0; cnt < Bases.at(0).getL(); cnt++) {
      tmp.at(cnt) += (RealType)nbi.at(cnt) * Vec(coff) * std::conj(Vec(coff));
    }
    coff++;
  }
  return tmp;
}

ComplexMatrixType NiNj( const std::vector<Basis> &Bases,
  const ComplexVectorType &Vec ){
  ComplexMatrixType tmp = ComplexMatrixType::Zero(Bases.at(0).getL(), Bases.at(0).getL());
  std::vector< std::vector<int> > b = Bases.at(0).getBStates();
  assert( b.size() == Vec.size() );
  int coff = 0;
  for ( auto &nbi : b ){
    for (size_t cnt1 = 0; cnt1 < Bases.at(0).getL(); cnt1++) {
      for (size_t cnt2 = 0; cnt2 < Bases.at(0).getL(); cnt2++) {
        tmp(cnt1, cnt2) += (RealType)nbi.at(cnt1) * (RealType)nbi.at(cnt2) *
          Vec(coff) * std::conj(Vec(coff));
      }
    }
    coff++;
  }
  return tmp;
}

int main(int argc, char const *argv[]) {
  Eigen::setNbThreads(NumCores);
#ifdef MKL
  mkl_set_num_threads(NumCores);
#endif
  INFO("Eigen3 uses " << Eigen::nbThreads() << " threads.");
  const int OBC = 0;
  int L, N, Tsteps, TBloc;
  RealType jaa, jab, Uin, Vin, dt, gamma;
  LoadParameters( "conf.h5", L, N, jaa, jab, Uin, Vin);
  INFO("Build Lattice - ");
  std::vector<ComplexType> JAB, JAA;
  JAB = std::vector<ComplexType>(L, ComplexType(jab, 0.0));
  JAA = std::vector<ComplexType>(L/2, ComplexType(jaa, 0.0));
  INFO("");
  const std::vector< Node<ComplexType>* > lattice = SawTooth(L, JAB, JAA, OBC);
  for ( auto &lt : lattice ){
    if ( !(lt->VerifySite()) ) RUNTIME_ERROR("Wrong lattice setup!");
  }
  INFO("DONE!");
  INFO("Build Basis - ");
  Basis B1(L, N);
  B1.Boson();
  INFO("DONE!");
  INFO_NONEWLINE("Build Hamiltonian - ");
  std::vector<Basis> Bases;
  Bases.push_back(B1);
  Hamiltonian<ComplexType> ham( Bases );
  std::vector< std::vector<ComplexType> > Vloc;
  std::vector<ComplexType> Vtmp(L, (ComplexType)Vin);
  Vloc.push_back(Vtmp);
  std::vector< std::vector<ComplexType> > Uloc;
  std::vector<ComplexType> Utmp(L, (ComplexType)Uin);
  Uloc.push_back(Utmp);
  ham.BuildLocalHamiltonian(Vloc, Uloc, Bases);
  ham.BuildHoppingHamiltonian(Bases, lattice);
  ham.BuildTotalHamiltonian();
  INFO("DONE!");
  INFO_NONEWLINE("Diagonalize Hamiltonian - ");
  std::vector<RealType> Val;
  Hamiltonian<ComplexType>::VectorType Vec;
  ham.eigh(Val, Vec);
  INFO("GS energy = " << Val.at(0));
  HDF5IO file("GS.h5");
  file.saveVector("GS", "EVec", Vec);
  file.saveStdVector("GS", "EVal", Val);
  INFO("DONE!");
  std::vector<ComplexType> Nbi = Ni( Bases, Vec );
  for (auto &n : Nbi ){
    INFO( n << " " );
  }
  ComplexMatrixType Nij = NiNj( Bases, Vec );
  INFO(Nij);
  INFO(Nij.diagonal());
  file.saveStdVector("Obs", "Nb", Nbi);
  file.saveMatrix("Obs", "Nij", Nij);
  return 0;
}
