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

std::vector<RealType> Ni( const std::vector<Basis> &Bases, const RealVectorType &Vec );
RealMatrixType NiNj( const std::vector<Basis> &Bases, const RealVectorType &Vec );

int main(int argc, char const *argv[]) {
  int L = 3;
  int OBC = false;
  int N = atoi(argv[1]);
  RealType Uin = -10.0;
  RealType Vin = 0.0;
  std::vector<RealType> J;
  for (size_t cnr = 0; cnr < L; cnr++) {
    J.push_back(0.1);
  }
  // J.at(0) *= 1.01;
  J.at(0) = 1.0;

  HDF5IO file("1D_Boson.h5");
  INFO("Build Lattice - ");

  for ( auto &val : J ){
    INFO_NONEWLINE(val << " ");
  }
  INFO("");
  const std::vector< Node<RealType>* > lattice = NN_1D_Chain(L, J, OBC);
  file.saveNumber("1DChain", "L", L);
  file.saveNumber("1DChain", "U", Uin);
  file.saveNumber("1DChain", "V", Vin);
  file.saveStdVector("1DChain", "J", J);
  for ( auto &lt : lattice ){
    if ( !(lt->VerifySite()) ) RUNTIME_ERROR("Wrong lattice setup!");
  }
  INFO("DONE!");
  INFO("Build Basis - ");
  // int N1 = (L+1)/2;
  Basis B1(L, N);
  B1.Boson();
  file.saveNumber("1DChain", "N", N);
  INFO("DONE!");
  INFO_NONEWLINE("Build Hamiltonian - ");
  std::vector<Basis> Bases;
  Bases.push_back(B1);
  Hamiltonian<RealType> ham( Bases );
  std::vector< std::vector<RealType> > Vloc;
  std::vector<RealType> Vtmp(L, Vin);
  Vloc.push_back(Vtmp);
  std::vector< std::vector<RealType> > Uloc;
  std::vector<RealType> Utmp(L, Uin);
  Uloc.push_back(Utmp);
  ham.BuildLocalHamiltonian(Vloc, Uloc, Bases);
  ham.BuildHoppingHamiltonian(Bases, lattice);
  ham.BuildTotalHamiltonian();
  INFO("DONE!");
  INFO_NONEWLINE("Diagonalize Hamiltonian - ");
  std::vector<RealType> Val;
  Hamiltonian<RealType>::VectorType Vec;
  ham.eigh(Val, Vec);
  INFO("GS energy = " << Val.at(0));
  file.saveVector("GS", "EVec", Vec);
  file.saveStdVector("GS", "EVal", Val);
  std::vector< std::vector<int> > st = B1.getBStates();
  std::vector< RealType > tg = B1.getBTags();
  for (size_t cnt = 0; cnt < tg.size(); cnt++) {
    INFO_NONEWLINE( std::setw(3) << cnt << " - ");
    for (auto &j : st.at(cnt)){
      INFO_NONEWLINE(j << " ");
    }
    INFO("- " << Vec(cnt));
  }
  INFO("DONE!");
  std::vector<RealType> Nbi = Ni( Bases, Vec );
  for (auto &n : Nbi ){
    INFO( n << " " );
  }
  std::cout << "" << std::endl;
  RealMatrixType Nij = NiNj( Bases, Vec );
  INFO(Nij);
  std::cout << "" << std::endl;
  INFO(Nij.diagonal());
  file.saveStdVector("Obs", "Nb", Nbi);
  file.saveMatrix("Obs", "Nij", Nij);
  return 0;
}

std::vector<RealType> Ni( const std::vector<Basis> &Bases,
  const RealVectorType &Vec ){
  std::vector<RealType> tmp(Bases.at(0).getL(), 0.0e0);
  std::vector< std::vector<int> > b = Bases.at(0).getBStates();
  assert( b.size() == Vec.size() );
  int coff = 0;
  for ( auto &nbi : b ){
    for (size_t cnt = 0; cnt < Bases.at(0).getL(); cnt++) {
      tmp.at(cnt) += (RealType)nbi.at(cnt) * Vec(coff) * Vec(coff);
    }
    coff++;
  }
  return tmp;
}

RealMatrixType NiNj( const std::vector<Basis> &Bases,
  const RealVectorType &Vec ){
  RealMatrixType tmp = RealMatrixType::Zero(Bases.at(0).getL(), Bases.at(0).getL());
  std::vector< std::vector<int> > b = Bases.at(0).getBStates();
  assert( b.size() == Vec.size() );
  int coff = 0;
  for ( auto &nbi : b ){
    for (size_t cnt1 = 0; cnt1 < Bases.at(0).getL(); cnt1++) {
      for (size_t cnt2 = 0; cnt2 < Bases.at(0).getL(); cnt2++) {
        tmp(cnt1, cnt2) += (RealType)nbi.at(cnt1) * (RealType)nbi.at(cnt2) *
          Vec(coff) * Vec(coff);
      }
    }
    coff++;
  }
  return tmp;
}
