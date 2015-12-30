#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>//sqrt
#include "src/EDType.hpp"
#include "src/Node/Node.hpp"
#include "src/Lattice/preset.hpp"
#include "src/Basis/Basis.hpp"
#include "src/Hamiltonian/Hamiltonian.hpp"

std::vector<RealType> Ni( const std::vector<Basis> &Bases, const RealVectorType &Vec );
RealMatrixType NiNj( const std::vector<Basis> &Bases, const RealVectorType &Vec );

int main(int argc, char const *argv[]) {
  const int L = 6;
  INFO("Build Lattice - ");
  // const bool OBC = true;
  // const std::vector<RealType> J(L-1, 1.0);
  // const std::vector< Node<RealType, int>* > lattice = NN_1D_Chain(L, J, OBC);
  const bool OBC = false;
  std::vector<RealType> JAB(L, -1.0e0*sqrt(2.0e0));
  JAB.at(0) = JAB.at(0) * (1.0e0 + 0.010e0);
  JAB.at(1) = JAB.at(1) * (1.0e0 + 0.010e0);
  std::vector<RealType> JAA(L/2, -1.0e0);
  // JAA.at(0) = JAA.at(0) * (1.0e0 + 0.010e0);
  const std::vector< Node<RealType, int>* > lattice = SawTooth(L, JAB, JAA, OBC);
  for ( auto &j : lattice ){
    INFO_NONEWLINE(j->data << " " << j->VerifySite() << " " << j->NumNeighbors());
    INFO_NONEWLINE(" Neighbors: ");
    for ( auto &p : j->getNeighbors() ){
      INFO_NONEWLINE( p->data << " ");
    }
    INFO(" ");
  }
  INFO("DONE!");
  INFO("Build Basis - ");
  const int N = 2;
  Basis B1(L, N);
  B1.Boson();
  std::vector< std::vector<int> > st = B1.getBStates();
  std::vector< RealType > tg = B1.getBTags();
  for (size_t cnt = 0; cnt < tg.size(); cnt++) {
    INFO_NONEWLINE( std::setw(3) << cnt << " - ");
    for (auto &j : st.at(cnt)){
      INFO_NONEWLINE(j << " ");
    }
    INFO("- " << tg.at(cnt));
  }
  INFO("DONE!");
  INFO_NONEWLINE("Build Hamiltonian - ");
  std::vector<Basis> Bases;
  Bases.push_back(B1);
  Hamiltonian<RealType,int> ham( Bases );
  std::vector< std::vector<RealType> > Vloc;
  std::vector<RealType> Vtmp(L, 0.0);
  Vloc.push_back(Vtmp);
  std::vector< std::vector<RealType> > Uloc;
  std::vector<RealType> Utmp(L, 0.0);
  Uloc.push_back(Utmp);
  ham.BuildLocalHamiltonian(Vloc, Uloc, Bases);
  ham.BuildHoppingHamiltonian(Bases, lattice);
  ham.BuildTotalHamiltonian();
  INFO("DONE!");
  INFO_NONEWLINE("Diagonalize Hamiltonian - ");
  RealType Val = 0.0;
  Hamiltonian<RealType,int>::VectorType Vec;
  ham.eigh(Val, Vec);
  INFO("GS energy = " << Val);
  INFO("DONE!");
  std::vector<RealType> Nbi = Ni( Bases, Vec );
  for (auto &n : Nbi ){
    INFO( n << " " );
  }
  RealMatrixType Nij = NiNj( Bases, Vec );
  INFO(Nij);
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
        tmp(cnt1, cnt2) += (RealType)nbi.at(cnt1) * (RealType)nbi.at(cnt2) * Vec(coff) * Vec(coff);
      }
    }
    coff++;
  }
  return tmp;
}
