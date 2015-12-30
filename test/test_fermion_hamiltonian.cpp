#include <iostream>
#include <iomanip>
#include <vector>
#include "src/EDType.hpp"
#include "src/bitwise.h"
#include "src/Node/Node.hpp"
#include "src/Lattice/preset.hpp"
#include "src/Basis/Basis.hpp"
#include "src/Hamiltonian/Hamiltonian.hpp"

std::vector<RealType> HarmonicPotential( const int L, const RealType alpha );
std::vector< std::vector<RealType> > Ni( const std::vector<Basis> &Bases,
  const RealVectorType &Vec, Hamiltonian<RealType,int> &ham );

int main(int argc, char const *argv[]) {
  const int L = 13;
  INFO("Build Lattice - ");
  const bool OBC = true;
  std::vector<RealType> J(L-1, 1.0);
  for (size_t cnt = 0; cnt < L-1; cnt+=2) {
    J.at(cnt) *= 0.1;
  }
  // const bool OBC = false;
  // const std::vector<double> J(L, 1.0);
  const std::vector< Node<RealType, int>* > lattice = NN_1D_Chain(L, J, OBC);
  // for ( auto &j : lattice ){
  //   INFO_NONEWLINE(j->data << " " << j->VerifySite() << " " << j->NumNeighbors());
  //   INFO_NONEWLINE(" Neighbors: ");
  //   for ( auto &p : j->getNeighbors() ){
  //     INFO_NONEWLINE( p->data << " ");
  //   }
  //   INFO(" ");
  // }
  INFO("DONE!");
  INFO("Build Basis - ");
  int N1 = (L+1)/2;
  Basis F1(L, N1, true);
  F1.Fermion();
  std::vector<int> st1 = F1.getFStates();
  std::vector<size_t> tg1 = F1.getFTags();
  INFO("Species - 1");
  // for (size_t cnt = 0; cnt < st1.size(); cnt++) {
  //   INFO_NONEWLINE( std::setw(3) << st1.at(cnt) << " - ");
  //   F1.printFermionBasis(st1.at(cnt));
  //   INFO("- " << tg1.at(st1.at(cnt)));
  // }
  int N2 = (L-1)/2;
  Basis F2(L, N2, true);
  F2.Fermion();
  std::vector<int> st2 = F2.getFStates();
  std::vector<size_t> tg2 = F2.getFTags();
  INFO("Species - 2");
  // for (size_t cnt = 0; cnt < st2.size(); cnt++) {
  //   INFO_NONEWLINE( std::setw(3) << st2.at(cnt) << " - ");
  //   F2.printFermionBasis(st2.at(cnt));
  //   INFO("- " << tg2.at(st2.at(cnt)));
  // }
  INFO("DONE!");
  INFO_NONEWLINE("Build Hamiltonian - ");
  std::vector<Basis> Bases;
  Bases.push_back(F1);
  Bases.push_back(F2);
  Hamiltonian<RealType,int> ham( Bases );
  std::vector< std::vector<RealType> > Vloc;
  std::vector<RealType> Vtmp = HarmonicPotential( L, 10.0e0 );
  // std::vector<RealType> Vtmp(L, 1.0);
  Vloc.push_back(Vtmp);
  Vloc.push_back(Vtmp);
  std::vector< std::vector<RealType> > Uloc;
  std::vector<RealType> Utmp(L, 10.0);
  Uloc.push_back(Utmp);
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
  std::vector< std::vector<RealType> > Nfi = Ni( Bases, Vec, ham );
  INFO(" Up Spin - ");
  for ( auto &b : Nfi.at(0) ){
    INFO( b << " ");
  }
  INFO(" Down Spin - ");
  for ( auto &b : Nfi.at(1) ){
    INFO( b << " ");
  }
  return 0;
}

std::vector<RealType> HarmonicPotential( const int L, const RealType alpha ){
  RealType center;
  center = (RealType)(L - 1) / 2.0e0;
  std::vector<RealType> v(L, 0.0e0);
  for (size_t l = 0; l < L; l++) {
    v.at(l) = 0.50e0 * alpha *
      pow( std::abs(center - (RealType)l), 2 ) *
      pow(1.0e0/(RealType)L, 2);
    // INFO(v.at(l));
  }
  return v;
}

std::vector< std::vector<RealType> > Ni( const std::vector<Basis> &Bases,
  const RealVectorType &Vec, Hamiltonian<RealType,int> &ham ){
  std::vector< std::vector<RealType> > out;
  std::vector<RealType> tmp1(Bases.at(0).getL(), 0.0e0);
  std::vector<RealType> tmp2(Bases.at(1).getL(), 0.0e0);
  std::vector< int > f1 = Bases.at(0).getFStates();
  std::vector< int > f2 = Bases.at(1).getFStates();
  size_t f1id = 0, f2id = 0;
  for ( int &nf2 : f2 ){
    std::vector<size_t> ids(2,f2id);
    f1id = 0;
    for ( int &nf1 : f1 ){
      ids.at(0) = f1id;
      size_t id = ham.DetermineTotalIndex(ids);
      for (size_t cnt = 0; cnt < Bases.at(0).getL(); cnt++) {
        if ( btest(nf1, cnt) ) tmp1.at(cnt) += Vec(id) * Vec(id);
      }
      for (size_t cnt = 0; cnt < Bases.at(1).getL(); cnt++) {
        if ( btest(nf2, cnt) ) tmp2.at(cnt) += Vec(id) * Vec(id);
      }
      f1id++;
    }
    f2id++;
  }
  out.push_back(tmp1);
  out.push_back(tmp2);
  return out;
}
