#include <iostream>
#include <vector>
#include "src/EDType.hpp"
#include "src/Node/Node.hpp"
#include "src/Lattice/preset.hpp"
#include "src/Basis/Basis.hpp"
#include "src/Hamiltonian/Hamiltonian.hpp"

int main(int argc, char const *argv[]) {
  const int L = 13;
  INFO("Build Lattice - ");
  const bool OBC = true;
  const std::vector<double> J(L-1, 1.0);
  const std::vector< Node<double, int>* > lattice = NN_1D_Chain(L, J, true);
  // const bool OBC = false;
  // const std::vector<double> J(L, 1.0);
  // const std::vector< Node<double, int>* > lattice = NN_1D_Chain(L, J, false);
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
  const int N = 11;
  Basis B1(L, N);
  B1.Boson();
  std::vector< std::vector<int> > st = B1.getBStates();
  std::vector< double > tg = B1.getBTags();
  // for (size_t cnt = 0; cnt < tg.size(); cnt++) {
  //   INFO_NONEWLINE(cnt << " ");
  //   for (auto &j : st.at(cnt)){
  //     INFO_NONEWLINE(j << " ");
  //   }
  //   INFO(tg.at(cnt));
  // }
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
  ham.BuildIntraLocalHamiltonian(Vloc, Uloc, Bases);
  ham.BuildIntraHoppingHamiltonian(Bases, lattice);
  ham.BuildTotalHamiltonian();
  ham.eigh();
  INFO("DONE!");
  return 0;
}
