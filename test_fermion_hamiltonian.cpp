#include <iostream>
#include <vector>
#include "src/EDType.hpp"
#include "src/Node/Node.hpp"
#include "src/Lattice/preset.hpp"
#include "src/Basis/Basis.hpp"
#include "src/Hamiltonian/Hamiltonian.hpp"

int main(int argc, char const *argv[]) {
  const int L = 4;
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
  int N = 2;
  Basis F1(L, N, true);
  F1.Fermion();
  std::vector<int> st = F1.getFStates();
  std::vector<size_t> tg = F1.getFTags();
  for (size_t cnt = 0; cnt < st.size(); cnt++) {
    INFO(cnt << " " << st.at(cnt) << " " << tg.at(st.at(cnt)));
  }
  INFO("DONE!");
  INFO_NONEWLINE("Build Hamiltonian - ");
  std::vector<Basis> Bases;
  Bases.push_back(F1);
  Bases.push_back(F1);
  Hamiltonian<RealType,int> ham( Bases );
  INFO("DONE!");
  return 0;
}
