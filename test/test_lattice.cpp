#include <iostream>
#include <vector>
#include "src/EDType.hpp"
#include "src/Lattice/preset.hpp"

int main(int argc, char const *argv[]) {
  int L = 4;
  bool OBC = true;
  std::vector<double> J(L-1, 1.0);
  std::vector< Node<double, int>* > lattice = NN_1D_Chain(L, J, true);
  // bool OBC = false;
  // std::vector<double> J(L, 1.0);
  // std::vector< Node<double, int>* > lattice = NN_1D_Chain(L, J, false);
  for ( auto &j : lattice ){
    INFO(j->VerifySite());
  }
}
