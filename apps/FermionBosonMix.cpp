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


int main(int argc, char const *argv[]) {
  const int FL = 1;
  const bool FOBC = true;
  const int BL = 8;
  const int BN = 4;
  const bool BOBC = false;
  std::vector<double> BJ;
  for (size_t cnt = 0; cnt < BL; cnt++) {
    BJ.push_back(1.0);
  }
  // For bosonic chain
  const std::vector< Node<RealType>* > BLattice = NN_1D_Chain(BL, BJ, BOBC);

  // Basis
  std::vector<Basis> Bs;
  Basis Bosons(BL, BN);
  Bosons.Boson(BN-4);

  // testing
  std::vector< std::vector<int> > BStates = Bosons.getBStates();
  std::vector<RealType> BTags = Bosons.getBTags();
  for( size_t i = 0; i < BTags.size(); ++i ){
    std::cout << BTags.at(i) << ": " << std::flush;
    Bosons.printBosonBasis( BStates.at(i) );
  }
  // Basis Fermion(FL, FN, true);
}
