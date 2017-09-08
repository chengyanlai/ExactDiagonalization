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
  /* Basis */
  std::vector<Basis> Bs;
  // For bosons
  const int BL = 8;
  const int BN = 4;

  Basis Bosons(BL, BN);
  Bosons.Boson(0, 1);
  Bs.push_back(Bosons);
  // testing
  std::vector< std::vector<int> > BStates = Bosons.getBStates();
  std::vector<RealType> BTags = Bosons.getBTags();
  for( size_t i = 0; i < BTags.size(); ++i ){
    std::cout << BTags.at(i) << ": " << std::flush;
    Bosons.printBosonBasis( BStates.at(i) );
  }

  // For fermions
  const int FL = 2;
  const int FN = 2;
  Basis Fermions(FL, FN, true);
  Fermions.Fermion(0);
  Bs.push_back(Fermions);
  // testing
  std::vector<int> FStates = Fermions.getFStates();
  std::vector<size_t> FTags = Fermions.getFTags();
  for( size_t i = 0; i < FTags.size(); ++i ){
    std::cout << FTags.at(i) << ": " << std::flush;
    Fermions.printFermionBasis( FStates.at(i) );
    std::cout << std::endl;
  }

  /* Hamiltonian */
  // For bosonic chain
  const bool BOBC = false;
  std::vector<double> BJ;
  for (size_t cnt = 0; cnt < BL; cnt++) {
    BJ.push_back(1.0);
  }
  const std::vector< Node<RealType>* > BLattice = NN_1D_Chain(BL, BJ, BOBC);
  // For fermionic chain
  const bool FOBC = true;
  std::vector<double> FJ;
  for (size_t cnt = 0; cnt < FL-1; cnt++) {
    FJ.push_back(1.0);
  }
  const std::vector< Node<RealType>* > FLattice = NN_1D_Chain(FL, FJ, FOBC);
}
