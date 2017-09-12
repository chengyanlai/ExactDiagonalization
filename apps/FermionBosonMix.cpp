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
  const int BL = 2;
  const int BN = 2;
  const int maxLocalB = 1;
  Basis Bosons(BL, BN);
  Bosons.Boson(maxLocalB);
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
  // Local potential
  std::vector<std::vector<double> > Vls;
  std::vector<double> Vloc(BL, 0.0e0);
  Vls.push_back(Vloc);
  Vloc.assign(FL, 0.0e0);
  Vls.push_back(Vloc);
  // Local interaction
  //    Here, U for fermion is NN density density terms.
  std::vector<std::vector<double> > Uls;
  std::vector<double> Uloc(BL, 0.0e0);
  Uls.push_back(Uloc);
  Uloc.assign(FL, 1.0e0);
  Uls.push_back(Uloc);
  Hamiltonian<double> Ham(Bs);
  Ham.BuildLocalHamiltonian(Vls, Uls, Bs);
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
    FJ.push_back(0.10);
  }
  const std::vector< Node<RealType>* > FLattice = NN_1D_Chain(FL, FJ, FOBC);
  std::vector< std::vector< Node<RealType>* > > LT;
  LT.push_back(BLattice);
  LT.push_back(FLattice);
  Ham.BuildHoppingHamiltonian(Bs, LT);
  Ham.BuildTotalHamiltonian();
  /* Print to check. Only for small matrix */
  RealSparseMatrixType w = Ham.getTotalHamiltonian();
  RealMatrixType h(w);
  std::cout << h << std::endl;
  const double DeltaDC = 0.20e0;
  std::vector< std::tuple<int, int, double> > DeltaTerm;
  for ( size_t i = 0; i < FL; i++){
    for (size_t j = 0; j < BL; j++){
      DeltaTerm.push_back(std::make_tuple(j, i, DeltaDC));
    }
  }
  Ham.AddHybridHamiltonian( 0, 1, DeltaTerm, Bs, maxLocalB);
  /* Print to check. Only for small matrix */
  RealSparseMatrixType w2 = Ham.getTotalHamiltonian();
  RealMatrixType h2(w2);
  std::cout << h2 << std::endl;
}
