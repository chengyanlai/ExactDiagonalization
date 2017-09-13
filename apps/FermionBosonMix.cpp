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

void density( const std::vector<Basis> &bs, const RealVectorType GS, Hamiltonian<double> Ham,
  std::vector<double> &nb, std::vector<double> &nf){
  nb.clear();
  nb.assign(bs.at(0).getL(), 0.0e0);
  nf.clear();
  nf.assign(bs.at(1).getL(), 0.0e0);
  size_t bid = 0;
  for ( std::vector<int> b : bs.at(0).getBStates() ){
    size_t fid = 0;
    for ( int f : bs.at(1).getFStates() ){
      for ( size_t fl = 0; fl < bs.at(1).getL(); fl++ ){
        if ( btest(f, fl) ){
          size_t id = Ham.DetermineTotalIndex(vec<size_t>(bid, fid));
          nf.at(fl) += GS(id) * GS(id);
        }
      }
      for ( size_t bl = 0; bl < bs.at(0).getL(); bl++ ){
        if ( b.at(bl) > 0 ){
          size_t id = Ham.DetermineTotalIndex(vec<size_t>(bid, fid));
          nb.at(bl) += GS(id) * GS(id);
        }
      }
      fid++;
    }
    bid++;
  }
}

int main(int argc, char const *argv[]) {
  /* Parameters */
  double Jbb = 2.00e0;
  double Jff = 10.00e0;
  double Vbb = 0.010e0;
  double Vff = 0.50e0;
  double Uff = 0.00e0;
  double DeltaDC = 0.20e0;
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
  std::vector<double> Vloc(BL, Vbb);
  Vls.push_back(Vloc);
  Vloc.assign(FL, Vff);
  Vls.push_back(Vloc);
  // Local interaction
  //    Here, U for fermion is NN density density terms.
  std::vector<std::vector<double> > Uls;
  std::vector<double> Uloc(BL, 0.0e0);
  Uls.push_back(Uloc);
  Uloc.assign(FL, Uff);
  Uls.push_back(Uloc);
  Hamiltonian<double> Ham(Bs);
  Ham.BuildLocalHamiltonian(Vls, Uls, Bs);
  // For bosonic chain
  bool BOBC = false;
  std::vector<double> BJ;
  for (size_t cnt = 0; cnt < BL; cnt++) {
    BJ.push_back(Jbb);
  }
  if ( BL == 2 ) {
    BOBC = true;
    BJ.pop_back();
  }
  const std::vector< Node<RealType>* > BLattice = NN_1D_Chain(BL, BJ, BOBC);
  // For fermionic chain
  const bool FOBC = true;
  std::vector<double> FJ;
  for (size_t cnt = 0; cnt < FL-1; cnt++) {
    FJ.push_back(Jff);
  }
  const std::vector< Node<RealType>* > FLattice = NN_1D_Chain(FL, FJ, FOBC);
  std::vector< std::vector< Node<RealType>* > > LT;
  LT.push_back(BLattice);
  LT.push_back(FLattice);
  Ham.BuildHoppingHamiltonian(Bs, LT);
  Ham.BuildTotalHamiltonian();
  /* Print to check. Only for small matrix */
  // RealSparseMatrixType w = Ham.getTotalHamiltonian();
  // RealMatrixType h(w);
  // std::cout << h << std::endl;
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

  // Diagonalization and get GS
  RealVectorType EigVal;
  RealMatrixType EigVec;
  Ham.diag(EigVal, EigVec);
  std::cout << EigVal(0) << " " << EigVal(1) << std::endl;
  RealVectorType GS = EigVec.row(0);
  // compare arpack
  std::vector<RealType> Val;
  RealVectorType Vec;
  Ham.eigh(Val, Vec);
  // test only
  // std::cout << "GS - " << std::endl;
  // for ( size_t i = 0; i < GS.rows(); i++){
  //   std::cout << GS(i) << " " << Vec(i) << std::endl;
  // }

  std::vector<double> Nf, Nb;
  density( Bs, GS, Ham, Nb, Nf);
  std::cout << "Nb - " << std::endl;
  for( auto val : Nb){
    std::cout << val << " " << std::flush;
  }
  std::cout << "\nNf - " << std::endl;
  for( auto val : Nf){
    std::cout << val << " " << std::flush;
  }
  std::cout << std::endl;
}
