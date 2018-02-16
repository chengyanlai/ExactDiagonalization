#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include "src/EDType.hpp"
#include "src/Node/Node.hpp"
#include "src/Lattice/preset.hpp"
#include "src/Basis/Basis.hpp"
#include "src/Hamiltonian/Holstein/Holstein.hpp"
#include "src/hdf5io/hdf5io.hpp"

#ifdef MKL
  #include "mkl.h"
#endif

ComplexVectorType NPhonon( const std::vector<Basis> &Bases, const ComplexVectorType &Vec, const Holstein<ComplexType>& Ham){
  int L = Bases.at(0).GetL() + 1;
  // int N = Bases.at(0).GetN();
  ComplexVectorType out(L-1, arma::fill::zeros);
  std::vector< std::vector<int> > b = Bases.at(0).GetBStates();
  assert( b.size() * L == Vec.size() );
  for ( size_t cnt = 0; cnt < L; cnt++ ){
    int coff = 0;
    for ( auto &nbi : b ){
      for (size_t cnt = 0; cnt < L-1; cnt++) {
        size_t idx = Ham.DetermineTotalIndex( vec<size_t>(cnt, coff) );
        out.at(cnt) += (RealType)nbi.at(cnt) * std::pow(std::abs(Vec(idx)), 2);
      }
      coff++;
    }
  }
  return out;
}

ComplexVectorType NFermion( const std::vector<Basis> &Bases, const ComplexVectorType &Vec, const Holstein<ComplexType>& Ham){
  int L = Bases.at(0).GetL() + 1;
  // int N = Bases.at(0).GetN();
  ComplexVectorType out(L, arma::fill::zeros);
  assert( Bases.at(0).GetHilbertSpace() * L == Vec.size() );
  for ( size_t cnt = 0; cnt < L; cnt++ ){
    int coff = 0;
    for ( size_t b = 0; b < Bases.at(0).GetHilbertSpace(); b++ ){
      size_t idx = Ham.DetermineTotalIndex( vec<size_t>(cnt, b) );
      out.at(cnt) += std::pow(std::abs(Vec(idx)), 2);
    }
  }
  return out;
}

void Equilibrium(const std::string prefix){
  std::ofstream LogOut;
  LogOut.open(prefix + "Holstein.K.eqm", std::ios::app);
  int L = 4;
  int N = 10 * L;
  const RealType Jin = 1.0;
  RealType Win = 10.0;
  RealType Gin = 10.0;

  /* k points */
  LogOut << "Build k-point - " << std::flush;
  std::vector<int> Kn;
  Kn.clear();
  Kn.push_back(0);
  for ( size_t n = 1; n <= L/2 - 1; n++ ){
    Kn.push_back(+n);
    Kn.push_back(-n);
  }
  Kn.push_back(L/2);
  // PrintVector(Kn, 3, " ");
  assert( Kn.size() == L );
  LogOut << Kn.size() << " k-points DONE!" << std::endl;

  // HDF5IO* file = new HDF5IO("HolsteinK.h5");
  // delete file;
  LogOut << "Build Basis - max. phonon quanta = " << N << ", ends up with Phonon Hilbert space = " << std::flush;
  Basis P1(L-1, N);// Get rid og k=0 phonon mode, so using L - 1!!
  P1.PhononK();
  LogOut << P1.GetHilbertSpace() << std::flush;
  // LogOut << P1 << std::endl;
  std::vector<Basis> Bases;
  Bases.push_back(P1);
  LogOut << " DONE!" << std::endl;

  LogOut << "Build Hamiltonian - " << std::flush;
  Holstein<ComplexType> Ham0( Kn.size(), Bases );
  Ham0.HolsteinModelK(Kn, Bases, Win, Gin, Jin, N);
  LogOut << "Hermitian = " << Ham0.CheckHermitian() << ", Hilbert space = " << Ham0.GetTotalHilbertSpace() << ", DONE!" << std::endl;

  LogOut << "Diagonalize Hamiltonian - " << std::flush;
  RealVectorType Vals;
  ComplexMatrixType Vecs;
  std::string Target = "SR";
  if ( Ham0.GetTotalHilbertSpace() > 1200 ){
    Ham0.eigh(Vals, Vecs, 20, true, Target);
  }else{
    Ham0.diag(Vals, Vecs);// Full spectrum
  }
  LogOut << "DONE!" << std::endl;
  LogOut << "\tGS energy = " << std::setprecision(12) << Vals[0] << std::endl;
  LogOut << "\tFES energy = " << std::setprecision(12) << Vals[1] << std::endl;
  for ( size_t i = 0; i < 20; i++ ){
    ComplexVectorType Vec = Vecs.col(i);
    ComplexVectorType Npi = NPhonon( Bases, Vec, Ham0);
    ComplexVectorType Nfi = NFermion( Bases, Vec, Ham0);
    LogOut << i << ": E = " << Vals[i] << ", Np = " << arma::accu(Npi) << ", Nf = " << Nfi.st() << std::endl;
  }
  LogOut.close();
}

/* main program */
int main(int argc, char *argv[]){
  if ( argc < 2 ) RUNTIME_ERROR(" Need at least one argument to run program. Use 0 to run Equilibrium.");
#ifdef MKL
  mkl_set_num_threads(NumCores);
#endif
  Equilibrium("");
  return 0;
}
