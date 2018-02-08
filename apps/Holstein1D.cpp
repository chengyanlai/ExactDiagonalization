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

#define DT ComplexType
#define DTV ComplexVectorType
#define DTM ComplexMatrixType
// #define DT RealType
// #define DTV RealVectorType
// #define DTM RealMatrixType

void Equilibrium(const std::string prefix){
  std::ofstream LogOut;
  LogOut.open(prefix + "Holstein.1d.eqm", std::ios::app);
  int L = 16;
  const int OBC = 0;
  int N = 18;
  RealType Momentum = 0.0;
  int dynamics = 0;
  int Tsteps = 3000;
  RealType dt = 0.005;
  std::vector<RealType> Jin, Win, Gin;
  Jin = std::vector<RealType>(L, 1.0);// t0 = 1
  // Win = std::vector<RealType>(L,10.00);// omega_0 = 10
  // Gin = std::vector<RealType>(L,20.00);// g = 20; \lambda = g^2 / (2 t0 \omega) = 20
  Win = std::vector<RealType>(L, 1.0);// testing
  Gin = std::vector<RealType>(L, 1.0);// testing

  LogOut << "Build Lattice - " << std::endl;
  std::vector<DT> JWork(Jin.begin(), Jin.end());
  const std::vector< Node<DT>* > lattice = NN_1D_Chain(L, JWork, OBC);
  for ( auto &lt : lattice ){
    if ( !(lt->VerifySite(LogOut)) ) RUNTIME_ERROR("Wrong lattice setup!");
  }
  LogOut << "DONE!" << std::endl;
  LogOut << "Build Basis - " << std::flush;
  Basis B1(L, N);
  B1.Phonon();
  LogOut << B1.GetHilbertSpace() << std::flush;
  std::string filename = "LFS-L";
  filename.append( std::to_string( (unsigned long)L) );
  filename.append( ".h5" );
  B1.Save(prefix + filename, std::to_string( (unsigned long)N) );
  std::vector<Basis> Bases;
  Bases.push_back(B1);
  LogOut << " DONE!" << std::endl;
  LogOut << B1 << std::endl;
  LogOut << "Build Hamiltonian - " << std::flush;
  Holstein<DT> Ham0( Bases );
  std::vector<ComplexType> Wloc(Win.begin(), Win.end());
  std::vector<ComplexType> Gloc(Gin.begin(), Gin.end());
  Ham0.HolsteinModel( Bases, Momentum, lattice,  Wloc,  Gloc );
  // ComplexMatrixType H(Ham0.GetTotalHamiltonian());
  // std::cout << H << std::endl;
  LogOut << "Hermitian = " << Ham0.CheckHermitian() << ", Hilbert space = " << Ham0.GetTotalHilbertSpace() << ", DONE!" << std::endl;
  LogOut << "Diagonalize Hamiltonian - " << std::flush;
  RealVectorType Vals;
  ComplexMatrixType Vecs;
  // Ham0.eigh(Vals, Vecs, 2);
  Ham0.diag(Vals, Vecs);// Full spectrum
  LogOut << "DONE!" << std::endl;
  LogOut << "\tGS energy = " << Vals[0] << std::endl;
  LogOut << "\tFES energy = " << Vals[1] << std::endl;
  // HDF5IO* file = new HDF5IO(prefix + "Holstein.h5");
  // delete file;
  LogOut.close();
}

/* main program */
int main(int argc, char *argv[]){
  if ( argc < 2 ) RUNTIME_ERROR(" Need at least one argument to run program. Use 0 to run Equilibrium.");
#if defined(MKL)
  mkl_set_num_threads(NumCores);
#endif
  Equilibrium("");
  return 0;
}
