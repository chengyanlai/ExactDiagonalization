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
  const int OBC = 0;
  int N = 20;
  int L = 2 * N;
  int dynamics = 0;
  int Tsteps = 3000;
  RealType dt = 0.005;

  std::vector<RealType> Momentum;
  Momentum.push_back(0.0);
  // for ( int i = 0; i < 100; i ++ ) Momentum.push_back( 0.01 * RealType(i) );

  std::vector<RealType> Jin, Win, Gin;
  Jin = std::vector<RealType>(L, 1.0);// t0 = 1
  // Win = std::vector<RealType>(L,10.00);// omega_0 = 10
  // Gin = std::vector<RealType>(L,20.00);// g = 20; \lambda = g^2 / (2 t0 \omega) = 20
  Win = std::vector<RealType>(L, 1.0);// testing
  Gin = std::vector<RealType>(L, 1.0);// testing - 1
  // Gin = std::vector<RealType>(L, sqrt(2.0));// testing - 1

  // LogOut << "Build Lattice - " << std::endl;
  // std::vector<DT> JWork(Jin.begin(), Jin.end());
  // const std::vector< Node<DT>* > lattice = NN_1D_Chain(L, JWork, OBC);
  // for ( auto &lt : lattice ){
  //   if ( !(lt->VerifySite(LogOut)) ) RUNTIME_ERROR("Wrong lattice setup!");
  // }
  // LogOut << "DONE!" << std::endl;
  LogOut << "Build Basis - " << std::flush;
  Basis B1(L, N);
  B1.Phonon();
  LogOut << B1.GetHilbertSpace() << std::flush;
  std::string filename = "LFS-L";
  filename.append( std::to_string( (unsigned long)L) );
  filename.append( "N" );
  filename.append( std::to_string( (unsigned long)N) );
  filename.append( ".h5" );
  B1.Save(prefix + filename, "LFS" );
  std::vector<Basis> Bases;
  Bases.push_back(B1);
  LogOut << " DONE!" << std::endl;
  LogOut << "Build Hamiltonian - " << std::flush;
  Holstein<DT> Ham0( Bases );
  std::vector<ComplexType> Wloc(Win.begin(), Win.end());
  std::vector<ComplexType> Gloc(Gin.begin(), Gin.end());

  for( int i = 0; i < Momentum.size(); i++ ){
    // Ham0.HolsteinModel( Bases, Momentum, lattice,  Wloc,  Gloc );
    Ham0.HolsteinModel( Bases, Momentum.at(i), Jin.at(0),  Wloc,  Gloc );
    // ComplexMatrixType H(Ham0.GetTotalHamiltonian());
    // std::cout << H << std::endl;
    LogOut << "Hermitian = " << Ham0.CheckHermitian() << ", Hilbert space = " << Ham0.GetTotalHilbertSpace() << ", DONE!" << std::endl;
    LogOut << "Diagonalize Hamiltonian - " << std::flush;
    RealVectorType Vals;
    ComplexMatrixType Vecs;
    if ( Ham0.GetTotalHilbertSpace() > 1200 ){
      Ham0.eigh(Vals, Vecs, 2);
    }else{
      Ham0.diag(Vals, Vecs);// Full spectrum
      std::ofstream SOut;
      SOut.open(prefix + "Holstein.1d.eqm.spectrum.L" + std::to_string((unsigned long)L) + "N" + std::to_string( (unsigned long)N), std::ios::app);
      SOut << std::setprecision(5) << Momentum.at(i) << "\t";
      for ( int j = 0; j < 30 ; j++){
        SOut << std::setprecision(12) << Vals[j] << "\t";
      }
      SOut<< std::endl;
      SOut.close();
    }
    LogOut << "DONE!" << std::endl;
    LogOut << "\tGS energy = " << std::setprecision(12) << Vals[0] << std::endl;
    LogOut << "\tFES energy = " << std::setprecision(12) << Vals[1] << std::endl;
  }

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
