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
#include "src/Hamiltonian/Hamiltonian.hpp"
#include "src/hdf5io/hdf5io.hpp"

#ifdef MKL
  #include "mkl.h"
#endif

// #define DT ComplexType
// #define DTV ComplexVectorType
// #define DTM ComplexMatrixType
#define DT RealType
#define DTV RealVectorType
#define DTM RealMatrixType

void Equilibrium(const std::string prefix){
  std::ofstream LogOut;
  LogOut.open(prefix + "Holstein.1d.eqm", std::ios::app);
  int L = 16;
  const int OBC = 0;
  int N = 10;
  int dynamics = 0;
  int Tsteps = 3000;
  RealType dt=0.005;
  std::vector<RealType> Jin, Uin, Vin;
  Jin = std::vector<RealType>(L, 1.0);// PBC
  Uin = std::vector<RealType>(L, 1.0);
  Vin = std::vector<RealType>(L, 0.0);

  // HDF5IO* file = new HDF5IO("Holstein.h5");
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
  std::cout << B1.GetHilbertSpace() << std::flush;
  // B1.PrintAllBosonBasis();
  std::vector<Basis> Bases;
  Bases.push_back(B1);
  LogOut << "DONE!" << std::endl;
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
