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

DTV Ni( const std::vector<Basis> &Bases, const DTV &Vec ){
  DTV tmp(Bases.at(0).getL(), arma::fill::zeros);
  std::vector< std::vector<int> > b = Bases.at(0).getBStates();
  assert( b.size() == Vec.size() );
  int coff = 0;
  for ( auto &nbi : b ){
    for (size_t cnt = 0; cnt < Bases.at(0).getL(); cnt++) {
      tmp.at(cnt) += (RealType)nbi.at(cnt) * std::pow(std::abs(Vec(coff)), 2);
    }
    coff++;
  }
  return tmp;
}

DTM NiNj( const std::vector<Basis> &Bases, const DTV &Vec ){
  DTM tmp(Bases.at(0).getL(), Bases.at(0).getL(), arma::fill::zeros);
  std::vector< std::vector<int> > b = Bases.at(0).getBStates();
  assert( b.size() == Vec.size() );
  int coff = 0;
  for ( auto &nbi : b ){
    for (size_t cnt1 = 0; cnt1 < Bases.at(0).getL(); cnt1++) {
      for (size_t cnt2 = 0; cnt2 < Bases.at(0).getL(); cnt2++) {
        tmp(cnt1, cnt2) += (RealType)nbi.at(cnt1) * (RealType)nbi.at(cnt2) * std::pow(std::abs(Vec(coff)), 2);
      }
    }
    coff++;
  }
  return tmp;
}

void Equilibrium(const std::string prefix){
  std::ofstream LogOut;
  LogOut.open(prefix + "bhm.1d.eqm", std::ios::app);
  int L = 10;
  int OBC = 1;
  int N = 10;
  int dynamics = 0;
  int Tsteps = 3000;
  RealType dt=0.005;
  std::vector<RealType> Jin, Uin, Vin;
  Jin = std::vector<RealType>(L-1, 1.0);// OBC
  Uin = std::vector<RealType>(L, 1.0);
  Vin = std::vector<RealType>(L, 0.0);

  HDF5IO* file = new HDF5IO("BHMChainData.h5");
  LogOut << "Build Lattice - " << std::endl;
  std::vector<DT> JWork(Jin.begin(), Jin.end());
  const std::vector< Node<DT>* > lattice = NN_1D_Chain(L, JWork, OBC);
  file->SaveNumber("1DChain", "L", L);
  file->SaveNumber("1DChain", "N", N);
  file->SaveStdVector("1DChain", "U", Uin);
  file->SaveStdVector("1DChain", "V", Vin);
  file->SaveStdVector("1DChain", "J", Jin);
  for ( auto &lt : lattice ){
    if ( !(lt->VerifySite(LogOut)) ) RUNTIME_ERROR("Wrong lattice setup!");
  }
  LogOut << "DONE!" << std::endl;
  LogOut << "Build Basis - " << std::flush;
  Basis B1(L, N);
  B1.Boson();
  std::vector<Basis> Bases;
  Bases.push_back(B1);
  LogOut << "DONE!" << std::endl;
  LogOut << "Build Hamiltonian - " << std::flush;
  Hamiltonian<DT> Ham0( Bases );
  std::vector<DT> UWork(Uin.begin(), Uin.end());
  std::vector<DT> VWork(Vin.begin(), Vin.end());
  Ham0.BoseHubbardModel(Bases, lattice, VWork, UWork);
  Ham0.CheckHermitian();
  LogOut << Ham0.GetTotalHilbertSpace() << " DONE!" << std::endl;
  LogOut << "Diagonalize Hamiltonian - " << std::flush;
  RealVectorType Vals;
  DTM Vecs;
  Ham0.eigh(Vals, Vecs, 2);
  // Ham0.diag(Vals, Vecs);// Full spectrum
  LogOut << "DONE!" << std::endl;
  LogOut << "\tGS energy = " << Vals[0] << std::endl;
  LogOut << "\tFES energy = " << Vals[1] << std::endl;
  DTV Vec = Vecs.col(0);
// RealSparseMatrixType H0 = Ham0.getTotalHamiltonian();
// std::cout << arma::cdot(Vec, H0 * Vec) << std::endl;
  file->SaveVector("GS", "EVec", Vec);
  file->SaveVector("GS", "EVal", Vals);
  DTV Nbi = Ni( Bases, Vec );
  DT NbT = arma::sum(Nbi);
  LogOut << " Nbi - " << std::endl;
  LogOut << Nbi << std::endl;
  LogOut << "Total Nb = " << NbT << std::endl;
  DTM Nij = NiNj( Bases, Vec );
  LogOut << " NiNj - " << std::endl;
  LogOut << Nij << std::endl;
  file->SaveVector("Obs", "Nb", Nbi);
  file->SaveNumber("Obs", "NbT", NbT);
  file->SaveMatrix("Obs", "Nij", Nij);
  delete file;
}

/* main program */
int main(int argc, char *argv[]){
  if ( argc < 2 ) RUNTIME_ERROR(" Need at least one argument to run program. Use 0 to run Equilibrium.");
#if defined(MKL)
  mkl_set_dynamic(0);
  mkl_set_num_threads(NumCores);
#endif
  Equilibrium("");
  return 0;
}
