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

void LoadParameters( const std::string filename, int &L, RealType &J12ratio, int &OBC,
  int &N, RealType &Uloc, std::vector<RealType> &Vloc, RealType &phi);
std::vector<ComplexType> Ni( const std::vector<Basis> &Bases, const ComplexVectorType &Vec );
ComplexMatrixType NiNj( const std::vector<Basis> &Bases, const ComplexVectorType &Vec );

int main(int argc, char const *argv[]) {
  int L;
  RealType J12ratio;
  int OBC;
  int N;
  RealType Uin, phi;
  std::vector<RealType> Vin;
  LoadParameters( "conf.h5", L, J12ratio, OBC, N, Uin, Vin, phi);
  HDF5IO file("SSH.h5");
  // const int L = 5;
  // const bool OBC = true;
  // const RealType J12ratio = 0.010e0;
  INFO("Build Lattice - ");
  std::vector<ComplexType> J;
  if ( OBC ){
    J = std::vector<ComplexType>(L - 1, ComplexType(1.0, 0.0));
    for (size_t cnt = 0; cnt < L-1; cnt+=2) {
      J.at(cnt) *= J12ratio;
    }
  } else{
    J = std::vector<ComplexType>(L, ComplexType(1.0, 0.0));
    for (size_t cnt = 0; cnt < L; cnt+=2) {
      J.at(cnt) *= J12ratio;
    }
    if ( std::abs(phi) > 1.0e-10 ){
      J.at(L-1) *= exp( ComplexType(0.0e0, 1.0e0) * phi );
      // INFO(exp( ComplexType(0.0e0, 1.0e0) * phi ));
    }
  }
  for ( auto &val : J ){
    INFO_NONEWLINE(val << " ");
  }
  INFO("");
  const std::vector< Node<ComplexType, int>* > lattice = NN_1D_Chain(L, J, OBC);
  file.saveNumber("1DChain", "L", L);
  file.saveNumber("1DChain", "U", Uin);
  file.saveStdVector("1DChain", "J", J);
  for ( auto &lt : lattice ){
    if ( !(lt->VerifySite()) ) RUNTIME_ERROR("Wrong lattice setup!");
  }
  INFO("DONE!");
  INFO("Build Basis - ");
  // int N1 = (L+1)/2;
  Basis B1(L, N);
  B1.Boson();
  // std::vector< std::vector<int> > st = B1.getBStates();
  // std::vector< RealType > tg = B1.getBTags();
  // for (size_t cnt = 0; cnt < tg.size(); cnt++) {
  //   INFO_NONEWLINE( std::setw(3) << cnt << " - ");
  //   for (auto &j : st.at(cnt)){
  //     INFO_NONEWLINE(j << " ");
  //   }
  //   INFO("- " << tg.at(cnt));
  // }
  file.saveNumber("1DChain", "N", N);
  // file.saveStdVector("Basis", "States", st);
  // file.saveStdVector("Basis", "Tags", tg);
  INFO("DONE!");
  INFO_NONEWLINE("Build Hamiltonian - ");
  std::vector<Basis> Bases;
  Bases.push_back(B1);
  Hamiltonian<ComplexType,int> ham( Bases );
  std::vector< std::vector<ComplexType> > Vloc;
  std::vector<ComplexType> Vtmp;//(L, 1.0);
  for ( RealType &val : Vin ){
    Vtmp.push_back((ComplexType)val);
  }
  Vloc.push_back(Vtmp);
  std::vector< std::vector<ComplexType> > Uloc;
  // std::vector<ComplexType> Utmp(L, ComplexType(10.0e0, 0.0e0) );
  std::vector<ComplexType> Utmp(L, (ComplexType)Uin);
  Uloc.push_back(Utmp);
  ham.BuildLocalHamiltonian(Vloc, Uloc, Bases);
  ham.BuildHoppingHamiltonian(Bases, lattice);
  ham.BuildTotalHamiltonian();
  INFO("DONE!");
  INFO_NONEWLINE("Diagonalize Hamiltonian - ");
  RealType Val = 0.0e0;
  Hamiltonian<ComplexType,int>::VectorType Vec;
  ham.eigh(Val, Vec);
  INFO("GS energy = " << Val);
  file.saveVector("GS", "EVec", Vec);
  file.saveNumber("GS", "EVal", Val);
  INFO("DONE!");
  std::vector<ComplexType> Nbi = Ni( Bases, Vec );
  for (auto &n : Nbi ){
    INFO( n << " " );
  }
  ComplexMatrixType Nij = NiNj( Bases, Vec );
  INFO(Nij);
  INFO(Nij.diagonal());
  file.saveStdVector("Obs", "Nb", Nbi);
  file.saveMatrix("Obs", "Nij", Nij);
  return 0;
}

void LoadParameters( const std::string filename, int &L, RealType &J12ratio,
  int &OBC, int &N, RealType &Uloc,
  std::vector<RealType> &Vloc, RealType &phi){
    HDF5IO file(filename);
    L = file.loadInt("Parameters", "L");
    J12ratio = file.loadReal("Parameters", "J12");
    OBC = file.loadInt("Parameters", "OBC");
    N = file.loadInt("Parameters", "N");
    Uloc = file.loadReal("Parameters", "U");
    phi = file.loadReal("Parameters", "phi");
    file.loadStdVector("Parameters", "V", Vloc);
}

std::vector<ComplexType> Ni( const std::vector<Basis> &Bases,
  const ComplexVectorType &Vec ){
  std::vector<ComplexType> tmp(Bases.at(0).getL(), 0.0e0);
  std::vector< std::vector<int> > b = Bases.at(0).getBStates();
  assert( b.size() == Vec.size() );
  int coff = 0;
  for ( auto &nbi : b ){
    for (size_t cnt = 0; cnt < Bases.at(0).getL(); cnt++) {
      tmp.at(cnt) += (RealType)nbi.at(cnt) * Vec(coff) * std::conj(Vec(coff));
    }
    coff++;
  }
  return tmp;
}

ComplexMatrixType NiNj( const std::vector<Basis> &Bases,
  const ComplexVectorType &Vec ){
  ComplexMatrixType tmp = ComplexMatrixType::Zero(Bases.at(0).getL(), Bases.at(0).getL());
  std::vector< std::vector<int> > b = Bases.at(0).getBStates();
  assert( b.size() == Vec.size() );
  int coff = 0;
  for ( auto &nbi : b ){
    for (size_t cnt1 = 0; cnt1 < Bases.at(0).getL(); cnt1++) {
      for (size_t cnt2 = 0; cnt2 < Bases.at(0).getL(); cnt2++) {
        tmp(cnt1, cnt2) += (RealType)nbi.at(cnt1) * (RealType)nbi.at(cnt2) *
          Vec(coff) * std::conj(Vec(coff));
      }
    }
    coff++;
  }
  return tmp;
}
