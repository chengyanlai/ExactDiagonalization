#include <iostream>
#include <sstream>
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
#ifdef MKL
  #include "mkl.h"
#endif

#ifndef NumCores
#define NumCores 2
#endif

void SaveObs( const std::string filename, const std::string gname,
  const std::vector<ComplexType> &v, const ComplexMatrixType &m);
void LoadEqmParameters( const std::string filename, int &L, int &OBC,
  int &N, RealType &Uloc, std::vector<RealType> &Vloc);
void LoadDynParameters( const std::string filename, RealType &dt, int &Tsteps,
  RealType &Uloc, std::vector<RealType> &Vloc);
std::vector<ComplexType> Ni( const std::vector<Basis> &Bases, const ComplexVectorType &Vec );
ComplexMatrixType NiNj( const std::vector<Basis> &Bases, const ComplexVectorType &Vec );

int main(int argc, char const *argv[]) {
  Eigen::setNbThreads(NumCores);
#ifdef MKL
  mkl_set_num_threads(NumCores);
#endif
  INFO("Eigen3 uses " << Eigen::nbThreads() << " threads.");
  int L;
  int OBC;
  int N;
  RealType Uin;
  std::vector<RealType> Vin;
  LoadEqmParameters( "Eqm.h5", L, OBC, N, Uin, Vin );
  INFO("Build Lattice - ");
  std::vector<ComplexType> J;
  if ( OBC ){
    J = std::vector<ComplexType>(L - 1, ComplexType(1.0, 0.0));
  } else{
    J = std::vector<ComplexType>(L, ComplexType(1.0, 0.0));
  }
  for ( auto &val : J ){
    INFO_NONEWLINE(val << " ");
  }
  INFO("");
  const std::vector< Node<ComplexType, int>* > lattice = NN_1D_Chain(L, J, OBC);
  for ( auto &lt : lattice ){
    if ( !(lt->VerifySite()) ) RUNTIME_ERROR("Wrong lattice setup!");
  }
  INFO("DONE!");
  INFO("Build Basis - ");
  Basis B1(L, N);
  B1.Boson();
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
  HDF5IO file("GS.h5");
  file.saveNumber("1DChain", "L", L);
  file.saveStdVector("1DChain", "J", J);
  file.saveNumber("Basis", "N", N);
  file.saveVector("GS", "EVec", Vec);
  file.saveNumber("GS", "EVal", Val);
  INFO("DONE!");
  std::vector<ComplexType> Nbi = Ni( Bases, Vec );
  INFO("N_i = ");
  for (auto &n : Nbi ){
    INFO_NONEWLINE( n.real() << " " );
  }
  INFO(" ");
  ComplexMatrixType Nij = NiNj( Bases, Vec );
  INFO("N_iN_j = ");
  INFO(Nij.real());
  INFO("N_i^2 = ");
  INFO(Nij.diagonal().real());
  SaveObs("SourceSink.h5", "Obs", Nbi, Nij);
  /* NOTE: Real-time dynamics */
  int Tsteps;
  RealType dt;
  Vin.clear();
  LoadDynParameters( "Dyn.h5", dt, Tsteps, Uin, Vin);
  Vtmp.clear();
  INFO("Quench the local potetnial to");
  for ( RealType &val : Vin ){
    Vtmp.push_back((ComplexType)val);
    INFO_NONEWLINE(" " << val);
  }
  Vloc.clear();
  Vloc.push_back(Vtmp);
  INFO("");
  INFO("Quench the global interaction to " << Uin);
  Utmp.assign(L, (ComplexType)Uin);
  Uloc.clear();
  Uloc.push_back(Utmp);
  INFO("Update the local/total Hamiltonian");
  ham.BuildLocalHamiltonian(Vloc, Uloc, Bases);
  ham.BuildTotalHamiltonian();
  INFO("Start Time-Evolution");
  for (size_t cntT = 1; cntT <= Tsteps; cntT++) {
    ComplexType Prefactor = ComplexType(0.0, dt);/* NOTE: hbar = 1 */
    ham.expH(Prefactor, Vec);
    /* NOTE: Expectation values */
    Nbi = Ni( Bases, Vec );
    Nij = NiNj( Bases, Vec );
    /* NOTE: H5 group name */
    std::string gname = "Obs-";
    gname.append(std::to_string((unsigned long long)cntT));
    gname.append("/");
    SaveObs("SourceSink.h5", gname, Nbi, Nij);
  }
  INFO("End Program!!");
  return 0;
}

void SaveObs( const std::string filename, const std::string gname,
  const std::vector<ComplexType> &v, const ComplexMatrixType &m){
  HDF5IO file(filename);
  file.saveStdVector(gname, "Nb", v);
  ComplexType Ntot = std::accumulate(v.begin(), v.end(), ComplexType(0.0, 0.0));
  file.saveNumber(gname, "Ntot", Ntot.real());
  file.saveMatrix(gname, "Nij", m);
}

void LoadEqmParameters( const std::string filename, int &L, int &OBC, int &N,
  RealType &Uloc, std::vector<RealType> &Vloc){
    HDF5IO file(filename);
    L = file.loadInt("Parameters", "L");
    OBC = file.loadInt("Parameters", "OBC");
    N = file.loadInt("Parameters", "N");
    Uloc = file.loadReal("Parameters", "U");
    file.loadStdVector("Parameters", "V", Vloc);
}

void LoadDynParameters( const std::string filename, RealType &dt, int &Tsteps,
  RealType &Uloc, std::vector<RealType> &Vloc){
    HDF5IO file(filename);
    dt = file.loadReal("Parameters", "dt");
    Tsteps = file.loadInt("Parameters", "Tsteps");
    Uloc = file.loadReal("Parameters", "U");
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
