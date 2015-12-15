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

// #define DTYPE 0//comment out this to compile complex

#ifndef DTYPE
#define DT ComplexType
#define DTV ComplexVectorType
#define DTM ComplexMatrixType
#else
#define DT RealType
#define DTV RealVectorType
#define DTM RealMatrixType
#endif

void LoadParameters( const std::string filename, int &L, RealType &J12ratio,
  int &OBC, int &N1, int &N2, RealType &Uloc,
  std::vector<RealType> &Vloc, RealType &phi);
std::vector< DTV > Ni( const std::vector<Basis> &Bases, const DTV &Vec,
  Hamiltonian<DT,int> &ham );
DTM NupNdn( const std::vector<Basis> &Bases, const DTV &Vec,
  Hamiltonian<DT,int> &ham );
DTM NupNup( const std::vector<Basis> &Bases, const DTV &Vec,
  Hamiltonian<DT,int> &ham );
DTM NdnNdn( const std::vector<Basis> &Bases, const DTV &Vec,
  Hamiltonian<DT,int> &ham );

int main(int argc, char const *argv[]) {
  int L;
  RealType J12ratio;
  int OBC;
  int N1;
  int N2;
  RealType Uin, phi;
  std::vector<RealType> Vin;
  LoadParameters( "conf.h5", L, J12ratio, OBC, N1, N2, Uin, Vin, phi);
  HDF5IO file("FSSH.h5");
  // const int L = 5;
  // const bool OBC = true;
  // const RealType J12ratio = 0.010e0;
  INFO("Build Lattice - ");
  std::vector<DT> J;
  if ( OBC ){
    // J = std::vector<DT>(L - 1, 0.0);// NOTE: Atomic limit testing
    J = std::vector<DT>(L - 1, 1.0);
    if ( J12ratio > 1.0e0 ) {
      for (size_t cnt = 1; cnt < L-1; cnt+=2) {
        J.at(cnt) /= J12ratio;
      }
    } else{
      for (size_t cnt = 0; cnt < L-1; cnt+=2) {
        J.at(cnt) *= J12ratio;
      }
    }
  } else{
    J = std::vector<DT>(L, 1.0);
    if ( J12ratio > 1.0e0 ) {
      for (size_t cnt = 1; cnt < L; cnt+=2) {
        J.at(cnt) /= J12ratio;
      }
    } else{
      for (size_t cnt = 0; cnt < L; cnt+=2) {
        J.at(cnt) *= J12ratio;
      }
    }
#ifndef DTYPE
    if ( std::abs(phi) > 1.0e-10 ){
      J.at(L-1) *= exp( DT(0.0, 1.0) * phi );
    }
#endif
  }
  for ( auto &val : J ){
    INFO_NONEWLINE(val << " ");
  }
  INFO("");
  const std::vector< Node<DT, int>* > lattice = NN_1D_Chain(L, J, OBC);
  file.saveNumber("1DChain", "L", L);
  file.saveStdVector("1DChain", "J", J);
  for ( auto &lt : lattice ){
    if ( !(lt->VerifySite()) ) RUNTIME_ERROR("Wrong lattice setup!");
  }
  INFO("DONE!");
  INFO("Build Basis - ");
  // int N1 = (L+1)/2;
  Basis F1(L, N1, true);
  F1.Fermion();
  std::vector<int> st1 = F1.getFStates();
  std::vector<size_t> tg1 = F1.getFTags();
  // for (size_t cnt = 0; cnt < st1.size(); cnt++) {
  //   INFO_NONEWLINE( std::setw(3) << st1.at(cnt) << " - ");
  //   F1.printFermionBasis(st1.at(cnt));
  //   INFO("- " << tg1.at(st1.at(cnt)));
  // }
  // int N2 = (L-1)/2;
  Basis F2(L, N2, true);
  F2.Fermion();
  std::vector<int> st2 = F2.getFStates();
  std::vector<size_t> tg2 = F2.getFTags();
  // for (size_t cnt = 0; cnt < st2.size(); cnt++) {
  //   INFO_NONEWLINE( std::setw(3) << st2.at(cnt) << " - ");
  //   F2.printFermionBasis(st2.at(cnt));
  //   INFO("- " << tg2.at(st2.at(cnt)));
  // }
  file.saveNumber("Basis", "N1", N1);
  file.saveStdVector("Basis", "F1States", st1);
  file.saveStdVector("Basis", "F1Tags", tg1);
  file.saveNumber("Basis", "N2", N2);
  file.saveStdVector("Basis", "F2States", st2);
  file.saveStdVector("Basis", "F2Tags", tg2);
  INFO("DONE!");
  INFO_NONEWLINE("Build Hamiltonian - ");
  std::vector<Basis> Bases;
  Bases.push_back(F1);
  Bases.push_back(F2);
  Hamiltonian<DT,int> ham( Bases );
  std::vector< std::vector<DT> > Vloc;
  std::vector<DT> Vtmp;//(L, 1.0);
  for ( RealType &val : Vin ){
    Vtmp.push_back((DT)val);
  }
  Vloc.push_back(Vtmp);
  Vloc.push_back(Vtmp);
  std::vector< std::vector<DT> > Uloc;
  // std::vector<DT> Utmp(L, DT(10.0e0, 0.0e0) );
  std::vector<DT> Utmp(L, (DT)Uin);
  Uloc.push_back(Utmp);
  Uloc.push_back(Utmp);
  ham.BuildLocalHamiltonian(Vloc, Uloc, Bases);
  ham.BuildHoppingHamiltonian(Bases, lattice);
  ham.BuildTotalHamiltonian();
  INFO("DONE!");
  INFO_NONEWLINE("Diagonalize Hamiltonian - ");
  std::vector<RealType> Val;
  Hamiltonian<DT,int>::VectorType Vec;
  ham.eigh(Val, Vec);
  INFO("GS energy = " << Val.at(0));
  file.saveVector("GS", "EVec", Vec);
  file.saveStdVector("GS", "EVal", Val);
  INFO("DONE!");
  std::vector< DTV > Nfi = Ni( Bases, Vec, ham );
  INFO(" Up Spin - ");
  INFO(Nfi.at(0));
  INFO(" Down Spin - ");
  INFO(Nfi.at(1));
  INFO(" N_i - ");
  DTV Ni = Nfi.at(0) + Nfi.at(1);
  INFO(Ni);
  DTM Nud = NupNdn( Bases, Vec, ham );
  INFO(" Correlation NupNdn");
  INFO(Nud);
  DTM Nuu = NupNup( Bases, Vec, ham );
  INFO(" Correlation NupNup");
  INFO(Nuu);
  DTM Ndd = NdnNdn( Bases, Vec, ham );
  INFO(" Correlation NdnNdn");
  INFO(Ndd);
  INFO(" N_i^2 - ");
  DTM Ni2 = Nuu.diagonal() + Ndd.diagonal() + 2.0e0 * Nud.diagonal();
  INFO(Ni2);
  file.saveVector("Obs", "Nup", Nfi.at(0));
  file.saveVector("Obs", "Ndn", Nfi.at(1));
  file.saveMatrix("Obs", "NupNdn", Nud);
  file.saveMatrix("Obs", "NupNup", Nuu);
  file.saveMatrix("Obs", "NdnNdn", Ndd);
  return 0;
}

void LoadParameters( const std::string filename, int &L, RealType &J12ratio,
  int &OBC, int &N1, int &N2, RealType &Uloc,
  std::vector<RealType> &Vloc, RealType &phi){
    HDF5IO file(filename);
    L = file.loadInt("Parameters", "L");
    J12ratio = file.loadReal("Parameters", "J12");
    OBC = file.loadInt("Parameters", "OBC");
    N1 = file.loadInt("Parameters", "N1");
    N2 = file.loadInt("Parameters", "N2");
    Uloc = file.loadReal("Parameters", "U");
    phi = file.loadReal("Parameters", "phi");
    file.loadStdVector("Parameters", "V", Vloc);
}

std::vector< DTV > Ni( const std::vector<Basis> &Bases,
  const DTV &Vec, Hamiltonian<DT,int> &ham ){
  std::vector< DTV > out;
  DTV tmp1 = DTV::Zero(Bases.at(0).getL());//(Bases.at(0).getL(), 0.0e0);
  DTV tmp2 = DTV::Zero(Bases.at(1).getL());//(Bases.at(1).getL(), 0.0e0);
  std::vector< int > f1 = Bases.at(0).getFStates();
  std::vector< int > f2 = Bases.at(1).getFStates();
  size_t f1id = 0, f2id = 0;
  for ( int &nf2 : f2 ){
    std::vector<size_t> ids(2,f2id);
    f1id = 0;
    for ( int &nf1 : f1 ){
      ids.at(0) = f1id;
      size_t id = ham.DetermineTotalIndex(ids);
      for (size_t cnt = 0; cnt < Bases.at(0).getL(); cnt++) {
        if ( btest(nf1, cnt) ) tmp1(cnt) += std::pow(std::abs(Vec(id)), 2);//Vec(id) * std::conj( Vec(id) );
      }
      for (size_t cnt = 0; cnt < Bases.at(1).getL(); cnt++) {
        if ( btest(nf2, cnt) ) tmp2(cnt) += std::pow(std::abs(Vec(id)), 2);//Vec(id) * std::conj( Vec(id) );
      }
      f1id++;
    }
    f2id++;
  }
  out.push_back(tmp1);
  out.push_back(tmp2);
  return out;
}

DTM NupNdn( const std::vector<Basis> &Bases,
  const DTV &Vec, Hamiltonian<DT,int> &ham ){
  DTM out = DTM::Zero(
    Bases.at(0).getL(), Bases.at(1).getL() );
  std::vector< int > f1 = Bases.at(0).getFStates();
  std::vector< int > f2 = Bases.at(1).getFStates();
  size_t f1id = 0, f2id = 0;
  for ( int &nf2 : f2 ){
    std::vector<size_t> ids(2,f2id);
    f1id = 0;
    for ( int &nf1 : f1 ){
      ids.at(0) = f1id;
      size_t id = ham.DetermineTotalIndex(ids);
      for (size_t cnt1 = 0; cnt1 < Bases.at(0).getL(); cnt1++) {
        for (size_t cnt2 = 0; cnt2 < Bases.at(1).getL(); cnt2++) {
          if ( btest(nf2, cnt1) && btest(nf1, cnt2) ) {
            out(cnt1, cnt2) +=  std::pow(std::abs(Vec(id)), 2);//Vec(id) * std::conj( Vec(id) );
          }
        }
      }
      f1id++;
    }
    f2id++;
  }
  return out;
}

DTM NupNup( const std::vector<Basis> &Bases,
  const DTV &Vec, Hamiltonian<DT,int> &ham ){
  DTM out = DTM::Zero(
    Bases.at(0).getL(), Bases.at(0).getL() );
  std::vector< int > f1 = Bases.at(0).getFStates();
  std::vector< int > f2 = Bases.at(1).getFStates();
  size_t f1id = 0, f2id = 0;
  for ( int &nf2 : f2 ){
    std::vector<size_t> ids(2,f2id);
    f1id = 0;
    for ( int &nf1 : f1 ){
      ids.at(0) = f1id;
      size_t id = ham.DetermineTotalIndex(ids);
      for (size_t cnt1 = 0; cnt1 < Bases.at(0).getL(); cnt1++) {
        for (size_t cnt2 = 0; cnt2 < Bases.at(0).getL(); cnt2++) {
          if ( btest(nf1, cnt1) && btest(nf1, cnt2) ) {
            out(cnt1, cnt2) +=  std::pow(std::abs(Vec(id)), 2);//Vec(id) * std::conj( Vec(id) );
          }
        }
      }
      f1id++;
    }
    f2id++;
  }
  return out;
}

DTM NdnNdn( const std::vector<Basis> &Bases,
  const DTV &Vec, Hamiltonian<DT,int> &ham ){
  DTM out = DTM::Zero(
    Bases.at(0).getL(), Bases.at(1).getL() );
  std::vector< int > f1 = Bases.at(0).getFStates();
  std::vector< int > f2 = Bases.at(1).getFStates();
  size_t f1id = 0, f2id = 0;
  for ( int &nf2 : f2 ){
    std::vector<size_t> ids(2,f2id);
    f1id = 0;
    for ( int &nf1 : f1 ){
      ids.at(0) = f1id;
      size_t id = ham.DetermineTotalIndex(ids);
      for (size_t cnt1 = 0; cnt1 < Bases.at(1).getL(); cnt1++) {
        for (size_t cnt2 = 0; cnt2 < Bases.at(1).getL(); cnt2++) {
          if ( btest(nf2, cnt1) && btest(nf2, cnt2) ) {
            out(cnt1, cnt2) +=  std::pow(std::abs(Vec(id)), 2);//Vec(id) * std::conj( Vec(id) );
          }
        }
      }
      f1id++;
    }
    f2id++;
  }
  return out;
}
