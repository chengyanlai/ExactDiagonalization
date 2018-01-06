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

#ifdef MKL
  #include "mkl.h"
#endif

#ifndef DEBUG
#define DEBUG 1
#endif

#ifndef NumCores
#define NumCores 12
#endif

// #define DTYPE 0//comment out this to use complex

// #ifndef DTYPE
// #define DT ComplexType
// #define DTV ComplexVectorType
// #define DTM ComplexMatrixType
// #else
#define DT RealType
#define DTV RealVectorType
#define DTM RealMatrixType
// #endif

void LoadParameters( const std::string filename, int &L, int &OBC, int &N1, int &N2,
  std::vector<RealType> &Jls, std::vector<RealType> &Uls, std::vector<RealType> &Vls,
  int &dynamics, int &Tsteps, RealType &dt){
    HDF5IO file(filename);
    L = file.loadInt("Parameters", "L");
    OBC = file.loadInt("Parameters", "OBC");
    N1 = file.loadInt("Parameters", "N1");
    N2 = file.loadInt("Parameters", "N2");
    file.loadStdVector("Parameters", "J", Jls);
    file.loadStdVector("Parameters", "U", Uls);
    file.loadStdVector("Parameters", "V", Vls);
    dynamics = file.loadInt("Parameters", "dynamics");
    Tsteps = file.loadInt("Parameters", "Tsteps");
    dt = file.loadReal("Parameters", "dt");
}

std::vector< DTV > Ni( const std::vector<Basis> &Bases, const DTV &Vec,
  Hamiltonian<DT> &ham );
DTM NupNdn( const std::vector<Basis> &Bases, const DTV &Vec,
  Hamiltonian<DT> &ham );
DTM NupNup( const std::vector<Basis> &Bases, const DTV &Vec,
  Hamiltonian<DT> &ham );
DTM NdnNdn( const std::vector<Basis> &Bases, const DTV &Vec,
  Hamiltonian<DT> &ham );

int main(int argc, char const *argv[]) {
  Eigen::setNbThreads(NumCores);
#ifdef MKL
  mkl_set_num_threads(NumCores);
#endif
  INFO("Eigen3 uses " << Eigen::nbThreads() << " threads.");
  int L;
  int OBC;
  int N1, N2;
  int dynamics, Tsteps;
  RealType dt;
  std::vector<RealType> Jin, Uin, Vin;
  // Load parameters from file
  try{
    H5::Exception::dontPrint();
    H5::H5File::isHdf5("conf.h5");
    LoadParameters( "conf.h5", L, OBC, N1, N2, Jin, Uin, Vin, dynamics, Tsteps, dt);
    std::cout << "Parameters for calculation loaded!" << std::endl;
  }catch(H5::FileIException){
    L = 14;
    OBC = 1;
    N1 = 6;
    N2 = 6;
    Jin = std::vector<RealType>(L-1, 1.0);
    Uin = std::vector<RealType>(L, 9.0);
    Vin = std::vector<RealType>(L, 0.0);
    dynamics = 0;
    Tsteps = 0;
    dt = 0.005;
  }
  HDF5IO *file = new HDF5IO("FData.h5");
  INFO("Build Lattice - ");
  std::vector<DT> J;
  for ( size_t i = 0; i < Jin.size(); i++ ){
    J.push_back(DT( Jin.at(i) ));
  }
  for ( auto &val : J ){
    INFO_NONEWLINE(val << " ");
  }
  INFO("");
  const std::vector< Node<DT>* > lattice = NN_1D_Chain(L, J, OBC);
  file->saveNumber("1DChain", "L", L);
  file->saveStdVector("1DChain", "J", J);
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
  file->saveNumber("Basis", "N1", N1);
  file->saveStdVector("Basis", "F1States", st1);
  file->saveStdVector("Basis", "F1Tags", tg1);
  file->saveNumber("Basis", "N2", N2);
  file->saveStdVector("Basis", "F2States", st2);
  file->saveStdVector("Basis", "F2Tags", tg2);
  INFO("DONE!");
  INFO_NONEWLINE("Build Hamiltonian - ");
  std::vector<Basis> Bases;
  Bases.push_back(F1);
  Bases.push_back(F2);
  Hamiltonian<DT> ham( Bases );
  std::vector< std::vector<DT> > Vloc;
  std::vector<DT> Vtmp;//(L, 1.0);
  for ( RealType &val : Vin ){
    Vtmp.push_back((DT)val);
  }
  Vloc.push_back(Vtmp);
  Vloc.push_back(Vtmp);
  std::vector< std::vector<DT> > Uloc;
  // std::vector<DT> Utmp(L, DT(10.0e0, 0.0e0) );
  std::vector<DT> Utmp;
  for ( RealType &val : Uin ){
    Utmp.push_back((DT)val);
  }
  Uloc.push_back(Utmp);
  Uloc.push_back(Utmp);
  ham.BuildLocalHamiltonian(Vloc, Uloc, Bases);
  INFO(" - BuildLocalHamiltonian DONE!");
  ham.BuildHoppingHamiltonian(Bases, lattice);
  INFO(" - BuildHoppingHamiltonian DONE!");
  ham.BuildTotalHamiltonian();
  INFO("DONE!");
  INFO_NONEWLINE("Diagonalize Hamiltonian - ");
  RealVectorType Vals;
  RealMatrixType Vecs;
  ham.eigh(Vals, Vecs, 2);
  INFO("GS energy = " << Vals[0]);
  // file->saveVector("GS", "EVec", Vec);
  // file->saveStdVector("GS", "EVal", Val);
  // INFO("DONE!");
  // std::vector< DTV > Nfi = Ni( Bases, Vec, ham );
  // INFO(" Up Spin - ");
  // INFO(Nfi.at(0));
  // INFO(" Down Spin - ");
  // INFO(Nfi.at(1));
  // INFO(" N_i - ");
  // DTV Niall = Nfi.at(0) + Nfi.at(1);
  // INFO(Niall);
  // DTM Nud = NupNdn( Bases, Vec, ham );
  // INFO(" Correlation NupNdn");
  // INFO(Nud);
  // DTM Nuu = NupNup( Bases, Vec, ham );
  // INFO(" Correlation NupNup");
  // INFO(Nuu);
  // DTM Ndd = NdnNdn( Bases, Vec, ham );
  // INFO(" Correlation NdnNdn");
  // INFO(Ndd);
  // INFO(" N_i^2 - ");
  // DTM Ni2 = Nuu.diagonal() + Ndd.diagonal() + 2.0e0 * Nud.diagonal();
  // INFO(Ni2);
  // file->saveVector("Obs", "Nup", Nfi.at(0));
  // file->saveVector("Obs", "Ndn", Nfi.at(1));
  // file->saveMatrix("Obs", "NupNdn", Nud);
  // file->saveMatrix("Obs", "NupNup", Nuu);
  // file->saveMatrix("Obs", "NdnNdn", Ndd);
  // delete file;
  // if ( dynamics ){
  //   ComplexType Prefactor = ComplexType(0.0, -1.0e0*dt);/* NOTE: hbar = 1 */
  //   std::cout << "Begin dynamics......" << std::endl;
  //   std::cout << "Cut the boundary." << std::endl;
  //   J.pop_back();
  //   std::vector< Node<DT>* > lattice2 = NN_1D_Chain(L, J, true);// cut to open
  //   ham.BuildHoppingHamiltonian(Bases, lattice2);
  //   INFO(" - Update Hopping Hamiltonian DONE!");
  //   ham.BuildTotalHamiltonian();
  //   INFO(" - Update Total Hamiltonian DONE!");
  //   for (size_t cntT = 1; cntT <= Tsteps; cntT++) {
  //     ham.expH(Prefactor, Vec);
  //     if ( cntT % 2 == 0 ){
  //       HDF5IO file2("DYN.h5");
  //       std::string gname = "Obs-";
  //       gname.append(std::to_string((unsigned long long)cntT));
  //       gname.append("/");
  //       Nfi = Ni( Bases, Vec, ham );
  //       Nud = NupNdn( Bases, Vec, ham );
  //       Nuu = NupNup( Bases, Vec, ham );
  //       Ndd = NdnNdn( Bases, Vec, ham );
  //       file2.saveVector(gname, "Nup", Nfi.at(0));
  //       file2.saveVector(gname, "Ndn", Nfi.at(1));
  //       file2.saveMatrix(gname, "NupNdn", Nud);
  //       file2.saveMatrix(gname, "NupNup", Nuu);
  //       file2.saveMatrix(gname, "NdnNdn", Ndd);
  //     }
  //   }
  // }
  return 0;
}

std::vector< DTV > Ni( const std::vector<Basis> &Bases,
  const DTV &Vec, Hamiltonian<DT> &ham ){
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
  const DTV &Vec, Hamiltonian<DT> &ham ){
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
  const DTV &Vec, Hamiltonian<DT> &ham ){
  DTM out = DTM::Zero(Bases.at(0).getL(), Bases.at(0).getL() );
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
  const DTV &Vec, Hamiltonian<DT> &ham ){
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
