#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <numeric>
#include "src/EDType.hpp"
#include "src/bitwise.h"
#include "src/Node/Node.hpp"
#include "src/Lattice/preset.hpp"
#include "src/Basis/Basis.hpp"
#include "src/Hamiltonian/Hamiltonian.hpp"
#include "src/hdf5io/hdf5io.hpp"
#include "src/RungeKutta/RungeKutta.hpp"

#ifdef MKL
  #include "mkl.h"
#endif

#ifndef DEBUG
#define DEBUG 1
#endif

#ifndef NumCores
#define NumCores 2
#endif

RealType TraceRhos(const std::vector<ComplexMatrixType> &Rhos);
void LoadParameters( const std::string filename, int &L, int &OBC, int &N,
  RealType &Uloc, RealType &Vloc, RealType &dt, int &Tsteps, int &TBloc,
  RealType &Gamma);
void SaveObs( const std::string filename, const std::string gname,
  const std::vector<ComplexType> &v, const ComplexMatrixType &m);
void GetRho( const Basis &bs, ComplexMatrixType &Mat );
std::vector<std::vector<size_t> > IndexCollapse(const size_t TBloc,
  const std::vector<Basis> &bs);
std::vector<ComplexType> Ni( const std::vector<Basis> &Bases,
  const std::vector<ComplexMatrixType> &Vec );
ComplexMatrixType NiNj( const std::vector<Basis> &Bases,
  const std::vector<ComplexMatrixType> &Rhos );

int main(int argc, char const *argv[]) {
  Eigen::setNbThreads(NumCores);
#ifdef MKL
  mkl_set_num_threads(NumCores);
#endif
  INFO("Eigen3 uses " << Eigen::nbThreads() << " threads.");
  int L, OBC, N, Tsteps, TBloc;
  RealType Uin, Vin, dt, gamma;
  LoadParameters( "conf.h5", L, OBC, N, Uin, Vin, dt, Tsteps, TBloc, gamma );
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
  const std::vector< Node<ComplexType>* > lattice = NN_1D_Chain(L, J, OBC);
  for ( auto &lt : lattice ){
    if ( !(lt->VerifySite()) ) RUNTIME_ERROR("Wrong lattice setup!");
  }
  INFO("DONE!");
  INFO("Build Basis in each U(1) sector - ");
  std::vector<Basis> Bases;
  for (size_t cntN = 0; cntN <= N; cntN++) {
    Basis B1(L, cntN);
    B1.Boson();
    Bases.push_back(B1);
  }
  INFO("DONE!");
  INFO("Build Hamiltonian in each U(1) sector - ");
  std::vector< std::vector<ComplexType> > Vloc;
  std::vector<ComplexType> Vtmp(L, (ComplexType)Vin);
  Vloc.push_back(Vtmp);
  std::vector< std::vector<ComplexType> > Uloc;
  std::vector<ComplexType> Utmp(L, (ComplexType)Uin);
  Uloc.push_back(Utmp);
  std::vector<Hamiltonian<ComplexType> > Hams;
  for ( auto &bs : Bases ){
    Hamiltonian<ComplexType> ham( std::vector<Basis>(1, bs) );
    ham.BuildLocalHamiltonian(Vloc, Uloc, std::vector<Basis>(1, bs));
    ham.BuildHoppingHamiltonian(std::vector<Basis>(1, bs), lattice);
    ham.BuildTotalHamiltonian();
    Hams.push_back(ham);
  }
  INFO("Build Initial Density Matrix - ");
  std::vector<ComplexMatrixType> Rhos;
  for ( auto &bs : Bases ){
    ComplexMatrixType Rho = ComplexMatrixType::Zero(bs.getHilbertSpace(), bs.getHilbertSpace());
    if ( bs.getN() == N ) {
      GetRho(bs, Rho);
    }
    Rhos.push_back(Rho);
  }
  INFO("DONE!");
  std::vector<ComplexType> Nbi = Ni( Bases, Rhos );
  ComplexMatrixType Nij = NiNj( Bases, Rhos );
  std::string output_file = "TBLB";
  output_file.append(std::to_string(gamma));
  output_file.append(".h5");
  SaveObs(output_file, "Obs-0", Nbi, Nij);
  std::vector<std::vector<size_t> > CIdx = IndexCollapse(TBloc, Bases);
  INFO("Trace = " << TraceRhos(Rhos));
  for (size_t cntT = 1; cntT <= Tsteps; cntT++) {
    Lindblad_RK4( dt, gamma, TBloc, Bases, Hams, CIdx, Rhos);
    // Lindblad_Newton( dt, gamma, TBloc, Bases, Hams, CIdx, Rhos);
    /* NOTE: Check trace og Rhos */
    RealType tr = TraceRhos(Rhos);
    INFO("Trace = " << tr);
    /* NOTE: Expectation values */
    Nbi = Ni( Bases, Rhos );
    Nij = NiNj( Bases, Rhos );
    /* NOTE: H5 group name */
    std::string gname = "Obs-";
    gname.append(std::to_string((unsigned long long)cntT));
    gname.append("/");
    SaveObs(output_file, gname, Nbi, Nij);
    // std::cin.get();
  }
  return 0;
}

RealType TraceRhos(const std::vector<ComplexMatrixType> &Rhos){
  ComplexType tr = ComplexType(0.0e0, 0.0e0);
  for ( auto &rho: Rhos){
    tr += rho.trace();
  }
  assert( std::abs(tr.imag()) < 1.0e-12 );
  return tr.real();
}

std::vector<std::vector<size_t> > IndexCollapse(const size_t TBloc,
  const std::vector<Basis> &bs){
  std::vector<std::vector<size_t> > idx;
  for (size_t cnt = 0; cnt < bs.size() - 1; cnt++) {
    std::vector<size_t> tmp_idx;
    for ( auto &eb: bs.at(cnt).getBStates() ) {
      std::vector<int> v = eb;
      v.at(TBloc) += 1;
      tmp_idx.push_back(bs.at(cnt+1).getIndexFromTag(BosonBasisTag(v)));
    }
    idx.push_back(tmp_idx);
  }
  return idx;
}

void SaveObs( const std::string filename, const std::string gname,
  const std::vector<ComplexType> &v, const ComplexMatrixType &m){
  HDF5IO file(filename);
  file.saveStdVector(gname, "Nb", v);
  ComplexType Ntot = std::accumulate(v.begin(), v.end(), ComplexType(0.0, 0.0));
  file.saveNumber(gname, "Ntot", Ntot.real());
  file.saveMatrix(gname, "Nij", m);
}

void GetRho( const Basis &bs, ComplexMatrixType &Mat ){
  HDF5IO gsf("BSSH.h5");
  ComplexVectorType gswf;
  gsf.loadVector("GS", "EVec", gswf);
  if ( DEBUG ){
    int oldL = gsf.loadUlong("1DChain", "L");
    int oldN = gsf.loadUlong("1DChain", "N");
    assert( oldL == bs.getL() );
    assert( oldN == bs.getN() );
  }
  Mat = gswf * gswf.adjoint();
}

void LoadParameters( const std::string filename, int &L, int &OBC, int &N,
  RealType &Uloc, RealType &Vloc, RealType &dt, int &Tsteps, int &TBloc,
  RealType &Gamma){
    HDF5IO file(filename);
    L = file.loadInt("Parameters", "L");
    OBC = file.loadInt("Parameters", "OBC");
    N = file.loadInt("Parameters", "N");
    Uloc = file.loadReal("Parameters", "U");
    Vloc = file.loadReal("Parameters", "V");
    dt = file.loadReal("Parameters", "dt");
    Gamma = file.loadReal("Parameters", "Gamma");
    Tsteps = file.loadInt("Parameters", "Tsteps");
    TBloc = file.loadInt("Parameters", "TBloc");
}

std::vector<ComplexType> Ni( const std::vector<Basis> &Bases,
  const std::vector<ComplexMatrixType> &Rhos ){
  std::vector<ComplexType> tmp(Bases.at(0).getL(), 0.0e0);
  for (size_t cnt = 0; cnt < Bases.size(); cnt++) {
    std::vector< std::vector<int> > b = Bases.at(cnt).getBStates();
    assert( b.size() == Rhos.at(cnt).cols() );
    assert( b.size() == Rhos.at(cnt).rows() );
    int coff = 0;
    for ( auto &nbi : b ){
      for (size_t cntL = 0; cntL < Bases.at(cnt).getL(); cntL++) {
        tmp.at(cntL) += (RealType)nbi.at(cntL) * Rhos.at(cnt)(coff, coff);
      }
      coff++;
    }
  }
  return tmp;
}

ComplexMatrixType NiNj( const std::vector<Basis> &Bases,
  const std::vector<ComplexMatrixType> &Rhos ){
  ComplexMatrixType tmp = ComplexMatrixType::Zero(Bases.at(0).getL(), Bases.at(0).getL());
  for (size_t cnt = 0; cnt < Bases.size(); cnt++) {
    std::vector< std::vector<int> > b = Bases.at(cnt).getBStates();
    assert( b.size() == Rhos.at(cnt).cols() );
    assert( b.size() == Rhos.at(cnt).rows() );
    int coff = 0;
    for ( auto &nbi : b ){
      for (size_t cnt1 = 0; cnt1 < Bases.at(cnt).getL(); cnt1++) {
        for (size_t cnt2 = 0; cnt2 < Bases.at(cnt).getL(); cnt2++) {
          tmp(cnt1, cnt2) += (RealType)nbi.at(cnt1) * (RealType)nbi.at(cnt2) *
            Rhos.at(cnt)(coff, coff);
        }
      }
      coff++;
    }
  }
  return tmp;
}
