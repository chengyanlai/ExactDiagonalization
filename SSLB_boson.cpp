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
#include "src/Dephasing/Dephasing.hpp"

#ifdef MKL
  #include "mkl.h"
#endif

#ifndef DEBUG
#define DEBUG 1
#endif

#ifndef NumCores
#define NumCores 2
#endif

void LoadParameters( const std::string filename, int &L, int &OBC, int &N,
  RealType &Uloc, std::vector<RealType> &Veqm, std::vector<RealType> &Vdyn,
  RealType &dt, int &Tsteps, std::vector<size_t> &Sites, RealType &Gamma);
void SaveObs( const std::string filename, const std::string gname,
  const std::vector<ComplexType> &v, const ComplexMatrixType &m,
  const ComplexMatrixType &cm, const ComplexType val);
std::vector<ComplexType> Ni( const std::vector<Basis> &Bases,
  const ComplexMatrixType &Rhos );
ComplexMatrixType NiNj( const std::vector<Basis> &Bases,
  const ComplexMatrixType &Rhos );
ComplexMatrixType OPCM( const std::vector<Basis> &Bases,
  const ComplexMatrixType &Rhos);
ComplexType Energy( const Hamiltonian<ComplexType> &ham,
  const ComplexMatrixType &rho );

int main(int argc, char const *argv[]) {
  Eigen::setNbThreads(NumCores);
#ifdef MKL
  mkl_set_num_threads(NumCores);
#endif
  INFO("Eigen3 uses " << Eigen::nbThreads() << " threads.");
  int L, OBC, N, Tsteps;
  RealType Uin, dt, gamma;
  std::vector<RealType> Veqm, Vdyn;
  std::vector<size_t> Sites;
  LoadParameters( "conf.h5", L, OBC, N, Uin, Veqm, Vdyn, dt, Tsteps, Sites, gamma );
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
  INFO("Build Basis - ");
  std::vector<Basis> Bases;
  Basis B1(L, N);
  B1.Boson();
  Bases.push_back(B1);
  INFO("DONE!");
  INFO("Build Hamiltonian - ");
  std::vector< std::vector<ComplexType> > Vloc;
  INFO("Initial potential set to ");
  std::vector<ComplexType> Vtmp;
  for ( RealType &val : Veqm ){
    Vtmp.push_back((ComplexType)val);
    INFO_NONEWLINE(" " << val);
  }
  INFO(" ");
  Vloc.push_back(Vtmp);
  std::vector< std::vector<ComplexType> > Uloc;
  std::vector<ComplexType> Utmp(L, (ComplexType)Uin);
  Uloc.push_back(Utmp);
  Hamiltonian<ComplexType> ham( Bases );
  ham.BuildLocalHamiltonian(Vloc, Uloc, Bases);
  ham.BuildHoppingHamiltonian(Bases, lattice);
  ham.BuildTotalHamiltonian();
  INFO("DONE!");
  INFO("Diagonalize Hamiltonian - ");
  std::vector<RealType> Val;
  Hamiltonian<ComplexType>::VectorType Vec;
  ham.eigh(Val, Vec);
  INFO("E_0 = " << Val.at(0) << ", E_1 = " << Val.at(1));
  INFO("DONE!");
  INFO("Build Initial Density Matrix - ");
  ComplexMatrixType Rho = Vec * Vec.adjoint();
  INFO("Trace = " << Rho.trace());
  INFO("DONE!");
  std::vector<ComplexType> Nbi = Ni( Bases, Rho );
  ComplexMatrixType Nij = NiNj( Bases, Rho );
  ComplexMatrixType CM = OPCM( Bases, Rho );
  ComplexType Eavg = Energy(ham, Rho);
  std::string output_file = "SSLB.h5";
  SaveObs(output_file, "Obs-0", Nbi, Nij, CM, Eavg);

  INFO("Quench the local potetnial to");
  Vtmp.clear();
  for ( RealType &val : Vdyn ){
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
    Lindblad_RK4( dt, gamma, Sites, Bases, ham, Rho );
    /* NOTE: Check trace of Rhos */
    INFO("Time = " << cntT * dt << ", Trace = " << Rho.trace());
    /* NOTE: Expectation values */
    Nbi = Ni( Bases, Rho );
    Nij = NiNj( Bases, Rho );
    CM = OPCM( Bases, Rho );
    Eavg = Energy(ham, Rho);
    /* NOTE: H5 group name */
    std::string gname = "Obs-";
    gname.append(std::to_string((unsigned long long)cntT));
    gname.append("/");
    SaveObs(output_file, gname, Nbi, Nij, CM, Eavg);
  }
  return 0;
}

void SaveObs( const std::string filename, const std::string gname,
  const std::vector<ComplexType> &v, const ComplexMatrixType &m,
  const ComplexMatrixType &cm, const ComplexType val){
  HDF5IO file(filename);
  file.saveStdVector(gname, "Nb", v);
  ComplexType Ntot = std::accumulate(v.begin(), v.end(), ComplexType(0.0, 0.0));
  file.saveNumber(gname, "Ntot", Ntot.real());
  file.saveMatrix(gname, "Nij", m);
  file.saveMatrix(gname, "OPCM", cm);
  file.saveNumber(gname, "E", val);
}

void LoadParameters( const std::string filename, int &L, int &OBC, int &N,
  RealType &Uloc, std::vector<RealType> &Veqm, std::vector<RealType> &Vdyn,
  RealType &dt, int &Tsteps, std::vector<size_t> &Sites, RealType &Gamma){
    HDF5IO file(filename);
    L = file.loadInt("Parameters", "L");
    OBC = file.loadInt("Parameters", "OBC");
    N = file.loadInt("Parameters", "N");
    Uloc = file.loadReal("Parameters", "U");
    dt = file.loadReal("Parameters", "dt");
    Gamma = file.loadReal("Parameters", "Gamma");
    Tsteps = file.loadInt("Parameters", "Tsteps");
    file.loadStdVector("Parameters", "Veqm", Veqm);
    file.loadStdVector("Parameters", "Vdyn", Vdyn);
    file.loadStdVector("Parameters", "Sites", Sites);
}

std::vector<ComplexType> Ni( const std::vector<Basis> &Bases,
  const ComplexMatrixType &Rhos ){
  std::vector<ComplexType> tmp(Bases.at(0).getL(), 0.0e0);
  std::vector< std::vector<int> > b = Bases.at(0).getBStates();
  assert( b.size() == Rhos.cols() );
  assert( b.size() == Rhos.rows() );
  int coff = 0;
  for ( auto &nbi : b ){
    for (size_t cntL = 0; cntL < Bases.at(0).getL(); cntL++) {
      tmp.at(cntL) += (RealType)nbi.at(cntL) * Rhos(coff, coff);
    }
    coff++;
  }
  return tmp;
}

ComplexMatrixType NiNj( const std::vector<Basis> &Bases,
  const ComplexMatrixType &Rhos ){
  ComplexMatrixType tmp = ComplexMatrixType::Zero(Bases.at(0).getL(), Bases.at(0).getL());
  std::vector< std::vector<int> > b = Bases.at(0).getBStates();
  assert( b.size() == Rhos.cols() );
  assert( b.size() == Rhos.rows() );
  int coff = 0;
  for ( auto &nbi : b ){
    for (size_t cnt1 = 0; cnt1 < Bases.at(0).getL(); cnt1++) {
      for (size_t cnt2 = 0; cnt2 < Bases.at(0).getL(); cnt2++) {
        tmp(cnt1, cnt2) += (RealType)nbi.at(cnt1) * (RealType)nbi.at(cnt2) *
          Rhos(coff, coff);
      }
    }
    coff++;
  }
  return tmp;
}

ComplexMatrixType OPCM( const std::vector<Basis> &Bases,
  const ComplexMatrixType &Rhos){
  /* NOTE: Calculate <c^\dagger_i c_j> = Tr( \rho * c^\dagger_i c_j ) */
  ComplexMatrixType tmp = ComplexMatrixType::Zero(Bases.at(0).getL(), Bases.at(0).getL());
  std::vector< std::vector<int> > b = Bases.at(0).getBStates();
  assert( b.size() == Rhos.cols() );
  assert( b.size() == Rhos.rows() );
  int state_id1 = 0;
  for ( auto &nbi : b ){
    for (size_t site_i = 0; site_i < Bases.at(0).getL(); site_i++) {
      for (size_t site_j = site_i; site_j < Bases.at(0).getL(); site_j++) {
        if ( nbi.at(site_j) > 0 ) {
          std::vector<int> nbj = nbi;
          RealType val = (RealType)nbj.at(site_j);
          nbj.at(site_j) = nbj.at(site_j) - 1;
          nbj.at(site_i) = nbj.at(site_i) + 1;
          val *= (RealType)nbj.at(site_i);
          size_t state_id2 = Bases.at(0).getIndexFromTag( BosonBasisTag(nbj) );
          tmp(site_i,site_j) += sqrt(val) * Rhos(state_id1,state_id2);
        }
      }
    }
    state_id1 += 1;
  }
  return tmp + tmp.adjoint();
}

ComplexType Energy( const Hamiltonian<ComplexType> &ham,
  const ComplexMatrixType &rho ){
  ComplexSparseMatrixType h = ham.getTotalHamiltonian();
  return (rho * h).trace();
}
