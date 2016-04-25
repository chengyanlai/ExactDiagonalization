#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <numeric>
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

std::vector< std::vector<ComplexType> >  UpdateV(const std::vector<ComplexType> Veqm,
  const std::vector<ComplexType> Vfin, const int cntT, const int T0step, const RealType dt);
void SaveObs( const std::string filename, const std::string gname,
  const std::vector<ComplexType> &v, const ComplexMatrixType &m);
void LoadEqmParameters( const std::string filename, int &L, int &OBC,
  int &N, RealType &Uloc, std::vector<RealType> &Vloc);
void LoadDynParameters( const std::string filename, RealType &dt, int &Tsteps,
  RealType &Uloc, std::vector<RealType> &Vloc, int &T0step);
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
  const std::vector< Node<ComplexType>* > lattice = NN_1D_Chain(L, J, OBC);
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
  Hamiltonian<ComplexType> ham( Bases );
  std::vector< std::vector<ComplexType> > Veqm;
  std::vector<ComplexType> Vtmp;//(L, 1.0);
  for ( RealType &val : Vin ){
    Vtmp.push_back((ComplexType)val);
  }
  Veqm.push_back(Vtmp);
  std::vector< std::vector<ComplexType> > Uloc;
  std::vector<ComplexType> Utmp(L, (ComplexType)Uin);
  Uloc.push_back(Utmp);
  ham.BuildLocalHamiltonian(Veqm, Uloc, Bases);
  ham.BuildHoppingHamiltonian(Bases, lattice);
  ham.BuildTotalHamiltonian();
  INFO("DONE!");
  INFO_NONEWLINE("Diagonalize Hamiltonian - ");
  std::vector<RealType> Val;
  Hamiltonian<ComplexType>::VectorType Vec;
  ham.eigh(Val, Vec);
  INFO("GS energy = " << Val.at(0));
  HDF5IO file("GS.h5");
  file.saveNumber("1DChain", "L", L);
  file.saveStdVector("1DChain", "J", J);
  file.saveNumber("Basis", "N", N);
  file.saveVector("GS", "EVec", Vec);
  file.saveStdVector("GS", "EVal", Val);
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
  int T0step;
  RealType dt;
  Vin.clear();
  LoadDynParameters( "Dyn.h5", dt, Tsteps, Uin, Vin, T0step);
  Vtmp.clear();
  INFO("Final local potetnial set to");
  for ( RealType &val : Vin ){
    Vtmp.push_back((ComplexType)val);
    INFO_NONEWLINE(" " << val);
  }
  INFO("");
  INFO(" with time scale " << T0step * dt);
  std::vector< std::vector<ComplexType> > Vfin;
  Vfin.push_back(Vtmp);
  INFO("");
  INFO("Quench the global interaction to " << Uin);
  Utmp.assign(L, (ComplexType)Uin);
  Uloc.clear();
  Uloc.push_back(Utmp);
  INFO("Start Time-Evolution");
  ComplexType Prefactor = ComplexType(0.0, -1.0e0*dt);/* NOTE: hbar = 1 */
  for (size_t cntT = 1; cntT <= Tsteps; cntT++) {
    if ( cntT == T0step ){
      INFO("Update the local/total Hamiltonian (final).");
      ham.BuildLocalHamiltonian(Vfin, Uloc, Bases);
      ham.BuildTotalHamiltonian();
    }else if ( cntT < T0step ){
      INFO("Update the local/total Hamiltonian");
      std::vector< std::vector<ComplexType> > Vwork = UpdateV(Veqm.at(0), Vfin.at(0), cntT, T0step, dt);
      ham.BuildLocalHamiltonian(Vwork, Uloc, Bases);
      ham.BuildTotalHamiltonian();
    }
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

std::vector< std::vector<ComplexType> >  UpdateV(const std::vector<ComplexType> V1,
  const std::vector<ComplexType> V2, const int cntT, const int T0step, const RealType dt){
  std::vector< std::vector<ComplexType> > out;
  std::vector<ComplexType> tmp;
  RealType ratio = (RealType)cntT / (RealType)T0step;
  for (size_t i = 0; i < V1.size(); i++) {
    ComplexType val1 = V1.at(i) * (1.0 - ratio);
    ComplexType val2 = V2.at(i) * ratio;
    INFO_NONEWLINE(val1 + val2);
    tmp.push_back(val1 + val2);
  }
  INFO("");
  out.push_back(tmp);
  return out;
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
  RealType &Uloc, std::vector<RealType> &Vloc, int &T0step){
    HDF5IO file(filename);
    dt = file.loadReal("Parameters", "dt");
    Tsteps = file.loadInt("Parameters", "Tsteps");
    Uloc = file.loadReal("Parameters", "U");
    file.loadStdVector("Parameters", "V", Vloc);
    T0step = file.loadInt("Parameters", "T0step");
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