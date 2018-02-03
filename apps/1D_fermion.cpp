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
#define DT ComplexType
#define DTV ComplexVectorType
#define DTM ComplexMatrixType
// #else
// #define DT RealType
// #define DTV RealVectorType
// #define DTM RealMatrixType
// #endif

void LoadParameters( const std::string filename, int &L, int &OBC, int &N1, int &N2,
  RealType &J1, RealType &J2, RealType &Phase, std::vector<RealType> &Uls, std::vector<RealType> &Vls,
  int &dynamics, int &Tsteps, RealType &dt){
    HDF5IO file(filename);
    L = file.loadInt("Parameters", "L");
    OBC = file.loadInt("Parameters", "OBC");
    N1 = file.loadInt("Parameters", "N1");
    N2 = file.loadInt("Parameters", "N2");
    J1 = file.loadReal("Parameters", "J1");
    J2 = file.loadReal("Parameters", "J2");
    Phase = file.loadReal("Parameters", "Phase");
    file.loadStdVector("Parameters", "U", Uls);
    file.loadStdVector("Parameters", "V", Vls);
    dynamics = file.loadInt("Parameters", "dynamics");
    Tsteps = file.loadInt("Parameters", "Tsteps");
    dt = file.loadReal("Parameters", "dt");
}

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

ComplexMatrixType SingleParticleDensityMatrix( const int species, const std::vector<Basis> &Bases, const ComplexVectorType &Vec, Hamiltonian<ComplexType> &ham );
std::vector< DTV > Ni( const std::vector<Basis> &Bases, const DTV &Vec, Hamiltonian<DT> &ham );
DTM NupNdn( const std::vector<Basis> &Bases, const DTV &Vec, Hamiltonian<DT> &ham );
DTM NupNup( const std::vector<Basis> &Bases, const DTV &Vec, Hamiltonian<DT> &ham );
DTM NdnNdn( const std::vector<Basis> &Bases, const DTV &Vec, Hamiltonian<DT> &ham );

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
  RealType J1in, J2in, Phase;// Benzene
  // Load parameters from file
  try{
    H5::Exception::dontPrint();
    H5::H5File::isHdf5("conf.h5");
    // LoadParameters( "conf.h5", L, OBC, N1, N2, Jin, Uin, Vin, dynamics, Tsteps, dt);
    LoadParameters( "conf.h5", L, OBC, N1, N2, J1in, J2in, Phase, Uin, Vin, dynamics, Tsteps, dt);// Benzene
    std::cout << "Parameters for calculation loaded!" << std::endl;
  }catch(H5::FileIException){
    L = 6;
    OBC = 1;
    N1 = 3;
    N2 = 3;
    Jin = std::vector<RealType>(L-1, 1.0);
    Uin = std::vector<RealType>(L, 1.0);
    Vin = std::vector<RealType>(L, 0.0);
    // J1in = 1.00;// Benzene
    // J2in = 0.50;// Benzene
    // Phase = 5.0 / 6.0 - 0.035;// Benzene
    dynamics = 0;
    Tsteps = 0;
    dt = 0.005;
  }
  HDF5IO *file = new HDF5IO("FData.h5");
  INFO("Build Lattice - ");
  std::vector<DT> J( Jin.begin(), Jin.end() );
  INFO("");
  /* Temporory for Benzene */
  // DT J1 = J1in * exp(ComplexType(0.0, Phase)*PI);// Benzene
  // DT J2 = J2in * exp(ComplexType(0.0, Phase)*PI);// Benzene
  // DT J2 = J2in;// Benzene
  // std::vector< Node<DT>* > lattice = NNN_1D_Chain(L, J1, J2, OBC);// Benzene
  std::vector< Node<DT>* > lattice = NN_1D_Chain(L, J, OBC);
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
  ham.CheckTotalHamiltonian();
  INFO("DONE!");
  INFO_NONEWLINE("Diagonalize Hamiltonian - ");
  RealVectorType Vals;
  DTM Vecs;
  // ham.eigh(Vals, Vecs, 2);
  ham.diag(Vals, Vecs);
  INFO("GS energy = " << Vals[0]);
  INFO("FES energy = " << Vals[1]);
  DTV Vec = Vecs.row(0);
  file->saveVector("GS", "EVec", Vec);
  file->saveVector("GS", "EVal", Vals);
  INFO("DONE!");
  std::vector< DTV > Nfi = Ni( Bases, Vec, ham );
  INFO(" Up Spin - ");
  INFO(Nfi.at(0));
  INFO(" Down Spin - ");
  INFO(Nfi.at(1));
  INFO(" N_i - ");
  DTV Niall = Nfi.at(0) + Nfi.at(1);
  INFO(Niall);
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
  ComplexMatrixType CMUp = SingleParticleDensityMatrix( 0, Bases, Vec, ham );
  INFO(CMUp);
  ComplexMatrixType CMDn = SingleParticleDensityMatrix( 1, Bases, Vec, ham );
  INFO(CMDn);
  file->saveVector("Obs", "Nup", Nfi.at(0));
  file->saveVector("Obs", "Ndn", Nfi.at(1));
  file->saveMatrix("Obs", "NupNdn", Nud);
  file->saveMatrix("Obs", "NupNup", Nuu);
  file->saveMatrix("Obs", "NdnNdn", Ndd);
  file->saveMatrix("Obs", "CMUp", CMUp);
  file->saveMatrix("Obs", "CMDn", CMDn);
  delete file;
  return 0;
}

ComplexMatrixType SingleParticleDensityMatrix( const int species, const std::vector<Basis> &Bases, const ComplexVectorType &Vec, Hamiltonian<ComplexType> &ham ){
  size_t L = Bases.at(species).getL();
  ComplexMatrixType CM = ComplexMatrixType::Zero(L,L);
  std::vector< int > bs = Bases.at(species).getFStates();
  std::vector<size_t> tg = Bases.at(species).getFTags();
  for ( const int &b : bs ){
    size_t bid = tg.at(b);// Find their indices
    for ( size_t i=0; i < L; i++){
      for ( size_t j=i; j < L; j++){
        /* see if hopping exist */
        if ( btest(b, j) && !(btest(b, i)) ) {
          /* c^\dagger_i c_j if yes, no particle in i and one particle in j. */
          int CrossFermionNumber = 0;
          ComplexType tsign = (ComplexType)(1.0e0);
          if ( j - i > 1 ){
            // possible cross fermions and sign change.
            for ( int k = i+1; k < j; k++){
              CrossFermionNumber += btest(b, k);
            }
          }
          if (CrossFermionNumber % 2 == 1) tsign = (ComplexType)(-1.0e0);
          int p = ibset(b, i);
          p = ibclr(p, j);
          size_t pid = tg.at(p);// Find their indices
          size_t count;
          if ( species == 0 ) count = Bases.at(1).getHilbertSpace();
          else if ( species == 1 ) count = Bases.at(0).getHilbertSpace();
          std::vector<size_t> rids(2, bid);
          std::vector<size_t> cids(2, pid);
          for (size_t loop_id = 0; loop_id < count; loop_id++) {
            if ( species == 0 ){
              rids.at(1) = loop_id;
              cids.at(1) = loop_id;
            }else if ( species == 1 ){
              rids.at(0) = loop_id;
              cids.at(0) = loop_id;
            }
            size_t rid = ham.DetermineTotalIndex( rids );
            size_t cid = ham.DetermineTotalIndex( cids );
            CM(i, j) +=  tsign * Vec(cid) * std::conj( Vec(rid) );//Vec(id) * std::conj( Vec(id) );
          }
        }
      }
    }
  }
  return CM;
}

std::vector< DTV > Ni( const std::vector<Basis> &Bases, const DTV &Vec, Hamiltonian<DT> &ham ){
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
