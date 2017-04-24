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
#define NumCores 10
#endif

void LoadParameters( const std::string filename, int &L, RealType &SOC,
  int &N1, int &N2, RealType &Uloc,
  std::vector<RealType> &Vloc,
  int &dynamics, int &Tsteps, RealType &ComplexType){
    HDF5IO file(filename);
    L = file.loadInt("Parameters", "L");
    SOC = file.loadReal("Parameters", "SOC");
    N1 = file.loadInt("Parameters", "N1");
    N2 = file.loadInt("Parameters", "N2");
    Uloc = file.loadReal("Parameters", "U");
    file.loadStdVector("Parameters", "V", Vloc);
    dynamics = file.loadInt("Parameters", "dynamics");
    Tsteps = file.loadInt("Parameters", "Tsteps");
    ComplexType = file.loadReal("Parameters", "ComplexType");
}

ComplexMatrixType SingleParticleDensityMatrix( const int species, const std::vector<Basis> &Bases,
  const ComplexVectorType &Vec, Hamiltonian<ComplexType> &ham ){
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
          ComplexType tsign = ComplexType(1.0e0, 0.0e0);
          if ( j - i > 1 ){
            // possible cross fermions and sign change.
            for ( int k = i+1; k < j; k++){
              CrossFermionNumber += btest(b, k);
            }
          }
          if (CrossFermionNumber % 2 == 1) tsign = ComplexType(-1.0e0, 0.0e0);
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

std::vector<ComplexVectorType> Ni( const std::vector<Basis> &Bases, const ComplexVectorType &Vec,
  Hamiltonian<ComplexType> &ham );
ComplexMatrixType NupNdn( const std::vector<Basis> &Bases, const ComplexVectorType &Vec,
  Hamiltonian<ComplexType> &ham );
ComplexMatrixType NupNup( const std::vector<Basis> &Bases, const ComplexVectorType &Vec,
  Hamiltonian<ComplexType> &ham );
ComplexMatrixType NdnNdn( const std::vector<Basis> &Bases, const ComplexVectorType &Vec,
  Hamiltonian<ComplexType> &ham );

int main(int argc, char const *argv[]) {
  Eigen::setNbThreads(NumCores);
#ifdef MKL
  mkl_set_num_threads(NumCores);
#endif
  INFO("Eigen3 uses " << Eigen::nbThreads() << " threads.");
  const int L = 6;
  const RealType lambda = 1.0e0;
  const RealType phi = 0.20e0 * PI;
  const int N1 = 4, N2 = 4;
  const int dynamics = 0, Tsteps = 0;
  const RealType Uin = 0.0, dt = 0.0;
  const std::vector<RealType> Vin(L, 0.0e0);
  // LoadParameters( "conf.h5", L, SOC, N1, N2, Uin, Vin, dynamics, Tsteps, dt);
  HDF5IO *file = new HDF5IO("out.h5");
  INFO("Build Lattice - ");
  std::vector<ComplexType> J(L, 1.0);
  file->saveNumber("1DChain", "L", L);
  file->saveStdVector("1DChain", "J", J);
  std::vector< Node<ComplexType>* > latticeUP = NN_1D_Chain(L, J, false);// must be PBC
  // Add the SOC links
  const ComplexType SOC = lambda * exp(ComplexType(0.0e0, 1.0e0) * phi);
  // ComplexType vSOC = lambda * exp(ComplexType(0.0e0, 1.0e0) * phi);
  // ComplexType SOC = vSOC;
  // if ( std::abs(vSOC.real()) < 1.0e-12 ) SOC = ComplexType(0.0e0, vSOC.imag());
  // if ( std::abs(vSOC.imag()) < 1.0e-12 ) SOC = ComplexType(vSOC.real(), 0.0e0);
  latticeUP.at(0)->LinkTo(latticeUP.at(2), SOC);
  // latticeUP.at(2)->LinkTo(latticeUP.at(4), SOC);
  // latticeUP.at(4)->LinkTo(latticeUP.at(0), SOC);
  // latticeUP.at(1)->LinkTo(latticeUP.at(3), SOC);
  // latticeUP.at(3)->LinkTo(latticeUP.at(5), SOC);
  // latticeUP.at(5)->LinkTo(latticeUP.at(1), SOC);
  for ( auto &lt : latticeUP ){
    if ( !(lt->VerifySite()) ) RUNTIME_ERROR("Wrong latticeUP setup!");
  }
  std::vector< Node<ComplexType>* > latticeDN = NN_1D_Chain(L, J, false);// must be PBC
  // Add the SOC links
  latticeDN.at(0)->LinkTo(latticeDN.at(2),-SOC);
  // latticeDN.at(2)->LinkTo(latticeDN.at(4),-SOC);
  // latticeDN.at(4)->LinkTo(latticeDN.at(0),-SOC);
  // latticeDN.at(1)->LinkTo(latticeDN.at(3),-SOC);
  // latticeDN.at(3)->LinkTo(latticeDN.at(5),-SOC);
  // latticeDN.at(5)->LinkTo(latticeDN.at(1),-SOC);
  for ( auto &lt : latticeDN ){
    if ( !(lt->VerifySite()) ) RUNTIME_ERROR("Wrong latticeDN setup!");
  }
  std::vector< std::vector< Node<ComplexType>* > > lattice;
  lattice.push_back(latticeUP);
  lattice.push_back(latticeDN);
  INFO("DONE!");
  INFO("Build Basis - ");
  Basis FUP(L, N1, true);
  FUP.Fermion();
  std::vector<int> st1 = FUP.getFStates();
  std::vector<size_t> tg1 = FUP.getFTags();
  Basis FDN(L, N2, true);
  FDN.Fermion();
  std::vector<int> st2 = FDN.getFStates();
  std::vector<size_t> tg2 = FDN.getFTags();
  file->saveNumber("Basis", "N1", N1);
  file->saveStdVector("Basis", "F1States", st1);
  file->saveStdVector("Basis", "F1Tags", tg1);
  file->saveNumber("Basis", "N2", N2);
  file->saveStdVector("Basis", "F2States", st2);
  file->saveStdVector("Basis", "F2Tags", tg2);
  INFO("DONE!");
  INFO_NONEWLINE("Build Hamiltonian - ");
  std::vector<Basis> Bases;
  Bases.push_back(FUP);
  Bases.push_back(FDN);
  Hamiltonian<ComplexType> ham( Bases );
  std::vector< std::vector<ComplexType> > Vloc;
  std::vector<ComplexType> Vtmp;//(L, 1.0);
  for ( const RealType &val : Vin ){
    Vtmp.push_back((ComplexType)val);
  }
  Vloc.push_back(Vtmp);
  Vloc.push_back(Vtmp);
  std::vector< std::vector<ComplexType> > Uloc;
  std::vector<ComplexType> Utmp(L, (ComplexType)Uin);
  Uloc.push_back(Utmp);
  Uloc.push_back(Utmp);
  ham.BuildLocalHamiltonian(Vloc, Uloc, Bases);
  INFO(" - BuildLocalHamiltonian DONE!");
  ham.BuildHoppingHamiltonian(Bases, lattice);
  INFO(" - BuildHoppingHamiltonian DONE!");
  ham.BuildTotalHamiltonian();
  INFO("DONE!");
  INFO_NONEWLINE("Diagonalize Hamiltonian - ");
  std::vector<RealType> Val;
  Hamiltonian<ComplexType>::VectorType Vec;
  RealVectorType AllVal;
  ComplexMatrixType AllVec;
  ham.diag(AllVal, AllVec);
  file->saveVector("GS", "EValAll", AllVal);
  file->saveMatrix("GS", "EVecAll", AllVec);
  Val.push_back(AllVal(0));
  Vec = AllVec.row(0);
  // std::cout << "Eigenenergies = " << AllVal << std::endl;
  std::cout << "Eigenenergies = " << AllVal(0) << " " << AllVal(1) << std::endl;
  // INFO("GS energy = " << Val.at(0));
  file->saveVector("GS", "EVec", Vec);
  file->saveStdVector("GS", "EVal", Val);
  INFO("DONE!");
  // Observables
  std::vector< ComplexVectorType > Nfi = Ni( Bases, Vec, ham );
  INFO(" Up Spin - ");
  INFO(Nfi.at(0));
  INFO(" Down Spin - ");
  INFO(Nfi.at(1));
  INFO(" N_i - ");
  ComplexVectorType Niall = Nfi.at(0) + Nfi.at(1);
  INFO(Niall);
  ComplexMatrixType Nud = NupNdn( Bases, Vec, ham );
  INFO(" Correlation NupNdn");
  INFO(Nud);
  ComplexMatrixType Nuu = NupNup( Bases, Vec, ham );
  INFO(" Correlation NupNup");
  INFO(Nuu);
  ComplexMatrixType Ndd = NdnNdn( Bases, Vec, ham );
  INFO(" Correlation NdnNdn");
  INFO(Ndd);
  INFO(" N_i^2 - ");
  ComplexMatrixType Ni2 = Nuu.diagonal() + Ndd.diagonal() + 2.0e0 * Nud.diagonal();
  INFO(Ni2);
  INFO(" CMup - ");
  ComplexMatrixType CMup = SingleParticleDensityMatrix( 0, Bases, Vec, ham );
  INFO(CMup);
  INFO(" CMdn - ");
  ComplexMatrixType CMdn = SingleParticleDensityMatrix( 1, Bases, Vec, ham );
  INFO(CMdn);
  file->saveVector("Obs", "Nup", Nfi.at(0));
  file->saveVector("Obs", "Ndn", Nfi.at(1));
  file->saveMatrix("Obs", "NupNdn", Nud);
  file->saveMatrix("Obs", "NupNup", Nuu);
  file->saveMatrix("Obs", "NdnNdn", Ndd);
  file->saveMatrix("Obs", "CMup", CMup);
  file->saveMatrix("Obs", "CMdn", CMdn);
  delete file;
  return 0;
}

std::vector< ComplexVectorType > Ni( const std::vector<Basis> &Bases,
  const ComplexVectorType &Vec, Hamiltonian<ComplexType> &ham ){
  std::vector< ComplexVectorType > out;
  ComplexVectorType tmp1 = ComplexVectorType::Zero(Bases.at(0).getL());//(Bases.at(0).getL(), 0.0e0);
  ComplexVectorType tmp2 = ComplexVectorType::Zero(Bases.at(1).getL());//(Bases.at(1).getL(), 0.0e0);
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

ComplexMatrixType NupNdn( const std::vector<Basis> &Bases,
  const ComplexVectorType &Vec, Hamiltonian<ComplexType> &ham ){
  ComplexMatrixType out = ComplexMatrixType::Zero(
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

ComplexMatrixType NupNup( const std::vector<Basis> &Bases,
  const ComplexVectorType &Vec, Hamiltonian<ComplexType> &ham ){
  ComplexMatrixType out = ComplexMatrixType::Zero(Bases.at(0).getL(), Bases.at(0).getL() );
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

ComplexMatrixType NdnNdn( const std::vector<Basis> &Bases,
  const ComplexVectorType &Vec, Hamiltonian<ComplexType> &ham ){
  ComplexMatrixType out = ComplexMatrixType::Zero(
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
