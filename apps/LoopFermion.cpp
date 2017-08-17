#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include "src/EDType.hpp"
#include "src/bitwise.h"
#include "src/Node/Node.hpp"
#include "src/Basis/Basis.hpp"
#include "src/Hamiltonian/Hamiltonian.hpp"
#include "src/Lindblad-PE/lindblad.hpp"
#include "src/hdf5io/hdf5io.hpp"

#ifdef MKL
  #include "mkl.h"
#endif

#ifndef DEBUG
#define DEBUG 1
#endif

#ifndef NumCores
#define NumCores 1
#endif

#define FIXJ13 0

const int L = 3;
const int N1 = L;// Open system has no upper limit
const int N2 = L;
const RealType t12 = 1.0e0;

void LoadParameters( const std::string filename, RealType &t23,  RealType &t13,
  std::vector<RealType> &Uloc, std::vector<RealType> &Vloc,
  RealType &gammaL, RealType &gammaR, int &Tsteps, RealType &dt){
    HDF5IO file(filename);
    t23 = t12 * file.loadReal("Parameters", "t23");
    t13 = t12 * file.loadReal("Parameters", "t13");
    file.loadStdVector("Parameters", "U", Uloc);
    file.loadStdVector("Parameters", "V", Vloc);
    gammaL = file.loadReal("Parameters", "gammaL");
    gammaR = file.loadReal("Parameters", "gammaR");
    Tsteps = file.loadInt("Parameters", "Tsteps");
    dt = file.loadReal("Parameters", "dt");
}

RealType TraceRhos(const std::vector<ComplexMatrixType> &Rhos){
  ComplexType tr = ComplexType(0.0e0, 0.0e0);
  for ( auto &rho: Rhos){
    tr += rho.trace();
  }
  assert( std::abs(tr.imag()) < 1.0e-12 );
  return tr.real();
}

ComplexMatrixType SingleParticleDensityMatrix( const int spin, const std::vector<std::vector<Basis> > &Bases,
  const std::vector<ComplexMatrixType> &Rhos, std::vector<Hamiltonian<ComplexType> > &ham ){
  size_t L = Bases.at(0).at(spin).getL();
  ComplexMatrixType CM = ComplexMatrixType::Zero(L,L);
  for ( size_t cntB = 0; cntB < ham.size(); cntB++){
    std::vector<int> bs = Bases.at(cntB).at(spin).getFStates();
    std::vector<size_t> tg = Bases.at(cntB).at(spin).getFTags();
    ComplexMatrixType Rho = Rhos.at(cntB);
    for ( const int &b : bs ){
      size_t bid = tg.at(b);// Find their indices
      for ( size_t i=0; i < L; i++){
        for ( size_t j=i; j < L; j++){
          /* see if hopping exist */
          if ( (btest(b, j) && !(btest(b, i))) || (i == j && btest(b, i)) ) {
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
            int p;
            if ( i == j ){
              p = b;
            }else{
              p = ibset(b, i);
              p = ibclr(p, j);
            }
            size_t pid = tg.at(p);// Find their indices
            size_t count;
            if ( spin == 0 ) count = Bases.at(cntB).at(1).getHilbertSpace();
            else if ( spin == 1 ) count = Bases.at(cntB).at(0).getHilbertSpace();
            std::vector<size_t> rids(2, bid);
            std::vector<size_t> cids(2, pid);
            for (size_t loop_id = 0; loop_id < count; loop_id++) {
              if ( spin == 0 ){
                rids.at(1) = loop_id;
                cids.at(1) = loop_id;
              }else if ( spin == 1 ){
                rids.at(0) = loop_id;
                cids.at(0) = loop_id;
              }
              size_t rid = ham.at(cntB).DetermineTotalIndex( rids );
              size_t cid = ham.at(cntB).DetermineTotalIndex( cids );
              CM(i, j) +=  tsign * Rho(rid, cid);
            }
          }
        }
      }
    }
  }
  return CM;
}

ComplexMatrixType NiNj( const int spin, const std::vector<std::vector<Basis> > &Bases,
  const std::vector<ComplexMatrixType> &Rhos ){
  ComplexMatrixType tmp = ComplexMatrixType::Zero(Bases.at(0).at(spin).getL(), Bases.at(0).at(spin).getL());
  for (size_t cnt = 0; cnt < Bases.size(); cnt++) {
    std::vector<int> bs = Bases.at(cnt).at(spin).getFStates();
    size_t coff = 0;
    for ( auto &nf : bs ){
      for (size_t cnt1 = 0; cnt1 < Bases.at(cnt).at(spin).getL(); cnt1++) {
        for (size_t cnt2 = 0; cnt2 < Bases.at(cnt).at(spin).getL(); cnt2++) {
          if ( btest(nf, cnt1) && btest(nf, cnt2) ){
            tmp(cnt1, cnt2) += Rhos.at(cnt)(coff, coff);
          }
        }
      }
      coff++;
    }
  }
  return tmp;
}

int main(int argc, char const *argv[]) {
  Eigen::setNbThreads(NumCores);
#ifdef MKL
  mkl_set_num_threads(NumCores);
#endif
  INFO("Eigen3 uses " << Eigen::nbThreads() << " threads.");
  RealType t23, t13;
  int Tsteps;
  RealType dt, gammaL, gammaR;
  std::vector<RealType> Uin, Vin;
  LoadParameters( "conf.h5", t23, t13, Uin, Vin, gammaL, gammaR, Tsteps, dt);

  INFO("Build 3-site Triangle - ");
  std::vector< Node<ComplexType>* > LOOP;
  Node<ComplexType> *A = new Node<ComplexType>(0, NULL, t12, "A");
  LOOP.push_back(A);
  Node<ComplexType> *B = new Node<ComplexType>(1, LOOP[0], t12);
  LOOP.push_back(B);
  Node<ComplexType> *C = new Node<ComplexType>(2, LOOP[1], t23);
  LOOP.push_back(C);
  LOOP[0]->LinkTo( LOOP[2], t13);
  assert( LOOP.size() == L );
  for ( auto &lt : LOOP ){
    if ( !(lt->VerifySite()) ) RUNTIME_ERROR("Wrong lattice setup!");
  }
  INFO("DONE!");

  INFO("Build Basis in each U(1) sector - ");
  std::map<std::pair<int,int>, int > PairIndex1;
  std::map<int, std::pair<int,int> > PairIndex2;
  std::vector<std::vector<Basis> > Bases;
  size_t cntIdx = 0;
  for (size_t cntN1 = 0; cntN1 <= N1; cntN1++) {
    Basis F1(L, cntN1, true);
    F1.Fermion();
    for (size_t cntN2 = 0; cntN2 <= N2; cntN2++) {
      Basis F2(L, cntN2, true);
      // Basis F2(L, cntN1, true);
      F2.Fermion();
      std::vector<Basis> tmp;
      tmp.push_back(F1);
      tmp.push_back(F2);
      Bases.push_back(tmp);
      PairIndex1[std::make_pair(cntN1,cntN2)] = cntIdx;
      // std::cout << PairIndex1.at(std::make_pair(cntN1,cntN2)) << std::endl;
      PairIndex2[cntIdx] = std::make_pair(cntN1,cntN2);
      // std::cout << PairIndex2.at(cntIdx).first << " " << PairIndex2.at(cntIdx).second << std::endl;
      cntIdx += 1;
    }
  }
  INFO("DONE!");

  INFO("Build Hamiltonian in each Basis pair - ");
  std::vector<Hamiltonian<ComplexType> > Hams;
  std::vector< std::vector<ComplexType> > Vloc;
  std::vector<ComplexType> vtmp;
  for ( auto &val: Vin ){
    vtmp.push_back(ComplexType(val));
  }
  Vloc.push_back(vtmp);
  Vloc.push_back(vtmp);
  std::vector< std::vector<ComplexType> > Uloc;
  std::vector<ComplexType> utmp;
  for ( auto &val: Uin ){
    utmp.push_back(ComplexType(val));
  }
  Uloc.push_back(utmp);
  Uloc.push_back(utmp);
  // cntIdx = 0;
  for ( auto &bs : Bases ){
    // std::cout <<
    // PairIndex2.at(cntIdx).first << " " <<
    // PairIndex2.at(cntIdx).second << std::endl;
    Hamiltonian<ComplexType> ham( bs );
    ham.BuildLocalHamiltonian(Vloc, Uloc, bs);
    // std::cout << "1" << std::endl;
    ham.BuildHoppingHamiltonian(bs, LOOP);
    // std::cout << "2" << std::endl;
    ham.BuildTotalHamiltonian();
    // std::cout << "3" << std::endl;
    Hams.push_back(ham);
    // cntIdx += 1;
  }
  INFO("DONE!");

  INFO("Build Initial Density Matrix - ");
  std::vector<ComplexMatrixType> Rhos;
  int cnt = 0;
  ComplexMatrixType Rho;
  for ( auto &bs : Bases ){
    size_t hb = bs.at(0).getHilbertSpace() * bs.at(1).getHilbertSpace();
    std::cout << cnt << " " << bs.at(0).getHilbertSpace() << " " << bs.at(1).getHilbertSpace() << std::endl;
    if ( cnt == 0 ){
      Rho = ComplexMatrixType::Identity(hb, hb);
    }else{
      Rho = ComplexMatrixType::Zero(hb, hb);
    }
    Rhos.push_back(Rho);
    cnt += 1;
  }
  INFO("DONE!");

  std::cout << "Establish the index." << std::endl;
  std::vector<RealType> Gammas;
  std::vector<std::tuple<int,int,int> > SiteTypesSpin;
  std::vector<std::vector<std::pair<int,int> > > BasisIds;
  std::vector<std::vector<std::vector<std::pair<size_t, size_t> > > > CollapseIds;
  std::vector<std::pair<int,int> > w1;
  std::vector<std::vector<std::pair<size_t, size_t> > > w2;
  w1.clear();
  w2.clear();
  std::cout << "cfdag1" << std::endl;
  Cfdagger(0, 0, Bases, Hams, PairIndex1, PairIndex2, w1, w2);
  Gammas.push_back(gammaL);
  SiteTypesSpin.push_back(std::make_tuple(0, 1, 0));
  BasisIds.push_back(w1);
  CollapseIds.push_back(w2);
  w1.clear();
  w2.clear();
  std::cout << "cfdag2" << std::endl;
  Cfdagger(0, 1, Bases, Hams, PairIndex1, PairIndex2, w1, w2);
  Gammas.push_back(gammaL);
  SiteTypesSpin.push_back(std::make_tuple(0, 1, 1));
  BasisIds.push_back(w1);
  CollapseIds.push_back(w2);
  w1.clear();
  w2.clear();
  std::cout << "cf1" << std::endl;
  Cf(2, 0, Bases, Hams, PairIndex1, PairIndex2, w1, w2);
  Gammas.push_back(gammaR);
  SiteTypesSpin.push_back(std::make_tuple(2,-1, 0));
  BasisIds.push_back(w1);
  CollapseIds.push_back(w2);
  w1.clear();
  w2.clear();
  std::cout << "cf2" << std::endl;
  Cf(2, 1, Bases, Hams, PairIndex1, PairIndex2, w1, w2);
  Gammas.push_back(gammaR);
  SiteTypesSpin.push_back(std::make_tuple(2,-1, 1));
  BasisIds.push_back(w1);
  CollapseIds.push_back(w2);
  w1.clear();
  w2.clear();
  std::cout << "DONE!" << std::endl;

  std::cout << "Dynamics begins..." << std::endl;
  size_t MaxTry = 0;
#if FIXJ13
  RealType FinalJ13 = 1.0e10, PrevJ13 = 0.0e0;
  const RealType Target13 = 0.10e0;
  const RealType J13Tol = 1.0e-10;
  while ( (std::abs(TargetJ13 - FinalJ13) > J13Tol) && MaxTry < 1000 ){
#endif
  std::vector<RealType> tls;
  std::vector<ComplexType> Na, Nb, Nc, n12, n13, n23, NaUp;
  std::vector<ComplexType> j12, j23, j13;
  tls.push_back(0.0e0);
  ComplexMatrixType CM0 = SingleParticleDensityMatrix( 0, Bases, Rhos, Hams);
  ComplexMatrixType CM1 = SingleParticleDensityMatrix( 1, Bases, Rhos, Hams);
  ComplexMatrixType Nij0 = NiNj( 0, Bases, Rhos );
  ComplexMatrixType Nij1 = NiNj( 1, Bases, Rhos );
  NaUp.push_back(CM0(0,0));
  Na.push_back(CM0(0,0) + CM1(0,0));
  Nb.push_back(CM0(1,1) + CM1(1,1));
  Nc.push_back(CM0(2,2) + CM1(2,2));
  n12.push_back(Nij0(0,1) + Nij1(0,1));
  n23.push_back(Nij0(1,2) + Nij1(1,2));
  n13.push_back(Nij0(0,2) + Nij1(0,2));
  j12.push_back( t12 * (CM0(0,1) + CM1(0,1)) );
  j23.push_back( t23 * (CM0(1,2) + CM1(1,2)) );
  j13.push_back( t13 * (CM0(0,2) + CM1(0,2)) );
  INFO("Trace = " << TraceRhos(Rhos));
  for (size_t cntT = 1; cntT <= Tsteps; cntT++) {
    FRK4( dt, Gammas, SiteTypesSpin, BasisIds, CollapseIds, Bases, Hams, Rhos);
    if ( cntT % 20 == 0 ){
      /* NOTE: Check trace of Rhos */
      RealType tr = TraceRhos(Rhos);
      INFO("Trace = " << tr);
      /* NOTE: Expectation values */
      CM0 = SingleParticleDensityMatrix( 0, Bases, Rhos, Hams);
      CM1 = SingleParticleDensityMatrix( 1, Bases, Rhos, Hams);
      Nij0 = NiNj( 0, Bases, Rhos );
      Nij1 = NiNj( 1, Bases, Rhos );
      NaUp.push_back(CM0(0,0));
      Na.push_back(CM0(0,0) + CM1(0,0));
      Nb.push_back(CM0(1,1) + CM1(1,1));
      Nc.push_back(CM0(2,2) + CM1(2,2));
      n12.push_back(Nij0(0,1) + Nij1(0,1));
      n23.push_back(Nij0(1,2) + Nij1(1,2));
      n13.push_back(Nij0(0,2) + Nij1(0,2));
      j12.push_back( t12 * (CM0(0,1) + CM1(0,1)) );
      j23.push_back( t23 * (CM0(1,2) + CM1(1,2)) );
      j13.push_back( t13 * (CM0(0,2) + CM1(0,2)) );
      tls.push_back(cntT * dt);
    }
  }
#if FIXJ13
  PrevJ13  = FinalJ13;
  FinalJ13 = j13.at(j13.size()-1);
  if ( PrevJ13 > FinalJ13 ){//Operates at small gamma's
  // if ( PrevJ13 < FinalJ13 ){//Operates at large gamma's
    Gammas.clear();
    gammaL *= 0.990e0;
    gammaR *= 0.990e0;
    Gammas.push_back(gammaL);
    Gammas.push_back(gammaL);
    Gammas.push_back(gammaR);
    Gammas.push_back(gammaR);
  }else{
    Gammas.clear();
    gammaL *= 1.010e0;
    gammaR *= 1.010e0;
    Gammas.push_back(gammaL);
    Gammas.push_back(gammaL);
    Gammas.push_back(gammaR);
    Gammas.push_back(gammaR);
  }
  MaxTry++;
  }
#endif
  if ( MaxTry < 1000 ){
  /* NOTE: H5 group name */
  HDF5IO file = HDF5IO("TriFermi.h5");
  file.saveNumber("Obs", "GammaL", gammaL);
  file.saveNumber("Obs", "GammaR", gammaR);
  file.saveNumber("Obs", "t12", t12);
  file.saveNumber("Obs", "t23", t23);
  file.saveNumber("Obs", "t13", t13);
  file.saveStdVector("Obs", "tls", tls);
  file.saveStdVector("Obs", "NaUp", NaUp);
  file.saveStdVector("Obs", "Na", Na);
  file.saveStdVector("Obs", "Nb", Nb);
  file.saveStdVector("Obs", "Nc", Nc);
  file.saveStdVector("Obs", "n12", n12);
  file.saveStdVector("Obs", "n23", n23);
  file.saveStdVector("Obs", "n13", n13);
  file.saveStdVector("Obs", "j12", j12);
  file.saveStdVector("Obs", "j23", j23);
  file.saveStdVector("Obs", "j13", j13);
  std::cout << "DONE!" << std::endl;
  }
  return 0;
}
