/* Triangular Quantum Dot Metastructure */
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include "src/EDType.hpp"
#include "src/bitwise.h"
#include "src/Node/Node.hpp"
#include "src/Basis/Basis.hpp"
#include "src/ArmadilloMatrix.hpp"
#include "src/Hamiltonian/Hamiltonian.hpp"
#include "src/hdf5io/hdf5io.hpp"
#include "src/Lindblad/ParticleExchange/lindblad.hpp"

#ifdef MKL
  #include "mkl.h"
#endif

const int L = 3;
const int N1 = L;// Open system has no upper limit
const int N2 = L;
const RealType t12 = 1.0e0;

void LoadParameters( const std::string filename, RealType &t23,  RealType &t13, std::vector<RealType> &Uloc, std::vector<RealType> &Vloc, RealType &GammaL, RealType &GammaR, RealType &dt){
  HDF5IO file(filename);
  file.LoadNumber("Parameters", "t23", t23);
  t23 *= t12;
  file.LoadNumber("Parameters", "t13", t13);
  t13 *= t12;
  file.LoadStdVector("Parameters", "U", Uloc);
  file.LoadStdVector("Parameters", "V", Vloc);
  file.LoadNumber("Parameters", "GammaL", GammaL);
  file.LoadNumber("Parameters", "GammaR", GammaR);
  file.LoadNumber("Parameters", "dt", dt);
}

RealType TraceRhos(const std::vector<ComplexMatrixType> &Rhos){
  ComplexType tr = ComplexType(0.0e0, 0.0e0);
  for ( auto &rho: Rhos){
    tr += arma::trace(rho);
  }
  assert( std::abs(tr.imag()) < 1.0e-12 );
  return tr.real();
}

ComplexMatrixType SingleParticleDensityMatrix( const int spin, const std::vector<std::vector<Basis> >& Bases, const std::vector<ComplexMatrixType>& Rhos, const std::vector<Hamiltonian<ComplexType> >& ham ){
  size_t L = Bases.at(0).at(spin).getL();
  ComplexMatrixType CM(L, L, arma::fill::zeros);
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

ComplexMatrixType NiNj( const int spin, const std::vector<std::vector<Basis> > &Bases, const std::vector<ComplexMatrixType> &Rhos ){
  ComplexMatrixType tmp(Bases.at(0).at(spin).getL(), Bases.at(0).at(spin).getL(), arma::fill::zeros);
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

std::vector<ComplexMatrixType> SteadyState(const std::vector<std::vector<Basis> >& Bases, const std::vector<ComplexMatrixType>& OriginalRhos, const std::vector<Hamiltonian<ComplexType> >& Hams, const RealType& dt, const std::vector<RealType>& Gammas, const std::vector<std::tuple<int,int,int> >& SiteTypesSpin, const std::vector<std::vector<std::pair<int,int> > >& BasisIds, const std::vector<std::vector<std::vector<std::pair<size_t, size_t> > > >& CollapseIds, const int& Save = 0, const std::string prefix = "" ){
  std::vector<RealType> tls;
  std::vector<ComplexType> Na, Nb, Nc, n12, n13, n23, NaUp;
  std::vector<ComplexType> j12, j23, j13;
  std::vector<ComplexMatrixType> Rhos = OriginalRhos;
  tls.clear();
  NaUp.clear();
  Na.clear();
  Nb.clear();
  Nc.clear();
  n12.clear();
  n13.clear();
  n23.clear();
  j12.clear();
  j23.clear();
  j13.clear();
  ComplexMatrixType CM0 = SingleParticleDensityMatrix( 0, Bases, Rhos, Hams);
  ComplexMatrixType CM1 = SingleParticleDensityMatrix( 1, Bases, Rhos, Hams);
  ComplexMatrixType Nij0 = NiNj( 0, Bases, Rhos );
  ComplexMatrixType Nij1 = NiNj( 1, Bases, Rhos );
  if ( Save ){
    tls.push_back(0.0e0);
    NaUp.push_back(CM0(0,0));
    Na.push_back( CM0(0,0) + CM1(0,0) );
    Nb.push_back( CM0(1,1) + CM1(1,1) );
    Nc.push_back( CM0(2,2) + CM1(2,2) );
    n12.push_back( Nij0(0,1) + Nij1(0,1) );
    n23.push_back( Nij0(1,2) + Nij1(1,2) );
    n13.push_back( Nij0(0,2) + Nij1(0,2) );
    j12.push_back( CM0(0,1) + CM1(0,1) );
    j23.push_back( CM0(1,2) + CM1(1,2) );
    j13.push_back( CM0(0,2) + CM1(0,2) );
  }
  bool Converged = false;
  int cntT = 1;
  while ( !Converged && cntT < 200000 ){
    FRK4( dt, Gammas, SiteTypesSpin, BasisIds, CollapseIds, Bases, Hams, Rhos );
    if ( cntT % 20 == 0 ){
      CM0 = SingleParticleDensityMatrix( 0, Bases, Rhos, Hams);
      CM1 = SingleParticleDensityMatrix( 1, Bases, Rhos, Hams);
      if ( (std::abs(Na.back() - (CM0(0,0) + CM1(0,0))) < 1.0e-8) && (std::abs(Nb.back() - (CM0(1,1) + CM1(1,1))) < 1.0e-8) ) Converged = true;
      if ( Save ){
        Nij0 = NiNj( 0, Bases, Rhos );
        Nij1 = NiNj( 1, Bases, Rhos );
        tls.push_back(cntT * dt);
        NaUp.push_back(CM0(0,0));
        Na.push_back( CM0(0,0) + CM1(0,0) );
        Nb.push_back( CM0(1,1) + CM1(1,1) );
        Nc.push_back( CM0(2,2) + CM1(2,2) );
        n12.push_back( Nij0(0,1) + Nij1(0,1) );
        n23.push_back( Nij0(1,2) + Nij1(1,2) );
        n13.push_back( Nij0(0,2) + Nij1(0,2) );
        j12.push_back( CM0(0,1) + CM1(0,1) );
        j23.push_back( CM0(1,2) + CM1(1,2) );
        j13.push_back( CM0(0,2) + CM1(0,2) );
      }
    }
    cntT++;
  }
  /* NOTE: H5 group name */
  if ( Save ){
    HDF5IO* file = new HDF5IO(prefix + "TQDM.h5");
    file->SaveStdVector("Obs", "tls", tls);
    file->SaveStdVector("Obs", "NaUp", NaUp);
    file->SaveStdVector("Obs", "Na", Na);
    file->SaveStdVector("Obs", "Nb", Nb);
    file->SaveStdVector("Obs", "Nc", Nc);
    file->SaveStdVector("Obs", "n12", n12);
    file->SaveStdVector("Obs", "n23", n23);
    file->SaveStdVector("Obs", "n13", n13);
    file->SaveStdVector("Obs", "j12", j12);
    file->SaveStdVector("Obs", "j23", j23);
    file->SaveStdVector("Obs", "j13", j13);
    delete file;
  }
  return Rhos;
}

void Dynamics( const std::string prefix, const int SearchJ = 0 ){
  std::ofstream LogOut;
  LogOut.open(prefix + "TQDM.SS", std::ios::app);
  RealType t23, t13;
  RealType dt, GammaL, GammaR;
  std::vector<RealType> Uin, Vin;
  try{
    /* Load parameters from file */
    H5::Exception::dontPrint();
    H5::H5File::isHdf5(prefix + "conf.h5");
    LoadParameters( prefix + "conf.h5", t23, t13, Uin, Vin, GammaL, GammaR, dt);
  }catch(H5::FileIException){
    t23 = t12;
    t13 = t12;
    Uin = std::vector<RealType>(L, 0.0);
    Vin = std::vector<RealType>(L, 0.0);
    dt = 0.005;
  }
  LogOut << "Build 3-site Triangle Lattice " << std::flush;
  std::vector< Node<ComplexType>* > LOOP;
  Node<ComplexType> *A = new Node<ComplexType>(0, NULL, t12, "A");
  LOOP.push_back(A);
  Node<ComplexType> *B = new Node<ComplexType>(1, LOOP[0], t12);
  LOOP.push_back(B);
  Node<ComplexType> *C = new Node<ComplexType>(2, LOOP[1], t23);
  LOOP.push_back(C);
  LOOP[0]->LinkTo( LOOP[2], t13);
  assert( LOOP.size() == L );
  LogOut << "DONE!" << std::endl;

  LogOut << "Build Basis in each U(1) sector - " << std::flush;
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
  LogOut << "DONE!" << std::endl;

  LogOut << "Build Hamiltonian in each Basis pair - " << std::flush;
  std::vector<Hamiltonian<ComplexType> > Hams;
  std::vector<ComplexType> vtmp(Vin.begin(), Vin.end());
  std::vector<std::vector<ComplexType> > Vloc = vec(vtmp, vtmp);
  std::vector<ComplexType> Uloc(Uin.begin(), Uin.end());
  for ( auto &bs : Bases ){
    Hamiltonian<ComplexType> ham( bs );
    ham.FermiHubbardModel(bs, LOOP, Vloc, Uloc);
    Hams.push_back(ham);
  }
  LogOut << "DONE!" << std::endl;

  LogOut << "Build Initial Density Matrix - " << std::endl;
  std::vector<ComplexMatrixType> OriginalRhos;
  int cnt = 0;
  ComplexMatrixType Rho;
  for ( auto &bs : Bases ){
    size_t hb = bs.at(0).getHilbertSpace() * bs.at(1).getHilbertSpace();
    LogOut << "\t" << cnt << " " << bs.at(0).getHilbertSpace() << " " << bs.at(1).getHilbertSpace() << std::endl;
    if ( cnt == 0 ){
      Rho = ComplexMatrixType(hb, hb, arma::fill::eye);
    }else{
      Rho = ComplexMatrixType(hb, hb, arma::fill::zeros);
    }
    OriginalRhos.push_back(Rho);
    cnt += 1;
  }
  LogOut << "DONE!" << std::endl;

  LogOut << "Establish the index for Lindblad equation." << std::endl;
  std::vector<RealType> Gammas;
  std::vector<std::tuple<int,int,int> > SiteTypesSpin;
  std::vector<std::vector<std::pair<int,int> > > BasisIds;
  std::vector<std::vector<std::vector<std::pair<size_t, size_t> > > > CollapseIds;
  std::vector<std::pair<int,int> > w1;
  std::vector<std::vector<std::pair<size_t, size_t> > > w2;
  w1.clear();
  w2.clear();
  LogOut << "C^dagger_0,up" << std::endl;
  Cfdagger(0, 0, Bases, Hams, PairIndex1, PairIndex2, w1, w2);
  Gammas.push_back(GammaL);
  SiteTypesSpin.push_back(std::make_tuple(0, 1, 0));
  BasisIds.push_back(w1);
  CollapseIds.push_back(w2);
  w1.clear();
  w2.clear();
  std::cout << "C^dagger_0,dn" << std::endl;
  Cfdagger(0, 1, Bases, Hams, PairIndex1, PairIndex2, w1, w2);
  Gammas.push_back(GammaL);
  SiteTypesSpin.push_back(std::make_tuple(0, 1, 1));
  BasisIds.push_back(w1);
  CollapseIds.push_back(w2);
  w1.clear();
  w2.clear();
  std::cout << "C_2,up" << std::endl;
  Cf(2, 0, Bases, Hams, PairIndex1, PairIndex2, w1, w2);
  Gammas.push_back(GammaR);
  SiteTypesSpin.push_back(std::make_tuple(2,-1, 0));
  BasisIds.push_back(w1);
  CollapseIds.push_back(w2);
  w1.clear();
  w2.clear();
  std::cout << "C_2,dn" << std::endl;
  Cf(2, 1, Bases, Hams, PairIndex1, PairIndex2, w1, w2);
  Gammas.push_back(GammaR);
  SiteTypesSpin.push_back(std::make_tuple(2,-1, 1));
  BasisIds.push_back(w1);
  CollapseIds.push_back(w2);
  LogOut << "DONE!" << std::endl;

  LogOut << "Dynamics begins..." << std::endl;
  int Save = 1;
  LogOut << "Begin dynamics, trace = " << TraceRhos(OriginalRhos) << std::endl;
  std::vector<ComplexMatrixType>  FinalRhos = SteadyState( Bases, OriginalRhos, Hams, dt, Gammas, SiteTypesSpin, BasisIds, CollapseIds, Save, prefix );
  LogOut << "After dynamics, trace = " << TraceRhos(FinalRhos) << std::endl;
}

/* main program */
int main(int argc, char *argv[]){
  if ( argc < 2 ) RUNTIME_ERROR(" Need at least one argument to run program. Use 0 to run steady state current.");
#if defined(MKL)
  mkl_set_num_threads(NumCores);
#endif
  int SearchJ = std::atoi(argv[1]);
  Dynamics( "", SearchJ );
  return 0;
}
