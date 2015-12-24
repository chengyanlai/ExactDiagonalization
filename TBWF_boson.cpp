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

#ifdef MKL
  #include "mkl.h"
#endif

#ifndef DEBUG
#define DEBUG 1
#endif

#ifndef NumCores
#define NumCores 2
#endif

#ifndef METHOD
#define METHOD 1
#endif

void LoadParameters( const std::string filename, int &L, int &OBC, int &N,
  RealType &Uloc, RealType &Vloc, RealType &dt, int &Tsteps, int &TBloc,
  RealType &factor);
void SaveObs( const std::string filename, const std::string gname,
  const std::vector<ComplexType> &v, const ComplexMatrixType &m);
void TerminatorBeam( const size_t TBloc, const Basis &bs, ComplexVectorType &Vec,
  const bool HARD_CUT = true, const RealType factor = 0.0e0 );
void GetGS1( const size_t TBloc, const Basis &bs, ComplexVectorType &Vec,
  const bool HARD_CUT = true, const RealType factor = 0.0e0 );
void GetGS2( const size_t TBloc, const Basis &bs, ComplexVectorType &Vec );
std::vector<ComplexType> Ni( const std::vector<Basis> &Bases, const ComplexVectorType &Vec );
ComplexMatrixType NiNj( const std::vector<Basis> &Bases, const ComplexVectorType &Vec );

int main(int argc, char const *argv[]) {
  Eigen::setNbThreads(NumCores);
#ifdef MKL
  mkl_set_num_threads(NumCores);
#endif
  INFO("Eigen3 uses " << Eigen::nbThreads() << " threads.");
  int L, OBC, N, Tsteps, TBloc;
  RealType Uin, Vin, dt, factor;
  LoadParameters( "conf.h5", L, OBC, N, Uin, Vin, dt, Tsteps, TBloc, factor );
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
  bool HARD_CUT = true;
  assert( factor >= 0.0e0 && factor < 1.0e0 );
  if ( factor > 1.0e-10 ) {
    HARD_CUT = false;
  }
  B1.BosonTB(TBloc, HARD_CUT);
  // std::vector< std::vector<int> > st = B1.getBStates();
  // std::vector< RealType > tg = B1.getBTags();
  // for (size_t cnt = 0; cnt < tg.size(); cnt++) {
  //   INFO_NONEWLINE( std::setw(3) << cnt << " - ");
  //   for (auto &j : st.at(cnt)){
  //     INFO_NONEWLINE(j << " ");
  //   }
  //   INFO("- " << tg.at(cnt));
  // }
  INFO("DONE!");
  INFO("Build GS - ");
  ComplexVectorType Vec = ComplexVectorType::Zero(B1.getHilbertSpace());
  if ( METHOD == 1 ) {
    GetGS1(TBloc, B1, Vec);
  } else if ( METHOD == 2 ) {
    GetGS2(TBloc, B1, Vec);
  }
  INFO("DONE!");
  INFO("Build Hamiltonian - ");
  std::vector<Basis> Bases;
  Bases.push_back(B1);
  Hamiltonian<ComplexType> ham( Bases );
  std::vector< std::vector<ComplexType> > Vloc;
  std::vector<ComplexType> Vtmp(L, (ComplexType)Vin);
  Vloc.push_back(Vtmp);
  std::vector< std::vector<ComplexType> > Uloc;
  std::vector<ComplexType> Utmp(L, (ComplexType)Uin);
  Uloc.push_back(Utmp);
  ham.BuildLocalHamiltonian(Vloc, Uloc, Bases);
  ham.BuildHoppingHamiltonian(Bases, lattice);
  ham.BuildTotalHamiltonian();
  INFO("DONE!");
  std::vector<ComplexType> Nbi = Ni( Bases, Vec );
  ComplexMatrixType Nij = NiNj( Bases, Vec );
  SaveObs("TBWF.h5", "Obs-0", Nbi, Nij);
  for (size_t cntT = 1; cntT <= Tsteps; cntT++) {
    ComplexType Prefactor = ComplexType(0.0, -1.0e0*dt);/* NOTE: hbar = 1 */
    ham.expH(Prefactor, Vec);
    TerminatorBeam(TBloc, B1, Vec);
    // INFO(" ");
    // INFO(Vec.norm());
    // std::cin.get();
    /* NOTE: Expectation values */
    Nbi = Ni( Bases, Vec );
    Nij = NiNj( Bases, Vec );
    /* NOTE: H5 group name */
    std::string gname = "Obs-";
    gname.append(std::to_string((unsigned long long)cntT));
    gname.append("/");
    SaveObs("TBWF.h5", gname, Nbi, Nij);
  }
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

void TerminatorBeam( const size_t TBloc, const Basis &bs,
  ComplexVectorType &Vec, const bool HARD_CUT, const RealType factor ){
  assert( Vec.size() == bs.getHilbertSpace() );
  // ComplexVectorType TBVec = Vec;
  std::map<RealType, std::vector<int> > BsMap;
  std::vector< std::vector<int> > bSt = bs.getBStates();
  std::vector<RealType> bTg = bs.getBTags();
  std::transform( bTg.begin(), bTg.end(), bSt.begin(),
    std::inserter(BsMap, BsMap.end() ), std::make_pair<RealType const&,std::vector<int> const&> );
  for (std::map<RealType, std::vector<int> >::iterator it = BsMap.begin(); it != BsMap.end(); it++) {
    if ( it->second.at(TBloc) > 0 ) {
      std::vector<int> newbs = it->second;
      if ( HARD_CUT || newbs.at(TBloc) == 1 ) {
        newbs.at(TBloc) = 0;
      } else if ( newbs.at(TBloc) > 0 ) {
        newbs.at(TBloc) -= 1;
      } else {
        std::cout << "Some Wrong in TerminatorBeam" << std::endl;
      }
      RealType newtag = BosonBasisTag(newbs);
      size_t new_idx = bs.getIndexFromTag(newtag);
      assert( new_idx < bs.getHilbertSpace() );
      size_t old_idx = bs.getIndexFromTag(it->first);
      assert( old_idx < bs.getHilbertSpace() );
      if ( HARD_CUT ) {
        RealType val1 = (std::conj(Vec(old_idx)) * Vec(old_idx)).real();
        RealType val2 = (std::conj(Vec(new_idx)) * Vec(new_idx)).real();
        Vec(new_idx) = (ComplexType)std::sqrt( val1 + val2 );
        Vec(old_idx) = ComplexType(0.0e0, 0.0e0);
      } else {
        RealType val1 = (std::conj(Vec(old_idx)) * Vec(old_idx)).real();
        RealType val2 = (std::conj(Vec(new_idx)) * Vec(new_idx)).real();
        Vec(new_idx) = (ComplexType)std::sqrt( std::sqrt(factor) * val1 + val2 );
        Vec(old_idx) = (ComplexType)std::sqrt( std::sqrt(1.0e0 - factor) * val1 );
      }
    }
  }
}

void GetGS1( const size_t TBloc, const Basis &bs, ComplexVectorType &Vec,
  const bool HARD_CUT, const RealType factor ){
  HDF5IO gsf("BSSH.h5");
  ComplexVectorType gswf;
  gsf.loadVector("GS", "EVec", gswf);
  Basis GS(bs.getL(), bs.getN());
  GS.Boson();
  if ( DEBUG ){
    int oldL = gsf.loadUlong("1DChain", "L");
    int oldN = gsf.loadUlong("1DChain", "N");
    assert( oldL == bs.getL() );
    assert( oldN == bs.getN() );
  }
  std::vector< std::vector<int> > gsbs = GS.getBStates();
  std::vector<RealType> gsts = GS.getBTags();
  std::map<RealType, std::vector<int> > GSmap;
  std::transform( gsts.begin(), gsts.end(), gsbs.begin(),
    std::inserter(GSmap, GSmap.end() ), std::make_pair<RealType const&,std::vector<int> const&> );
  for (std::map<RealType, std::vector<int> >::iterator it = GSmap.begin(); it != GSmap.end(); it++) {
    if ( it->second.at(TBloc) > 0 ) {
      std::vector<int> newbs = it->second;
      if ( HARD_CUT ) {
        newbs.at(TBloc) = 0;
      } else {
        newbs.at(TBloc) -= 1;
      }
      RealType newtag = BosonBasisTag(newbs);
      if ( HARD_CUT ) {
        Vec(bs.getIndexFromTag(newtag)) = gswf(GS.getIndexFromTag(it->first));
      } else {
        size_t new_idx = bs.getIndexFromTag(newtag);
        assert( new_idx < bs.getHilbertSpace() );
        size_t old_idx = bs.getIndexFromTag(it->first);
        assert( old_idx < bs.getHilbertSpace() );
        RealType val1 = (std::conj(Vec(old_idx)) * Vec(old_idx)).real();
        RealType val2 = (std::conj(Vec(new_idx)) * Vec(new_idx)).real();
        RealType val3 = (std::conj(gswf(GS.getIndexFromTag(it->first))) *
                         gswf(GS.getIndexFromTag(it->first)) ).real();
        Vec(bs.getIndexFromTag(newtag)) = (ComplexType)std::sqrt( std::sqrt(factor) * val3 + val2 );
        Vec(bs.getIndexFromTag(it->first)) = (ComplexType)std::sqrt( std::sqrt(1.0e0 - factor) * val3 + val1 );
      }
    }else {
      Vec(bs.getIndexFromTag(it->first)) = gswf(GS.getIndexFromTag(it->first));
    }
  }
}

void GetGS2( const size_t TBloc, const Basis &bs, ComplexVectorType &Vec )
{
  HDF5IO gsf("BSSH.h5");
  ComplexVectorType gswf;
  gsf.loadVector("GS", "EVec", gswf);
  Basis GS(bs.getL() - 1, bs.getN());
  GS.Boson();
  if ( DEBUG ){
    int oldL = gsf.loadUlong("1DChain", "L");
    int oldN = gsf.loadUlong("1DChain", "N");
    assert( oldL == bs.getL() - 1 );
    assert( oldN == bs.getN() );
  }
  std::vector< std::vector<int> > gsbs = GS.getBStates();
  std::vector<RealType> gsts = GS.getBTags();
  std::map<RealType, std::vector<int> > GSmap;
  std::transform( gsts.begin(), gsts.end(), gsbs.begin(),
    std::inserter(GSmap, GSmap.end() ), std::make_pair<RealType const&,std::vector<int> const&> );
  for (std::map<RealType, std::vector<int> >::iterator it = GSmap.begin(); it != GSmap.end(); it++) {
    std::vector<int> newbs = it->second;
    newbs.push_back(0);
    RealType newtag = BosonBasisTag(newbs);
    Vec(bs.getIndexFromTag(newtag)) = gswf(GS.getIndexFromTag(it->first));
  }
}

void LoadParameters( const std::string filename, int &L, int &OBC, int &N,
  RealType &Uloc, RealType &Vloc, RealType &dt, int &Tsteps, int &TBloc,
  RealType &factor){
    HDF5IO file(filename);
    L = file.loadInt("Parameters", "L");
    OBC = file.loadInt("Parameters", "OBC");
    N = file.loadInt("Parameters", "N");
    Uloc = file.loadReal("Parameters", "U");
    Vloc = file.loadReal("Parameters", "V");
    dt = file.loadReal("Parameters", "dt");
    Tsteps = file.loadInt("Parameters", "Tsteps");
    TBloc = file.loadInt("Parameters", "TBloc");
    factor = file.loadReal("Parameters", "FACTOR");
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
