#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include "src/EDType.hpp"
#include "src/bitwise.h"
#include "src/Node/Node.hpp"
#include "src/Lattice/preset.hpp"
#include "src/Basis/Basis.hpp"
#include "src/Hamiltonian/FHM/FermiHubbard.hpp"
#include "src/hdf5io/hdf5io.hpp"

#ifdef MKL
  #include "mkl.h"
#endif

#define DT ComplexType
#define DTV ComplexVectorType
#define DTM ComplexMatrixType
/* The real type define cause the compilation errors on dynamics functions */
// #define DT RealType
// #define DTV RealVectorType
// #define DTM RealMatrixType

// const RealType t = 1.0e0;

void LoadEqmParameters( const std::string filename, int &L, int &OBC, int &N1, int &N2, std::vector<RealType> &Jls, std::vector<RealType> &Uls, std::vector<RealType> &Vls, std::vector<RealType> &Wls){
  HDF5IO file(filename);
  file.LoadNumber("Parameters", "L", L);
  file.LoadNumber("Parameters", "OBC", OBC);
  file.LoadNumber("Parameters", "N1", N1);
  file.LoadNumber("Parameters", "N2", N2);
  file.LoadStdVector("Parameters", "J", Jls);
  file.LoadStdVector("Parameters", "U", Uls);
  file.LoadStdVector("Parameters", "V", Vls);
  file.LoadStdVector("Parameters", "W", Wls);
}

void LoadPumpParameters( const std::string filename, std::vector<RealType> &At, int& TSteps, RealType& dt){
  HDF5IO file(filename);
  file.LoadStdVector("Parameters", "At", At);
  file.LoadNumber("Parameters", "TStepsA", TSteps);
  file.LoadNumber("Parameters", "dtA", dt);
}

void LoadXASParameters( const std::string filename, std::vector<RealType> &Ufls, std::vector<RealType> &Vfls, int& TSteps, RealType& dt, int& CoreHole, int& Species, int& Type){
  HDF5IO file(filename);
  file.LoadStdVector("Parameters", "Uf", Ufls);
  file.LoadStdVector("Parameters", "Vf", Vfls);
  file.LoadNumber("Parameters", "TStepsX", TSteps);
  file.LoadNumber("Parameters", "dtX", dt);
  file.LoadNumber("Parameters", "CoreHole", CoreHole);
  file.LoadNumber("Parameters", "Species", Species);
  file.LoadNumber("Parameters", "Type", Type);
}

ComplexVectorType Operate( const ComplexVectorType& Vin, const int CoreHole, const int Species, const int Type, const std::vector<Basis>& OldBases, const std::vector<Basis>& NewBases, const Hamiltonian<ComplexType>& OldHam, const Hamiltonian<ComplexType>& NewHam  ){
  size_t NewHilbertSpace = NewBases.at(0).GetHilbertSpace() * NewBases.at(1).GetHilbertSpace();
  ComplexVectorType Vout(NewHilbertSpace, arma::fill::zeros);
  std::vector<int> OldFup = OldBases.at(0).GetFStates();
  std::vector<int> OldFdn = OldBases.at(1).GetFStates();
  std::vector<size_t> OldFupTag = OldBases.at(0).GetFTags();
  std::vector<size_t> OldFdnTag = OldBases.at(1).GetFTags();
  std::vector<size_t> NewFupTag = NewBases.at(0).GetFTags();
  std::vector<size_t> NewFdnTag = NewBases.at(1).GetFTags();
  size_t NewFupIdx, NewFdnIdx, OldFupIdx, OldFdnIdx;
  for ( auto OldFupState : OldFup ){
    RealType FermionSignUp = 1.0e0;
    if (Species == 1){
      NewFupIdx = NewFupTag.at(OldFupState);// Find their indices
      OldFupIdx = OldFupTag.at(OldFupState);// Find their indices
    }else if ( Species == 0 && Type == 1 && !(btest(OldFupState, CoreHole)) ){// c^\dagger_up
      int NewFupState = ibset(OldFupState, CoreHole);
      NewFupIdx = NewFupTag.at(NewFupState);// Find their indices
      OldFupIdx = OldFupTag.at(OldFupState);// Find their indices
      int NCross = 0;
      // for ( int i = 0; i < CoreHole; i++ ) NCross += btest(OldFupState, i);// Same as below
      for ( int i = CoreHole+1; i < OldBases.at(0).GetL(); i++ ) NCross += btest(OldFupState, i);
      if ( NCross % 2 == 1 ) FermionSignUp = -1.0e0;
    }else if ( Species == 0 && Type == -1 && btest(OldFupState, CoreHole) ){// c_up
      int NewFupState = ibclr(OldFupState, CoreHole);
      NewFupIdx = NewFupTag.at(NewFupState);// Find their indices
      OldFupIdx = OldFupTag.at(OldFupState);// Find their indices
    }else{
      continue;// next OldFupState
    }
    for ( auto OldFdnState : OldFdn ){
      // RealType FermionSignDn = 1.0;
      if ( Species == 0 ){
        NewFdnIdx = NewFdnTag.at(OldFdnState);// Find their indices
        OldFdnIdx = OldFdnTag.at(OldFdnState);// Find their indices
      }else if ( Species == 1 && Type == 1 && !(btest(OldFdnState, CoreHole)) ){// c^\dagger_down
        int NewFdnState = ibset(OldFdnState, CoreHole);
        NewFdnIdx = NewFdnTag.at(NewFdnState);// Find their indices
        OldFdnIdx = OldFdnTag.at(OldFdnState);// Find their indices
      }else if ( Species == 1 && Type == -1 && btest(OldFdnState, CoreHole) ){// c_down
        int NewFdnState = ibclr(OldFdnState, CoreHole);
        NewFdnIdx = NewFdnTag.at(NewFdnState);// Find their indices
        OldFdnIdx = OldFdnTag.at(OldFdnState);// Find their indices
      }else{
        continue;// next OldFdnState
      }
      std::vector<size_t> OldIdx, NewIdx;
      OldIdx.clear();
      OldIdx.push_back(OldFupIdx);
      OldIdx.push_back(OldFdnIdx);
      NewIdx.clear();
      NewIdx.push_back(NewFupIdx);
      NewIdx.push_back(NewFdnIdx);
      size_t oid = OldHam.DetermineTotalIndex( OldIdx );
      size_t nid = NewHam.DetermineTotalIndex( NewIdx );
      Vout(nid) = FermionSignUp * Vin(oid);
    }
  }
  return Vout;//arma::normalise(Vout);
}

DTM SingleParticleDensityMatrix( const int species, const std::vector<Basis> &Bases, const DTV &Vec, Hamiltonian<DT> &Ham0 ){
  size_t L = Bases.at(species).GetL();
  DTM CM(L, L, arma::fill::zeros);
  std::vector< int > bs = Bases.at(species).GetFStates();
  std::vector<size_t> tg = Bases.at(species).GetFTags();
  for ( const int &b : bs ){
    size_t bid = tg.at(b);// Find their indices
    for ( size_t i=0; i < L; i++){
      for ( size_t j=i; j < L; j++){
        /* see if hopping exist */
        if ( btest(b, j) && !(btest(b, i)) ) {
          /* c^\dagger_i c_j if yes, no particle in i and one particle in j. */
          int CrossFermionNumber = 0;
          DT tsign = (DT)(1.0e0);
          if ( j - i > 1 ){
            // possible cross fermions and sign change.
            for ( int k = i+1; k < j; k++){
              CrossFermionNumber += btest(b, k);
            }
          }
          if (CrossFermionNumber % 2 == 1) tsign = (DT)(-1.0e0);
          int p = ibset(b, i);
          p = ibclr(p, j);
          size_t pid = tg.at(p);// Find their indices
          size_t count;
          if ( species == 0 ) count = Bases.at(1).GetHilbertSpace();
          else if ( species == 1 ) count = Bases.at(0).GetHilbertSpace();
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
            size_t rid = Ham0.DetermineTotalIndex( rids );
            size_t cid = Ham0.DetermineTotalIndex( cids );
            CM(i, j) +=  tsign * Vec(cid) * Conjg( Vec(rid) );//Vec(id) * std::conj( Vec(id) );
          }
        }
      }
    }
  }
  return CM;
}

std::vector<DTM> BOWCorrelation( const std::vector<Basis> &Bases, const DTV &Vec, FHM<DT> &Ham0 ){
  size_t L = Bases.at(0).GetL();
  std::vector<DTM> out;
  DTM UpUp(Bases.at(0).GetL(), Bases.at(1).GetL(), arma::fill::zeros );
  DTM UpDn(Bases.at(0).GetL(), Bases.at(1).GetL(), arma::fill::zeros );
  DTM DnUp(Bases.at(0).GetL(), Bases.at(1).GetL(), arma::fill::zeros );
  DTM DnDn(Bases.at(0).GetL(), Bases.at(1).GetL(), arma::fill::zeros );
  std::vector<arma::SpMat<DT> > OBOW0 = Ham0.NNHoppingOp(0, Bases.at(0));
  std::vector<arma::SpMat<DT> > OBOW1 = Ham0.NNHoppingOp(1, Bases.at(1));
  for ( int iL = 0; iL < L; iL++ ){
    arma::SpMat<DT> Op1Up = OBOW0.at(iL);
    arma::SpMat<DT> Op1Dn = OBOW1.at(iL);
    std::cout << arma::approx_equal(Op1Up, Op1Up.t(), "absdiff", 1.0e-5) << std::endl;
    std::cout << arma::approx_equal(Op1Dn, Op1Dn.t(), "absdiff", 1.0e-5) << std::endl;
    for ( int jL = iL; jL < L; jL++ ){
      arma::SpMat<DT> Op2Up = OBOW0.at(jL);
      arma::SpMat<DT> Op2Dn = OBOW1.at(jL);
    std::cout << arma::approx_equal(Op2Up, Op2Up.t(), "absdiff", 1.0e-5) << std::endl;
    std::cout << arma::approx_equal(Op2Dn, Op2Dn.t(), "absdiff", 1.0e-5) << std::endl;

      // DTV tmp = Op1Up * ( Op2Up * Vec );
      DTV tmp = ( Op2Up * Vec );
      UpUp(iL, jL) = arma::cdot(Vec, tmp);

      // tmp = Op1Up * ( Op2Dn * Vec );
      tmp = ( Op2Dn * Vec );
      UpDn(iL, jL) = arma::cdot(Vec, tmp);

      // tmp = Op1Dn * ( Op2Up * Vec );
      tmp = ( Op2Up * Vec );
      DnUp(iL, jL) = arma::cdot(Vec, tmp);

      // tmp = Op1Dn * ( Op2Dn * Vec );
      tmp = ( Op2Dn * Vec );
      DnDn(iL, jL) = arma::cdot(Vec, tmp);
    }
  }
  out.push_back( UpUp );
  out.push_back( UpDn );
  out.push_back( DnUp );
  out.push_back( DnDn );
  return out;
}

std::vector<RealVectorType> Ni( const std::vector<Basis> &Bases, const DTV &Vec, Hamiltonian<DT> &Ham0 ){
  std::vector<RealVectorType> out;
  DTV tmp1(Bases.at(0).GetL(), arma::fill::zeros);//(Bases.at(0).GetL(), 0.0e0);
  DTV tmp2(Bases.at(1).GetL(), arma::fill::zeros);//(Bases.at(1).GetL(), 0.0e0);
  std::vector<int> f1 = Bases.at(0).GetFStates();
  std::vector<int> f2 = Bases.at(1).GetFStates();
  size_t f1id = 0, f2id = 0;
  for ( int &nf2 : f2 ){
    std::vector<size_t> ids(2,f2id);
    f1id = 0;
    for ( int &nf1 : f1 ){
      ids.at(0) = f1id;
      size_t id = Ham0.DetermineTotalIndex(ids);
      for (size_t cnt = 0; cnt < Bases.at(0).GetL(); cnt++) {
        if ( btest(nf1, cnt) ) tmp1(cnt) += Vec(id) * std::conj( Vec(id) );
        if ( btest(nf2, cnt) ) tmp2(cnt) += Vec(id) * std::conj( Vec(id) );
      }
      f1id++;
    }
    f2id++;
  }
  out.push_back( arma::real(tmp1) );
  out.push_back( arma::real(tmp2) );
  return out;
}

RealMatrixType NiNj( const std::vector<Basis> &Bases, const DTV &Vec, Hamiltonian<DT> &Ham0, const int species1, const int species2 ){
  DTM out(Bases.at(0).GetL(), Bases.at(1).GetL(), arma::fill::zeros );
  std::vector<int> f1 = Bases.at(0).GetFStates();
  std::vector<int> f2 = Bases.at(1).GetFStates();
  size_t f1id = 0, f2id = 0;
  for ( int &nf2 : f2 ){
    std::vector<size_t> ids(2,f2id);
    f1id = 0;
    for ( int &nf1 : f1 ){
      ids.at(0) = f1id;
      size_t id = Ham0.DetermineTotalIndex(ids);
      int Ns1 = (species1 == 0)? nf1 : nf2;
      int Ns2 = (species2 == 0)? nf1 : nf2;
      for (size_t cnt1 = 0; cnt1 < Bases.at(0).GetL(); cnt1++) {
        for (size_t cnt2 = 0; cnt2 < Bases.at(1).GetL(); cnt2++) {
          if ( btest(Ns1, cnt1) && btest(Ns2, cnt2) ) {
            out(cnt1, cnt2) +=  RealPart( Vec(id) * std::conj( Vec(id) ) );
          }
        }
      }
      f1id++;
    }
    f2id++;
  }
  return arma::real(out);
}

void SpectralPeaks( const double EG, const ComplexVectorType AS, const RealVectorType &EigVal, const ComplexMatrixType &EigVec,
  std::vector<double> &PeakLocation, std::vector<double> &PeakWeight, const int MaxNumPeak){
  PeakLocation.clear();
  PeakWeight.clear();
  for ( size_t i = 0; i < EigVec.n_cols; i++){
    PeakLocation.push_back(EigVal(i) - EG);
    ComplexVectorType An = EigVec.col(i);
    ComplexType val = arma::cdot(An, AS);
    PeakWeight.push_back( RealPart(Conjg(val) * val) );
    if ( i > MaxNumPeak ) break;
  }
  assert( PeakLocation.size() == PeakWeight.size() );
}

void Equilibrium(const std::string prefix){
  std::ofstream LogOut;
  LogOut.open(prefix + "FHMEqm.log", std::ios::app);
  int L;
  int OBC;
  int N1, N2;
  std::vector<RealType> Jin, Uin, Vin, Win;
  try{
    /* Load parameters from file */
    H5::Exception::dontPrint();
    H5::H5File::isHdf5("conf.h5");
    LoadEqmParameters( "conf.h5", L, OBC, N1, N2, Jin, Uin, Vin, Win);
  }catch(H5::FileIException){
    L = 12;
    OBC = 0;
    N1 = 6;
    N2 = 6;
    Jin = std::vector<RealType>(L, 1.0);// PBC
    Uin = std::vector<RealType>(L, 6.0);
    Vin = std::vector<RealType>(L, 0.0);
    Win = std::vector<RealType>(L, 0.0);
  }
  HDF5IO *file = new HDF5IO("FHMChainData.h5");
  LogOut << "Build Lattice - " << std::endl;
  std::vector<DT> J(Jin.begin(), Jin.end());
  std::vector< Node<DT>* > lattice = NN_1D_Chain(L, J, OBC);
  file->SaveNumber("1DChain", "L", L);
  file->SaveStdVector("1DChain", "J", J);
  for ( auto &lt : lattice ){
    if ( !(lt->VerifySite(LogOut)) ) RUNTIME_ERROR("Wrong lattice setup!");
  }
  LogOut << "DONE!" << std::endl;
  LogOut << "Build Basis - " << std::flush;
  // int N1 = (L+1)/2;
  Basis F1(L, N1, true);
  F1.Fermion();
  Basis F2(L, N2, true);
  F2.Fermion();
  std::vector<Basis> Bases;
  Bases.push_back(F1);
  Bases.push_back(F2);
  LogOut << "DONE!" << std::endl;
  LogOut << "Build Hamiltonian - " << std::flush;
  FHM<DT> Ham0( Bases );
  // Potential
  std::vector<DT> Vtmp(Vin.begin(), Vin.end());
  std::vector< std::vector<DT> > Vloc = vec(Vtmp, Vtmp);
  // Interaction
  std::vector<DT> Uloc(Uin.begin(), Uin.end());
  // Ham0.FermiHubbardModel(Bases, lattice, Vloc, Uloc);
  // Exntended Hubbard
  std::vector<DT> Wloc(Win.begin(), Win.end());
  Ham0.ExtendedFermiHubbardModel(Bases, lattice, Vloc, Uloc, Wloc);
  Ham0.CheckHermitian();
  LogOut << Ham0.GetTotalHilbertSpace() << " DONE!" << std::endl;
  LogOut << "Diagonalize Hamiltonian - " << std::flush;
  RealVectorType Vals;
  DTM Vecs;
  if ( Ham0.GetTotalHilbertSpace() > 1000 ){
    Ham0.eigh(Vals, Vecs, 2);
  }else{
    Ham0.diag(Vals, Vecs);// Full spectrum
  }
  LogOut << "DONE!" << std::endl;
  LogOut << "\tGS energy = " << Vals[0] << std::endl;
  LogOut << "\tFES energy = " << Vals[1] << std::endl;
  DTV Vec = Vecs.col(0);
  // ComplexSparseMatrixType H0 = Ham0.GetTotalHamiltonian();
  // std::cout << arma::cdot(Vec, H0 * Vec) << std::endl;
  file->SaveVector("GS", "EVec", Vec);
  file->SaveVector("GS", "EVal", Vals);
  std::vector<RealVectorType> Nfi = Ni( Bases, Vec, Ham0 );
  LogOut << " Up Spin - " << std::endl;
  LogOut << Nfi.at(0) << std::endl;
  DT NupT = arma::sum(Nfi.at(0));
  LogOut << "Total N_up = " << NupT << std::endl;
  LogOut << " Down Spin - " << std::endl;
  LogOut << Nfi.at(1) << std::endl;
  DT NdnT = arma::sum(Nfi.at(1));
  LogOut << "Total N_down = " << NdnT << std::endl;
  LogOut << " N_i - " << std::endl;
  RealVectorType Niall = Nfi.at(0) + Nfi.at(1);
  LogOut << Niall << std::endl;
  RealMatrixType NdnNup = NiNj( Bases, Vec, Ham0, 1, 0 );
  LogOut << " Correlation NdnNup" << std::endl;
  LogOut << NdnNup << std::endl;
  RealMatrixType NupNdn = NiNj( Bases, Vec, Ham0, 0, 1 );
  LogOut << " Correlation NupNdn" << std::endl;
  LogOut << NupNdn << std::endl;
  RealMatrixType NupNup = NiNj( Bases, Vec, Ham0, 0, 0 );
  LogOut << " Correlation NupNup" << std::endl;
  LogOut << NupNup << std::endl;
  RealMatrixType NdnNdn = NiNj( Bases, Vec, Ham0, 1, 1 );
  LogOut << " Correlation NdnNdn" << std::endl;
  LogOut << NdnNdn << std::endl;
  LogOut << " N_i^2 - " << std::endl;
  RealMatrixType Ni2 = NupNup.diag() + NdnNdn.diag() + 2.0e0 * NupNdn.diag();
  LogOut << Ni2 << std::endl;
  DTM CMUp = SingleParticleDensityMatrix( 0, Bases, Vec, Ham0 );
  LogOut << CMUp << std::endl;
  DTM CMDn = SingleParticleDensityMatrix( 1, Bases, Vec, Ham0 );
  LogOut << CMDn << std::endl;
  std::vector<ComplexMatrixType> BOW = BOWCorrelation( Bases, Vec, Ham0 );
  file->SaveVector("Obs", "Nup", Nfi.at(0));
  file->SaveNumber("Obs", "NupT", NupT);
  file->SaveVector("Obs", "Ndn", Nfi.at(1));
  file->SaveNumber("Obs", "NdnT", NdnT);
  file->SaveMatrix("Obs", "NdnNup", NdnNup);
  file->SaveMatrix("Obs", "NupNdn", NupNdn);
  file->SaveMatrix("Obs", "NupNup", NupNup);
  file->SaveMatrix("Obs", "NdnNdn", NdnNdn);
  file->SaveMatrix("Obs", "CMUp", CMUp);
  file->SaveMatrix("Obs", "CMDn", CMDn);
  file->SaveMatrix("Obs", "HopUpUp", BOW.at(0));
  file->SaveMatrix("Obs", "HopUpDn", BOW.at(1));
  file->SaveMatrix("Obs", "HopDnUp", BOW.at(2));
  file->SaveMatrix("Obs", "HopDnDn", BOW.at(3));
  delete file;
  LogOut.close();
}

void Spectral(const std::string prefix){
  std::ofstream LogOut;
  LogOut.open(prefix + "Spectral.log", std::ios::app);
  int L;
  int OBC;
  int N1, N2;
  std::vector<RealType> Jeqm, Ueqm, Veqm, Weqm, Uch, Vch;
  int TSteps;
  RealType dt;
  int CoreHole, Species, Type;
  try{
    /* Load parameters from file */
    H5::Exception::dontPrint();
    H5::H5File::isHdf5(prefix + "conf.h5");
    LoadEqmParameters( prefix + "conf.h5", L, OBC, N1, N2, Jeqm, Ueqm, Veqm, Weqm);
    LoadXASParameters( prefix + "conf.h5", Uch, Vch, TSteps, dt, CoreHole, Species, Type);
  }catch(H5::FileIException){
    L = 12;
    OBC = 0;
    N1 = 6;
    N2 = 6;
    Jeqm = std::vector<RealType>(L, 1.0);// OBC
    Ueqm = std::vector<RealType>(L, 6.0);
    Veqm = std::vector<RealType>(L, 0.0);
    Weqm = std::vector<RealType>(L, 0.0);
    Uch = std::vector<RealType>(L, 0.0);
    Vch = std::vector<RealType>(L, 0.0);
    CoreHole = L / 2;
    // Vch.at(CoreHole) = -5.0;
    Species = 0;
    Type = 1;
  }
  //* Eqm or Pump
  LogOut << "Build Eqm/Pump Basis - " << std::flush;
  Basis F1EQM(L, N1, true);
  F1EQM.Fermion();
  Basis F2EQM(L, N2, true);
  F2EQM.Fermion();
  std::vector<Basis> EqmBases;
  EqmBases.push_back(F1EQM);
  EqmBases.push_back(F2EQM);
  LogOut << "DONE!" << std::endl;
  LogOut << "Build Eqm Lattice (shared with core Hole) - " << std::flush;
  std::vector<ComplexType> J(Jeqm.begin(), Jeqm.end());
  std::vector< Node<ComplexType>* > Lattice = NN_1D_Chain(L, J, OBC);
  LogOut << "DONE!" << std::endl;
  LogOut << "Build Eqm Hamiltonian - " << std::flush;
  FHM<ComplexType> EqmHam( EqmBases );
  //* Potential
  std::vector<ComplexType> Vw(Veqm.begin(), Veqm.end());
  std::vector< std::vector<ComplexType> > EqmVloc = vec(Vw, Vw);
  //* Interaction
  std::vector<ComplexType> EqmUloc(Ueqm.begin(), Ueqm.end());
  // EqmHam.FermiHubbardModel(EqmBases, Lattice, EqmVloc, EqmUloc);
  //* Exntended Hubbard
  std::vector<ComplexType> EqmWloc(Weqm.begin(), Weqm.end());
  EqmHam.ExtendedFermiHubbardModel(EqmBases, Lattice, EqmVloc, EqmUloc, EqmWloc);
  LogOut << "Hermitian = " << EqmHam.CheckHermitian() << ", Hilbert space = " << EqmHam.GetTotalHilbertSpace() << ", DONE!" << std::endl;
  //* Create excitation
  LogOut << "Build Excited Basis - " << std::flush;
  int N1CH, N2CH;
  if ( Species == 0 ){
    if ( Type == 1 ) N1CH = N1 + 1;
    else if ( Type == -1 ) N1CH = N1 - 1;
    N2CH = N2;
  }else if ( Species == 1 ){
    N1CH = N1;
    if ( Type == 1 ) N2CH = N2 + 1;
    else if ( Type == -1 ) N2CH = N2 - 1;
  }
  Basis F1CH(L, N1CH, true);
  F1CH.Fermion();
  Basis F2CH(L, N2CH, true);
  F2CH.Fermion();
  std::vector<Basis> CoreHoleBases;
  CoreHoleBases.push_back(F1CH);
  CoreHoleBases.push_back(F2CH);
  LogOut << "DONE!" << std::endl;
  LogOut << "Build Excited Hamiltonian - " << std::flush;
  FHM<ComplexType> CoreHoleHam( CoreHoleBases );
  //* Potential
  std::vector<ComplexType> CoreHoleVw(Vch.begin(), Vch.end());
  std::vector< std::vector<ComplexType> > CoreHoleVloc = vec(CoreHoleVw, CoreHoleVw);
  //* Interaction
  std::vector<ComplexType> CoreHoleUloc(Uch.begin(), Uch.end());
  // CoreHoleHam.FermiHubbardModel(CoreHoleBases, Lattice, CoreHoleVloc, CoreHoleUloc);
  CoreHoleHam.ExtendedFermiHubbardModel(CoreHoleBases, Lattice, CoreHoleVloc, CoreHoleUloc, EqmWloc);
  LogOut << "Hermitian = " << CoreHoleHam.CheckHermitian() << ", Hilbert space = " << CoreHoleHam.GetTotalHilbertSpace() << ", DONE!" << std::endl;

  //* Load Wavefunction
  RealVectorType ValInput;
  ComplexVectorType VecInput;
  try{
    H5::Exception::dontPrint();
    H5::H5File::isHdf5(prefix + "PumpWF.h5");
    LogOut << "Load Pump Wavefunction - " << std::flush;
    HDF5IO *WFFile = new HDF5IO(prefix + "PumpWF.h5");
    WFFile->LoadVector("WF", "Vec", VecInput);
    delete WFFile;
  }catch(H5::FileIException){
    try{
      H5::Exception::dontPrint();
      H5::H5File::isHdf5(prefix + "FHMChainData.h5");
      LogOut << "Load Eqm Wavefunction - " << std::flush;
      HDF5IO *WFFile = new HDF5IO(prefix + "FHMChainData.h5");
      WFFile->LoadVector("GS", "EVec", VecInput);
      WFFile->LoadVector("GS", "EVal", ValInput);
      delete WFFile;
    }catch(H5::FileIException){
      RUNTIME_ERROR(" No input file with initial wavefunction, Pump.h5 or FHMChainData.h5!!");
    }
  }
  LogOut << VecInput.n_rows << " DONE!" << std::endl;
  LogOut << "Operate on Wave function on site @ " << CoreHole << ", on species - " << Species << ", type - " << Type << std::flush;
  ComplexVectorType VecInit = Operate( VecInput, CoreHole, Species, Type, EqmBases, CoreHoleBases, EqmHam, CoreHoleHam );
  VecInit = arma::normalise(VecInit);
  LogOut << ", DONE!" << std::endl;

  LogOut << "Get spectrum ... " << std::flush;
  const size_t MaxNumPeak = 40;
  std::vector<double> PeakLocations, PeakWeights;
  DTM Vecs(VecInit.n_rows, MaxNumPeak);
  RealVectorType Vals;
  CoreHoleHam.SpectralH(Vals, Vecs, VecInit, MaxNumPeak);
  for ( int iL = 0; iL < L; iL++ ){
    ComplexVectorType KetState = Operate( VecInput, iL, Species, Type, EqmBases, CoreHoleBases, EqmHam, CoreHoleHam );
    SpectralPeaks( ValInput(0), KetState, Vals, Vecs, PeakLocations, PeakWeights, MaxNumPeak);
    LogOut << " Total Weights = " << std::accumulate(PeakWeights.begin(), PeakWeights.end(), 0.0e0) << std::flush;
    LogOut << " DONE!" << std::endl;
    HDF5IO* file1 = new HDF5IO("Spectral.h5");
    std::string gname = "C";
    gname.append( std::to_string((unsigned long long)iL ));
    gname.append( "-C" );
    gname.append( std::to_string((unsigned long long)CoreHole ));
    gname.append( "-" );
    gname.append( std::to_string((unsigned long long)Type) );
    file1->SaveStdVector(gname, "ei", PeakLocations);
    file1->SaveStdVector(gname, "wi", PeakWeights);
    delete file1;
  }
}

void PumpDynamics(const std::string prefix, const int MeasureEvery = 10, const int SaveWFEvery = 1000000 ){
  std::ofstream LogOut;
  LogOut.open(prefix + "Pump.log", std::ios::app);
  int L;
  int OBC;
  int N1, N2;
  std::vector<RealType> Jeqm, Ueqm, Veqm, Weqm, At;
  int TSteps;
  RealType dt;
  try{
    /* Load parameters from file */
    H5::Exception::dontPrint();
    H5::H5File::isHdf5("conf.h5");
    LoadEqmParameters( "conf.h5", L, OBC, N1, N2, Jeqm, Ueqm, Veqm, Weqm);
    LoadPumpParameters( "conf.h5", At, TSteps, dt);
  }catch(H5::FileIException){
    L = 4;
    OBC = 0;
    N1 = 2;
    N2 = 2;
    Jeqm = std::vector<RealType>(L, 1.0);// PBC
    Ueqm = std::vector<RealType>(L, 6.0);
    Veqm = std::vector<RealType>(L, 0.0);
    Weqm = std::vector<RealType>(L-1, 3.0);
    TSteps = 10;
    dt = 0.005;
    At = std::vector<RealType>(TSteps, 0.25 * PI);
  }
  LogOut << "Build Eqm Lattice - " << std::flush;
  std::vector<ComplexType> J(Jeqm.begin(), Jeqm.end());
  std::vector< Node<ComplexType>* > EqmLattice = NN_1D_Chain(L, J, OBC);
  LogOut << "DONE!" << std::endl;
  LogOut << "Build Eqm Basis (shared with pump) - " << std::flush;
  Basis F1(L, N1, true);
  F1.Fermion();
  Basis F2(L, N2, true);
  F2.Fermion();
  std::vector<Basis> Bases;
  Bases.push_back(F1);
  Bases.push_back(F2);
  LogOut << "DONE!" << std::endl;
  LogOut << "Build Eqm Hamiltonian - " << std::flush;
  FHM<ComplexType> Ham0( Bases );
  // Potential
  std::vector<ComplexType> Vtmp(Veqm.begin(), Veqm.end());
  std::vector< std::vector<ComplexType> > Vloc = vec(Vtmp, Vtmp);
  // Interaction
  std::vector<ComplexType> Uloc(Ueqm.begin(), Ueqm.end());
  // Ham0.FermiHubbardModel(Bases, EqmLattice, Vloc, Uloc);
  //* Exntended Hubbard
  std::vector<ComplexType> Wloc(Weqm.begin(), Weqm.end());
  Ham0.ExtendedFermiHubbardModel(Bases, EqmLattice, Vloc, Uloc, Wloc);
  LogOut << "Hermitian = " << Ham0.CheckHermitian() << ", Hilbert space = " << Ham0.GetTotalHilbertSpace() << ", DONE!" << std::endl;
  // Load Wavefunction
  LogOut << "Load Eqm Wavefunction - " << std::flush;
  HDF5IO *EqmFile = new HDF5IO("FHMChainData.h5");
  ComplexVectorType VecInit;
  EqmFile->LoadVector("GS", "EVec", VecInit);
  delete EqmFile;
  LogOut << VecInit.n_rows << " DONE!" << std::endl;

  /* Pumping */
  ComplexType Prefactor = ComplexType(0.0, -1.0e0*dt);/* NOTE: hbar = 1 */
  ComplexVectorType VecPump = VecInit;
  HDF5IO* file2 = new HDF5IO("Pump.h5");
  std::string gname = "Obs-0/";
  ComplexType Lecho = arma::cdot(VecInit, VecPump);
  std::vector<RealVectorType> Nfi = Ni( Bases, VecPump, Ham0 );
  file2->SaveNumber(gname, "Lecho", Lecho);
  file2->SaveVector(gname, "Nup", Nfi.at(0));
  file2->SaveVector(gname, "Ndn", Nfi.at(1));
  delete file2;
  LogOut << "Begin pumping......" << std::endl;
  for (size_t cntP = 1; cntP <= TSteps; cntP++) {
    if ( OBC ){
      J = std::vector<ComplexType>(L - 1, exp( ComplexType(0.0, 1.0) * At.at(cntP-1)) );
    } else{
      J = std::vector<ComplexType>(L, exp( ComplexType(0.0, 1.0) * At.at(cntP-1)) );
    }
    std::vector< Node<ComplexType>* > PLattice = NN_1D_Chain(L, J, OBC);
    // Ham0.FermiHubbardModel(Bases, PLattice, Vloc, Uloc);
    Ham0.ExtendedFermiHubbardModel(Bases, PLattice, Vloc, Uloc, Wloc);
    if ( !(Ham0.CheckHermitian()) ) LogOut << "non-Hermitian Hamiltonian at PStep = " << cntP << std::endl;
    // Evolve the state
    // std::cout << cntP << " " << Ham0.CheckHermitian() << " " << arma::norm(VecPump) << std::flush;
    Ham0.expH(Prefactor, VecPump);
    // std::cout << " " << arma::norm(VecPump) << std::endl;
    VecPump = arma::normalise(VecPump);
    if ( cntP % MeasureEvery == 0 ){
      file2 = new HDF5IO("Pump.h5");
      std::string gname = "Obs-";
      gname.append( std::to_string((unsigned long long)cntP ));
      gname.append("/");
      ComplexType Lecho = arma::cdot(VecInit, VecPump);
      Nfi = Ni( Bases, VecPump, Ham0 );
      file2->SaveNumber(gname, "Lecho", Lecho);
      file2->SaveVector(gname, "Nup", Nfi.at(0));
      file2->SaveVector(gname, "Ndn", Nfi.at(1));
      delete file2;
    }
    if ( cntP % SaveWFEvery == 0 ){
      HDF5IO* file3 = new HDF5IO("PumpWF.h5");
      gname = "WF-";
      gname.append( std::to_string((unsigned long long)cntP ));
      gname.append("/");
      file3->SaveVector(gname, "Vec", VecPump);
      delete file3;
    }
  }
  // Ham0.FermiHubbardModel(Bases, EqmLattice, Vloc, Uloc);
  Ham0.ExtendedFermiHubbardModel(Bases, EqmLattice, Vloc, Uloc, Wloc);
  ComplexSparseMatrixType H0 = Ham0.GetTotalHamiltonian();
  ComplexType Energy = arma::cdot(VecPump, H0 * VecPump);
  HDF5IO* file3 = new HDF5IO("PumpWF.h5");
  gname = "WF";
  file3->SaveVector(gname, "Vec", VecPump);
  file3->SaveNumber(gname, "Energy", Energy);
  delete file3;
  LogOut << "Finished pumping!!" << std::endl;
  LogOut.close();
}

void StateDynamics(const std::string prefix, const int MeasureEvery = 50, const int SaveWFEvery = 1000000 ){
  std::ofstream LogOut;
  LogOut.open(prefix + "QuenchState.log", std::ios::app);
  int L;
  int OBC;
  int N1, N2;
  std::vector<RealType> Jeqm, Ueqm, Veqm, Weqm;
  int TSteps;
  RealType dt;
  try{
    /* Load parameters from file */
    H5::Exception::dontPrint();
    H5::H5File::isHdf5( prefix + "conf.h5" );
    LoadEqmParameters( prefix + "conf.h5", L, OBC, N1, N2, Jeqm, Ueqm, Veqm, Weqm);
    // std::vector<RealType> At;
    // LoadPumpParameters( "conf.h5", At, TSteps, dt);
    std::vector<RealType> Uch, Vch;
    int CoreHole, Species, Type;
    LoadXASParameters( prefix + "conf.h5", Uch, Vch, TSteps, dt, CoreHole, Species, Type);
  }catch(H5::FileIException){
    L = 14;
    OBC = 1;
    N1 = 7;
    N2 = 7;
    Jeqm = std::vector<RealType>(L-1, 1.0);// OBC
    Ueqm = std::vector<RealType>(L, 1.0);
    Veqm = std::vector<RealType>(L, 0.0);
    TSteps = 10;
    dt = 0.005;
  }
  LogOut << "Build Lattice - " << std::flush;
  std::vector<ComplexType> J(Jeqm.begin(), Jeqm.end());
  std::vector< Node<ComplexType>* > EqmLattice = NN_1D_Chain(L, J, OBC);
  LogOut << "DONE!" << std::endl;
  LogOut << "Build Basis - " << std::flush;
  Basis F1(L, N1, true);
  F1.Fermion();
  Basis F2(L, N2, true);
  F2.Fermion();
  std::vector<Basis> Bases;
  Bases.push_back(F1);
  Bases.push_back(F2);
  LogOut << "DONE!" << std::endl;
  LogOut << "Build Hamiltonian - " << std::flush;
  FHM<ComplexType> Ham0( Bases );
  // Potential
  std::vector<ComplexType> Vtmp(Veqm.begin(), Veqm.end());
  std::vector< std::vector<ComplexType> > Vloc = vec(Vtmp, Vtmp);
  // Interaction
  std::vector<ComplexType> Uloc(Ueqm.begin(), Ueqm.end());
  // Ham0.FermiHubbardModel(Bases, EqmLattice, Vloc, Uloc);
  //* Exntended Hubbard
  std::vector<ComplexType> Wloc(Weqm.begin(), Weqm.end());
  Ham0.ExtendedFermiHubbardModel(Bases, EqmLattice, Vloc, Uloc, Wloc);
  LogOut << "Hermitian = " << Ham0.CheckHermitian() << ", Hilbert space = " << Ham0.GetTotalHilbertSpace() << ", DONE!" << std::endl;
  // Load Wavefunction
  ComplexVectorType VecInit;
  try{
    H5::Exception::dontPrint();
    H5::H5File::isHdf5(prefix + "PumpWF.h5");
    LogOut << "Load Pump Wavefunction - " << std::flush;
    HDF5IO *WFFile = new HDF5IO(prefix + "PumpWF.h5");
    WFFile->LoadVector("WF", "Vec", VecInit);
    delete WFFile;
  }catch(H5::FileIException){
    try{
      H5::Exception::dontPrint();
      H5::H5File::isHdf5(prefix + "FHMChainData.h5");
      LogOut << "Load Eqm Wavefunction - " << std::flush;
      HDF5IO *WFFile = new HDF5IO(prefix + "FHMChainData.h5");
      WFFile->LoadVector("GS", "EVec", VecInit);
      delete WFFile;
    }catch(H5::FileIException){
      RUNTIME_ERROR(" No input file with initial wavefunction, Pump.h5 or FHMChainData.h5!!");
    }
  }
  LogOut << VecInit.n_rows << " DONE!" << std::endl;

  /* Pumping */
  ComplexType Prefactor = ComplexType(0.0, -1.0e0*dt);/* NOTE: hbar = 1 */
  ComplexVectorType VecDyn = VecInit;
  HDF5IO* file2 = new HDF5IO("QuenchState.h5");
  std::string gname = "Obs-0/";
  ComplexType Lecho = arma::cdot(VecInit, VecDyn);
  std::vector<RealVectorType> Nfi = Ni( Bases, VecDyn, Ham0 );
  file2->SaveNumber(gname, "Lecho", Lecho);
  file2->SaveVector(gname, "Nup", Nfi.at(0));
  file2->SaveVector(gname, "Ndn", Nfi.at(1));
  delete file2;
  LogOut << "Begin dynamics......" << std::endl;
  for (size_t cntP = 1; cntP <= TSteps; cntP++) {
    // Evolve the state
    Ham0.expH(Prefactor, VecDyn);
    if ( cntP % MeasureEvery == 0 ){
      file2 = new HDF5IO("QuenchState.h5");
      std::string gname = "Obs-";
      gname.append( std::to_string((unsigned long long)cntP ));
      gname.append("/");
      ComplexType Lecho = arma::cdot(VecInit, VecDyn);
      Nfi = Ni( Bases, VecDyn, Ham0 );
      file2->SaveNumber(gname, "Lecho", Lecho);
      file2->SaveVector(gname, "Nup", Nfi.at(0));
      file2->SaveVector(gname, "Ndn", Nfi.at(1));
      delete file2;
    }
    if ( cntP % SaveWFEvery == 0 ){
      HDF5IO* file3 = new HDF5IO("QuenchStateWF.h5");
      gname = "WF-";
      gname.append( std::to_string((unsigned long long)cntP ));
      gname.append("/");
      file3->SaveVector(gname, "Vec", VecDyn);
      delete file3;
    }
  }
  // Ham0.FermiHubbardModel(Bases, EqmLattice, Vloc, Uloc);
  Ham0.ExtendedFermiHubbardModel(Bases, EqmLattice, Vloc, Uloc, Wloc);
  ComplexSparseMatrixType H0 = Ham0.GetTotalHamiltonian();
  ComplexType Energy = arma::cdot(VecDyn, H0 * VecDyn);
  HDF5IO* file3 = new HDF5IO("QuenchStateWF.h5");
  gname = "WF";
  file3->SaveVector(gname, "Vec", VecDyn);
  file3->SaveNumber(gname, "Energy", Energy);
  delete file3;
  LogOut << "Finished dynamics!!" << std::endl;
  LogOut.close();
}

void SpectralDynamics(const std::string prefix, const int MeasureEvery = 2, const int SaveWFEvery = 1000000 ){
  std::ofstream LogOut;
  LogOut.open(prefix + "SpectralD.log", std::ios::app);
  int L;
  int OBC;
  int N1, N2;
  std::vector<RealType> Jeqm, Ueqm, Veqm, Weqm, Uch, Vch;
  int TSteps;
  RealType dt;
  int CoreHole, Species, Type;
  try{
    /* Load parameters from file */
    H5::Exception::dontPrint();
    H5::H5File::isHdf5(prefix + "conf.h5");
    LoadEqmParameters( prefix + "conf.h5", L, OBC, N1, N2, Jeqm, Ueqm, Veqm, Weqm);
    LoadXASParameters( prefix + "conf.h5", Uch, Vch, TSteps, dt, CoreHole, Species, Type);
  }catch(H5::FileIException){
    L = 12;
    OBC = 0;
    N1 = 6;
    N2 = 6;
    Jeqm = std::vector<RealType>(L, 1.0);// OBC
    Ueqm = std::vector<RealType>(L, 6.0);
    Veqm = std::vector<RealType>(L, 0.0);
    Weqm = std::vector<RealType>(L, 3.0);
    Uch = std::vector<RealType>(L, 0.0);
    Vch = std::vector<RealType>(L, 0.0);
    CoreHole = L/2;
    // Vch.at(CoreHole) = -5.0;
    TSteps = 3000;
    dt = 0.005;
    Species = 0;
    Type = 1;
  }
  /* Eqm or Pump */
  LogOut << "Build Eqm/Pump Basis - " << std::flush;
  Basis F1EQM(L, N1, true);
  F1EQM.Fermion();
  Basis F2EQM(L, N2, true);
  F2EQM.Fermion();
  std::vector<Basis> EqmBases;
  EqmBases.push_back(F1EQM);
  EqmBases.push_back(F2EQM);
  LogOut << "DONE!" << std::endl;
  LogOut << "Build Eqm Lattice (shared with core Hole) - " << std::flush;
  std::vector<ComplexType> J(Jeqm.begin(), Jeqm.end());
  std::vector< Node<ComplexType>* > Lattice = NN_1D_Chain(L, J, OBC);
  LogOut << "DONE!" << std::endl;
  LogOut << "Build Eqm Hamiltonian - " << std::flush;
  FHM<ComplexType> EqmHam( EqmBases );
  // Potential
  std::vector<ComplexType> Vw(Veqm.begin(), Veqm.end());
  std::vector< std::vector<ComplexType> > EqmVloc = vec(Vw, Vw);
  // Interaction
  std::vector<ComplexType> EqmUloc(Ueqm.begin(), Ueqm.end());
  //* Hubbard
  // EqmHam.FermiHubbardModel(EqmBases, Lattice, EqmVloc, EqmUloc);
  //* Exntended Hubbard
  std::vector<ComplexType> EqmWloc(Weqm.begin(), Weqm.end());
  EqmHam.ExtendedFermiHubbardModel(EqmBases, Lattice, EqmVloc, EqmUloc, EqmWloc);
  LogOut << "Hermitian = " << EqmHam.CheckHermitian() << ", Hilbert space = " << EqmHam.GetTotalHilbertSpace() << ", DONE!" << std::endl;
  /* Core Hole */
  LogOut << "Build Core Hole Basis - " << std::flush;
  int N1CH, N2CH;
  if ( Species == 0 ){
    if ( Type == 1 ) N1CH = N1 + 1;
    else if ( Type == -1 ) N1CH = N1 - 1;
    N2CH = N2;
  }else if ( Species == 1 ){
    N1CH = N1;
    if ( Type == 1 ) N2CH = N2 + 1;
    else if ( Type == -1 ) N2CH = N2 - 1;
  }
  Basis F1CH(L, N1CH, true);
  F1CH.Fermion();
  Basis F2CH(L, N2CH, true);
  F2CH.Fermion();
  std::vector<Basis> CoreHoleBases;
  CoreHoleBases.push_back(F1CH);
  CoreHoleBases.push_back(F2CH);
  LogOut << "DONE!" << std::endl;
  LogOut << "Build Core Hole Hamiltonian - " << std::flush;
  FHM<ComplexType> CoreHoleHam( CoreHoleBases );
  // Potential
  std::vector<ComplexType> CoreHoleVw(Vch.begin(), Vch.end());
  std::vector< std::vector<ComplexType> > CoreHoleVloc = vec(CoreHoleVw, CoreHoleVw);
  // Interaction
  std::vector<ComplexType> CoreHoleUloc(Uch.begin(), Uch.end());
  // CoreHoleHam.FermiHubbardModel(CoreHoleBases, Lattice, CoreHoleVloc, CoreHoleUloc);
  CoreHoleHam.ExtendedFermiHubbardModel(CoreHoleBases, Lattice, CoreHoleVloc, CoreHoleUloc, EqmWloc);
  LogOut << "Hermitian = " << CoreHoleHam.CheckHermitian() << ", Hilbert space = " << CoreHoleHam.GetTotalHilbertSpace() << ", DONE!" << std::endl;
  // Load Wavefunction
  ComplexVectorType VecInput;
  try{
    H5::Exception::dontPrint();
    H5::H5File::isHdf5(prefix + "PumpWF.h5");
    LogOut << "Load Pump Wavefunction - " << std::flush;
    HDF5IO *WFFile = new HDF5IO(prefix + "PumpWF.h5");
    WFFile->LoadVector("WF", "Vec", VecInput);
    delete WFFile;
  }catch(H5::FileIException){
    try{
      H5::Exception::dontPrint();
      H5::H5File::isHdf5(prefix + "FHMChainData.h5");
      LogOut << "Load Eqm Wavefunction - " << std::flush;
      HDF5IO *WFFile = new HDF5IO(prefix + "FHMChainData.h5");
      WFFile->LoadVector("GS", "EVec", VecInput);
      delete WFFile;
    }catch(H5::FileIException){
      RUNTIME_ERROR(" No input file with initial wavefunction, Pump.h5 or FHMChainData.h5!!");
    }
  }
  LogOut << VecInput.n_rows << " DONE!" << std::endl;
  // Create Core Hole on wf
  LogOut << "Operate on Wave function with core hole @ " << CoreHole << ", on species - " << Species << ", type - " << Type << std::flush;
  ComplexVectorType VecInit = Operate( VecInput, CoreHole, Species, Type, EqmBases, CoreHoleBases, EqmHam, CoreHoleHam );
  ComplexVectorType VecXAS = arma::normalise(VecInit);
  HDF5IO* file2 = new HDF5IO("SpectralDWF.h5");
  std::string gname = "WF-0";
  file2->SaveVector(gname, "Vec", VecXAS);
  file2->SaveNumber(gname, "Norm", arma::norm(VecInit));
  delete file2;
  LogOut << ", DONE!" << std::endl;

  ComplexType Prefactor = ComplexType(0.0, -1.0e0*dt);/* NOTE: hbar = 1 */
  HDF5IO* file1 = new HDF5IO("SpectralD.h5");
  gname = "Obs-0/";
  ComplexType Lecho = arma::cdot(VecInit, VecXAS);
  std::vector<RealVectorType> Nfi = Ni( CoreHoleBases, VecXAS, CoreHoleHam );
  file1->SaveNumber(gname, "Lecho", Lecho);
  file1->SaveVector(gname, "Nup", Nfi.at(0));
  file1->SaveVector(gname, "Ndn", Nfi.at(1));
  delete file1;
  LogOut << "Begin dynamics to XAS......" << std::endl;
  for (size_t cntT = 1; cntT <= TSteps; cntT++) {
    // Evolve the state
    CoreHoleHam.expH(Prefactor, VecXAS);
    VecXAS = arma::normalise(VecXAS);/* NOTE: Does not change much. */
    if ( cntT % MeasureEvery == 0 ){
      file1 = new HDF5IO("SpectralD.h5");
      std::string gname = "Obs-";
      gname.append( std::to_string((unsigned long long)cntT ));
      gname.append("/");
      ComplexType Lecho = arma::cdot(VecInit, VecXAS);
      Nfi = Ni( CoreHoleBases, VecXAS, CoreHoleHam );
      file1->SaveNumber(gname, "Lecho", Lecho);
      file1->SaveVector(gname, "Nup", Nfi.at(0));
      file1->SaveVector(gname, "Ndn", Nfi.at(1));
      delete file1;
    }
    if ( cntT % SaveWFEvery == 0 ){
      file2 = new HDF5IO("SpectralDWF.h5");
      gname = "WF-";
      gname.append( std::to_string((unsigned long long)cntT ));
      gname.append("/");
      file2->SaveVector(gname, "Vec", VecXAS);
      delete file2;
    }
  }
  file2 = new HDF5IO("SpectralDWF.h5");
  gname = "WF";
  file2->SaveVector(gname, "Vec", VecXAS);
  delete file2;
  LogOut << "Finished Spectral from dynamics!!" << std::endl;
  LogOut.close();
}

void CalculateObs(const std::string prefix, const int Every){
  std::ofstream LogOut;
  LogOut.open(prefix + "Obs.log", std::ios::app);
  int L;
  int OBC;
  int N1, N2;
  std::vector<RealType> Jeqm, Ueqm, Veqm, Weqm, At;
  int TSteps;
  RealType dt;
  try{
    /* Load parameters from file */
    H5::Exception::dontPrint();
    H5::H5File::isHdf5(prefix + "conf.h5");
    LoadEqmParameters( prefix + "conf.h5", L, OBC, N1, N2, Jeqm, Ueqm, Veqm, Weqm);
    LoadPumpParameters( "conf.h5", At, TSteps, dt);
  }catch(H5::FileIException){
    L = 12;
    OBC = 0;
    N1 = 6;
    N2 = 6;
    Jeqm = std::vector<RealType>(L, 1.0);// OBC
    Ueqm = std::vector<RealType>(L, 6.0);
    Veqm = std::vector<RealType>(L, 0.0);
    Weqm = std::vector<RealType>(L, 3.0);
    TSteps = 3000;
    dt = 0.005;
  }
  LogOut << "Build Eqm/Pump Basis - " << std::flush;
  Basis F1EQM(L, N1, true);
  F1EQM.Fermion();
  Basis F2EQM(L, N2, true);
  F2EQM.Fermion();
  std::vector<Basis> EqmBases;
  EqmBases.push_back(F1EQM);
  EqmBases.push_back(F2EQM);
  LogOut << "DONE!" << std::endl;
  LogOut << "Build Eqm Hamiltonian - " << std::flush;
  FHM<ComplexType> EqmHam( EqmBases );
  LogOut << " (only for determine index) DONE!" << std::endl;
  ComplexVectorType Vec;
  for (size_t cnt = Every; cnt <= TSteps; cnt += Every) {
    HDF5IO* file1 = new HDF5IO(prefix + "PumpWF.h5");
    std::string gname = "WF-";
    gname.append( std::to_string((unsigned long long)cnt ));
    gname.append("/");
    file1->LoadVector(gname, "Vec", Vec);
    std::vector<RealVectorType> Nfi = Ni( EqmBases, Vec, EqmHam );
    RealMatrixType NupNdn = NiNj( EqmBases, Vec, EqmHam, 0, 1 );
    RealMatrixType NdnNup = NiNj( EqmBases, Vec, EqmHam, 1, 0 );
    RealMatrixType NupNup = NiNj( EqmBases, Vec, EqmHam, 0, 0 );
    RealMatrixType NdnNdn = NiNj( EqmBases, Vec, EqmHam, 1, 1 );
    ComplexMatrixType CMUp = SingleParticleDensityMatrix( 0, EqmBases, Vec, EqmHam );
    ComplexMatrixType CMDn = SingleParticleDensityMatrix( 1, EqmBases, Vec, EqmHam );
    std::vector<ComplexMatrixType> BOW = BOWCorrelation( EqmBases, Vec, EqmHam );
    HDF5IO* file2 = new HDF5IO(prefix + "PumpObs.h5");
    gname = "obs-";
    gname.append( std::to_string((unsigned long long)cnt ));
    file2->SaveVector(gname, "Nup", Nfi.at(0));
    file2->SaveVector(gname, "Ndn", Nfi.at(1));
    file2->SaveMatrix(gname, "NupNdn", NupNdn);
    file2->SaveMatrix(gname, "NdnNup", NdnNup);
    file2->SaveMatrix(gname, "NupNup", NupNup);
    file2->SaveMatrix(gname, "NdnNdn", NdnNdn);
    file2->SaveMatrix(gname, "HopUpUp", BOW.at(0));
    file2->SaveMatrix(gname, "HopUpDn", BOW.at(1));
    file2->SaveMatrix(gname, "HopDnUp", BOW.at(2));
    file2->SaveMatrix(gname, "HopDnDn", BOW.at(3));
    file2->SaveMatrix(gname, "CMUp", CMUp);
    file2->SaveMatrix(gname, "CMDn", CMDn);
    delete file2;
    delete file1;
    LogOut << cnt << ", " << std::flush;
  }
  LogOut << "DONE." << std::endl;
  LogOut.close();
}

void CalculateObs2(const std::string prefix, const int Every){
  std::ofstream LogOut;
  LogOut.open(prefix + "Obs2.log", std::ios::app);
  int L;
  int OBC;
  int N1, N2;
  std::vector<RealType> Jeqm, Ueqm, Veqm, Weqm, At;
  int TSteps;
  RealType dt;
  try{
    /* Load parameters from file */
    H5::Exception::dontPrint();
    H5::H5File::isHdf5(prefix + "conf.h5");
    LoadEqmParameters( prefix + "conf.h5", L, OBC, N1, N2, Jeqm, Ueqm, Veqm, Weqm);
    std::vector<RealType> Uch, Vch;
    int CoreHole, Species, Type;
    LoadXASParameters( prefix + "conf.h5", Uch, Vch, TSteps, dt, CoreHole, Species, Type);
  }catch(H5::FileIException){
    L = 12;
    OBC = 0;
    N1 = 6;
    N2 = 6;
    Jeqm = std::vector<RealType>(L, 1.0);// OBC
    Ueqm = std::vector<RealType>(L, 6.0);
    Veqm = std::vector<RealType>(L, 0.0);
    Weqm = std::vector<RealType>(L, 3.0);
    TSteps = 3000;
    dt = 0.005;
  }
  LogOut << "Build Eqm/Pump Basis - " << std::flush;
  Basis F1EQM(L, N1, true);
  F1EQM.Fermion();
  Basis F2EQM(L, N2, true);
  F2EQM.Fermion();
  std::vector<Basis> EqmBases;
  EqmBases.push_back(F1EQM);
  EqmBases.push_back(F2EQM);
  LogOut << "DONE!" << std::endl;
  LogOut << "Build Eqm Hamiltonian - " << std::flush;
  FHM<ComplexType> EqmHam( EqmBases );
  LogOut << " (only for determine index) DONE!" << std::endl;
  ComplexVectorType Vec;
  for (size_t cnt = Every; cnt <= TSteps; cnt += Every) {
    HDF5IO* file1 = new HDF5IO(prefix + "QuenchStateWF.h5");
    std::string gname = "WF-";
    gname.append( std::to_string((unsigned long long)cnt ));
    gname.append("/");
    file1->LoadVector(gname, "Vec", Vec);
    std::vector<RealVectorType> Nfi = Ni( EqmBases, Vec, EqmHam );
    RealMatrixType NupNdn = NiNj( EqmBases, Vec, EqmHam, 0, 1 );
    RealMatrixType NdnNup = NiNj( EqmBases, Vec, EqmHam, 1, 0 );
    RealMatrixType NupNup = NiNj( EqmBases, Vec, EqmHam, 0, 0 );
    RealMatrixType NdnNdn = NiNj( EqmBases, Vec, EqmHam, 1, 1 );
    ComplexMatrixType CMUp = SingleParticleDensityMatrix( 0, EqmBases, Vec, EqmHam );
    ComplexMatrixType CMDn = SingleParticleDensityMatrix( 1, EqmBases, Vec, EqmHam );
    std::vector<ComplexMatrixType> BOW = BOWCorrelation( EqmBases, Vec, EqmHam );
    HDF5IO* file2 = new HDF5IO(prefix + "QuenchStateObs.h5");
    gname = "obs-";
    gname.append( std::to_string((unsigned long long)cnt ));
    file2->SaveVector(gname, "Nup", Nfi.at(0));
    file2->SaveVector(gname, "Ndn", Nfi.at(1));
    file2->SaveMatrix(gname, "NupNdn", NupNdn);
    file2->SaveMatrix(gname, "NdnNup", NdnNup);
    file2->SaveMatrix(gname, "NupNup", NupNup);
    file2->SaveMatrix(gname, "NdnNdn", NdnNdn);
    file2->SaveMatrix(gname, "HopUpUp", BOW.at(0));
    file2->SaveMatrix(gname, "HopUpDn", BOW.at(1));
    file2->SaveMatrix(gname, "HopDnUp", BOW.at(2));
    file2->SaveMatrix(gname, "HopDnDn", BOW.at(3));
    file2->SaveMatrix(gname, "CMUp", CMUp);
    file2->SaveMatrix(gname, "CMDn", CMDn);
    delete file2;
    delete file1;
    LogOut << cnt << ", " << std::flush;
  }
  LogOut << "DONE." << std::endl;
  LogOut.close();
}

void CalculateAt(const std::string prefix, const int Every, const int From = 0){
  std::ofstream LogOut;
  LogOut.open(prefix + "At.log", std::ios::app);
  int L;
  int OBC;
  int N1, N2;
  std::vector<RealType> Jeqm, Ueqm, Veqm, Weqm, Uch, Vch;
  int TSteps;
  RealType dt;
  int CoreHole, Species, Type;
  try{
    /* Load parameters from file */
    H5::Exception::dontPrint();
    H5::H5File::isHdf5(prefix + "conf.h5");
    LoadEqmParameters( prefix + "conf.h5", L, OBC, N1, N2, Jeqm, Ueqm, Veqm, Weqm);
    LoadXASParameters( prefix + "conf.h5", Uch, Vch, TSteps, dt, CoreHole, Species, Type);
  }catch(H5::FileIException){
    L = 12;
    OBC = 0;
    N1 = 6;
    N2 = 6;
    Jeqm = std::vector<RealType>(L, 1.0);// OBC
    Ueqm = std::vector<RealType>(L, 6.0);
    Veqm = std::vector<RealType>(L, 0.0);
    Weqm = std::vector<RealType>(L, 3.0);
    Uch = std::vector<RealType>(L, 0.0);
    Vch = std::vector<RealType>(L, 0.0);
    CoreHole = L/2;
    // Vch.at(CoreHole) = -5.0;
    TSteps = 3000;
    dt = 0.005;
    Species = 0;
    Type = 1;
  }
  LogOut << "Build Eqm/Pump Basis - " << std::flush;
  Basis F1EQM(L, N1, true);
  F1EQM.Fermion();
  Basis F2EQM(L, N2, true);
  F2EQM.Fermion();
  std::vector<Basis> EqmBases;
  EqmBases.push_back(F1EQM);
  EqmBases.push_back(F2EQM);
  LogOut << "DONE!" << std::endl;
  // LogOut << "Build Eqm Lattice (shared with core Hole) - " << std::flush;
  // std::vector<ComplexType> J(Jeqm.begin(), Jeqm.end());
  // std::vector< Node<ComplexType>* > Lattice = NN_1D_Chain(L, J, OBC);
  // LogOut << "DONE!" << std::endl;
  LogOut << "Build Eqm Hamiltonian - " << std::flush;
  FHM<ComplexType> EqmHam( EqmBases );
  LogOut << " (only for determine index) DONE!" << std::endl;
  LogOut << "Build Core Hole Basis - " << std::flush;
  int N1CH, N2CH;
  if ( Species == 0 ){
    if ( Type == 1 ) N1CH = N1 + 1;
    else if ( Type == -1 ) N1CH = N1 - 1;
    N2CH = N2;
  }else if ( Species == 1 ){
    N1CH = N1;
    if ( Type == 1 ) N2CH = N2 + 1;
    else if ( Type == -1 ) N2CH = N2 - 1;
  }
  Basis F1CH(L, N1CH, true);
  F1CH.Fermion();
  Basis F2CH(L, N2CH, true);
  F2CH.Fermion();
  std::vector<Basis> CoreHoleBases;
  CoreHoleBases.push_back(F1CH);
  CoreHoleBases.push_back(F2CH);
  LogOut << "DONE!" << std::endl;
  LogOut << "Build Core Hole Hamiltonian - " << std::flush;
  FHM<ComplexType> CoreHoleHam( CoreHoleBases );
  LogOut << " (only for determine index) DONE!" << std::endl;

  for ( int cntL = From; cntL < L; cntL++ ){
    ComplexVectorType BraState, KetState, OpKetState;
    std::vector<ComplexType> At;
    At.clear();
    for (size_t cnt = 0; cnt <= TSteps; cnt += Every) {
      HDF5IO* file1 = new HDF5IO("SpectralDWF.h5");
      HDF5IO* file2 = new HDF5IO("QuenchStateWF.h5");
      std::string gname = "WF-";
      gname.append( std::to_string((unsigned long long)cnt ));
      gname.append("/");
      file1->LoadVector(gname, "Vec", BraState);
      if ( cnt == 0 ){
        OpKetState = BraState;
      }else{
        file2->LoadVector(gname, "Vec", KetState);
        OpKetState = Operate( KetState, CoreHole, Species, Type, EqmBases, CoreHoleBases, EqmHam, CoreHoleHam );
      }
      ComplexType val = arma::cdot(OpKetState, BraState);
      At.push_back(val);
      delete file2;
      delete file1;
    }
    HDF5IO* file3 = new HDF5IO("At.h5");
    std::string name = "CH";
    name.append( std::to_string((unsigned long long)cntL ));
    name.append("/");
    file3->SaveStdVector(name, "At", At);
    delete file3;
  }
  LogOut << "Finished calculating At!!" << std::endl;
  LogOut.close();
}

/* main program */
int main(int argc, char *argv[]){
  if ( argc < 2 ) RUNTIME_ERROR(" Need at least one argument to run program. Use 0 to run Equilibrium.");
#ifdef MKL
  mkl_set_num_threads(NumCores);
#endif
  if ( std::atoi(argv[1]) == 0 ){
    Equilibrium("");
  }else if ( std::atoi(argv[1]) == 1 ){
    PumpDynamics("", 20, 20);
  }else if ( std::atoi(argv[1]) == 2 ){
    SpectralDynamics("", 20, 20);
  }else if ( std::atoi(argv[1]) == 3 ){
    StateDynamics("", 20, 20);
  }else if ( std::atoi(argv[1]) == 4 ){
    Spectral("");
  }else if ( std::atoi(argv[1]) == 5 ) {
    int From = std::atoi(argv[2]);
    CalculateAt("", 20, From);
  }else if ( std::atoi(argv[1]) == 6 ) {
    CalculateObs("", 20);
  }else if ( std::atoi(argv[1]) == 7 ) {
    CalculateObs2("", 20);
  }
  return 0;
}
