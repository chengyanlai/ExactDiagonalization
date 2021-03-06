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

void LoadParameters( const std::string filename, int &L,
  int &OBC, int &N1, int &N2, std::vector<RealType> &Uinit,
  std::vector<RealType> &Uloc, std::vector<RealType> &Vloc,
  int &CHloc, int &dynamics, int &Tsteps, RealType &RealType){
    HDF5IO file(filename);
    L = file.loadInt("Parameters", "L");
    OBC = file.loadInt("Parameters", "OBC");
    N1 = file.loadInt("Parameters", "N1");
    N2 = file.loadInt("Parameters", "N2");
    CHloc = file.loadInt("Parameters", "CHloc");
    file.loadStdVector("Parameters", "Uinit", Uinit);
    file.loadStdVector("Parameters", "U", Uloc);
    file.loadStdVector("Parameters", "V", Vloc);
    dynamics = file.loadInt("Parameters", "dynamics");
    Tsteps = file.loadInt("Parameters", "Tsteps");
    RealType = file.loadReal("Parameters", "dt");
}

void peaks( const double EG, const ComplexVectorType AS, const RealVectorType &EigVal, const ComplexMatrixType &EigVec,
  std::vector<double> &PeakLocation, std::vector<double> &PeakWeight, const int MaxNumPeak){
  PeakLocation.clear();
  PeakWeight.clear();
  for ( size_t i = 0; i < EigVec.rows(); i++){
    PeakLocation.push_back(EigVal(i) - EG);
    ComplexVectorType An = EigVec.row(i);
    ComplexType val = An.dot(AS);
    PeakWeight.push_back( RealPart(val * Conjg(val)) );
    if ( i > MaxNumPeak ) break;
  }
  assert( PeakLocation.size() == PeakWeight.size() );
}

ComplexVectorType OperateCdagger( const std::vector<Basis> OldBases, const ComplexVectorType Vin,
  const std::vector<Basis> NewBases, const int CHloc, const int species,
  Hamiltonian<ComplexType> OldHam, Hamiltonian<ComplexType> NewHam  ){
  size_t NewHilbertSpace = NewBases.at(0).GetHilbertSpace() * NewBases.at(1).GetHilbertSpace();
  ComplexVectorType Vout = ComplexVectorType::Zero(NewHilbertSpace);
  std::vector<int> OldFup = OldBases.at(0).GetFStates();
  std::vector<int> OldFdn = OldBases.at(1).GetFStates();
  std::vector<size_t> OldFupTag = OldBases.at(0).GetFTags();
  std::vector<size_t> OldFdnTag = OldBases.at(1).GetFTags();
  std::vector<size_t> NewFupTag = NewBases.at(0).GetFTags();
  std::vector<size_t> NewFdnTag = NewBases.at(1).GetFTags();
  size_t NewFupIdx, NewFdnIdx, OldFupIdx, OldFdnIdx;
  // int Update = false;
  for ( auto OldFupState : OldFup ){
    if (species == 1){
      NewFupIdx = NewFupTag.at(OldFupState);// Find their indices
      OldFupIdx = OldFupTag.at(OldFupState);// Find their indices
    }else if ( (species == 0) && !(btest(OldFupState, CHloc)) ){
      int NewFupState = ibset(OldFupState, CHloc);
      NewFupIdx = NewFupTag.at(NewFupState);// Find their indices
      OldFupIdx = OldFupTag.at(OldFupState);// Find their indices
    }else{
      continue;
    }
    for ( auto OldFdnState : OldFdn ){
      if ( species == 0 ){
        NewFdnIdx = NewFdnTag.at(OldFdnState);// Find their indices
        OldFdnIdx = OldFdnTag.at(OldFdnState);// Find their indices
      }else if ( (species == 1) && !(btest(OldFdnState, CHloc)) ){
        int NewFdnState = ibset(OldFdnState, CHloc);
        NewFdnIdx = NewFdnTag.at(NewFdnState);// Find their indices
        OldFdnIdx = OldFdnTag.at(OldFdnState);// Find their indices
      }else{
        continue;
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
      Vout(nid) = Vin(oid);
    }
  }
  return Vout;
}

ComplexVectorType OperateC( const std::vector<Basis> OldBases, const ComplexVectorType Vin,
  const std::vector<Basis> NewBases, const int CHloc, const int species,
  Hamiltonian<ComplexType> OldHam, Hamiltonian<ComplexType> NewHam  ){
  size_t NewHilbertSpace = NewBases.at(0).GetHilbertSpace() * NewBases.at(1).GetHilbertSpace();
  ComplexVectorType Vout = ComplexVectorType::Zero(NewHilbertSpace);
  std::vector<int> OldFup = OldBases.at(0).GetFStates();
  std::vector<int> OldFdn = OldBases.at(1).GetFStates();
  std::vector<size_t> OldFupTag = OldBases.at(0).GetFTags();
  std::vector<size_t> OldFdnTag = OldBases.at(1).GetFTags();
  std::vector<size_t> NewFupTag = NewBases.at(0).GetFTags();
  std::vector<size_t> NewFdnTag = NewBases.at(1).GetFTags();
  size_t NewFupIdx, NewFdnIdx, OldFupIdx, OldFdnIdx;
  // int Update = false;
  for ( auto OldFupState : OldFup ){
    if (species == 1){
      NewFupIdx = NewFupTag.at(OldFupState);// Find their indices
      OldFupIdx = OldFupTag.at(OldFupState);// Find their indices
    }else if ( (species == 0) && (btest(OldFupState, CHloc)) ){
      int NewFupState = ibclr(OldFupState, CHloc);
      NewFupIdx = NewFupTag.at(NewFupState);// Find their indices
      OldFupIdx = OldFupTag.at(OldFupState);// Find their indices
    }else{
      continue;
    }
    for ( auto OldFdnState : OldFdn ){
      if ( species == 0 ){
        NewFdnIdx = NewFdnTag.at(OldFdnState);// Find their indices
        OldFdnIdx = OldFdnTag.at(OldFdnState);// Find their indices
      }else if ( (species == 1) && (btest(OldFdnState, CHloc)) ){
        int NewFdnState = ibclr(OldFdnState, CHloc);
        NewFdnIdx = NewFdnTag.at(NewFdnState);// Find their indices
        OldFdnIdx = OldFdnTag.at(OldFdnState);// Find their indices
      }else{
        continue;
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
      Vout(nid) = Vin(oid);
    }
  }
  return Vout;
}

std::vector< RealVectorType > Ni( const std::vector<Basis> &Bases, const ComplexVectorType &Vec, Hamiltonian<ComplexType> &ham ){
  std::vector< RealVectorType > out;
  RealVectorType tmp1 = RealVectorType::Zero(Bases.at(0).GetL());//(Bases.at(0).GetL(), 0.0e0);
  RealVectorType tmp2 = RealVectorType::Zero(Bases.at(1).GetL());//(Bases.at(1).GetL(), 0.0e0);
  std::vector< int > f1 = Bases.at(0).GetFStates();
  std::vector< int > f2 = Bases.at(1).GetFStates();
  size_t f1id = 0, f2id = 0;
  for ( int &nf2 : f2 ){
    std::vector<size_t> ids(2,f2id);
    f1id = 0;
    for ( int &nf1 : f1 ){
      ids.at(0) = f1id;
      size_t id = ham.DetermineTotalIndex(ids);
      for (size_t cnt = 0; cnt < Bases.at(0).GetL(); cnt++) {
        if ( btest(nf1, cnt) ) tmp1(cnt) += std::pow(std::abs(Vec(id)), 2);//Vec(id) * std::conj( Vec(id) );
      }
      for (size_t cnt = 0; cnt < Bases.at(1).GetL(); cnt++) {
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

int main(int argc, char const *argv[]) {
  Eigen::setNbThreads(NumCores);
#ifdef MKL
  mkl_set_num_threads(NumCores);
#endif
  INFO("Eigen3 uses " << Eigen::nbThreads() << " threads.");
  int L, CHloc;
  int OBC;
  int N1, N2;
  int dynamics, Tsteps;
  RealType dt;
  std::vector<RealType> Uinit, Vin, Uin;
  LoadParameters( "conf.h5", L, OBC, N1, N2, Uinit, Uin, Vin, CHloc, dynamics, Tsteps, dt);
  HDF5IO *file = new HDF5IO("XASEQM.h5");
  INFO("Build Lattice - ");
  std::vector<ComplexType> J;
  if ( OBC ){
    J = std::vector<ComplexType>(L - 1, 1.0);
  } else{
    J = std::vector<ComplexType>(L, 1.0);
  }
  for ( auto &val : J ){
    INFO_NONEWLINE(val << " ");
  }
  INFO("");
  const std::vector< Node<ComplexType>* > lattice = NN_1D_Chain(L, J, OBC);
  file->saveNumber("1DChain", "L", L);
  file->saveStdVector("1DChain", "J", J);
  for ( auto &lt : lattice ){
    if ( !(lt->VerifySite()) ) RUNTIME_ERROR("Wrong lattice setup!");
  }
  INFO("DONE!");
  INFO("Build Basis - ");
  Basis F1(L, N1, true);
  F1.Fermion();
  std::vector<int> st1 = F1.GetFStates();
  std::vector<size_t> tg1 = F1.GetFTags();
  Basis F2(L, N2, true);
  F2.Fermion();
  std::vector<int> st2 = F2.GetFStates();
  std::vector<size_t> tg2 = F2.GetFTags();
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
  Hamiltonian<ComplexType> ham( Bases );
  std::vector<ComplexType> Vtmp(L, ComplexType(0.0e0, 0.0e0) );
  std::vector<std::vector<ComplexType> > Vloc;
  Vloc.push_back(Vtmp);
  Vloc.push_back(Vtmp);
  std::vector<ComplexType> Utmp;
  for ( auto val : Uinit ){
    Utmp.push_back(ComplexType(val));
  }
  std::vector<std::vector<ComplexType> > Uloc;
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
  ComplexMatrixType Vecs;
  ham.eigh(Vals, Vecs, 1);
  ComplexVectorType Vec = Vecs.row(0);
  RealType EG = Vals[0];
  INFO("GS energy = " << Vals[0]);
  file->saveVector("GS", "EVec", Vec);
  file->saveVector("GS", "EVal", Vals);
  INFO("DONE!");
  std::vector< RealVectorType > Nfi = Ni( Bases, Vec, ham );
  RealVectorType Niall = Nfi.at(0) + Nfi.at(1);
  file->saveVector("Obs", "Nup", Nfi.at(0));
  file->saveVector("Obs", "Ndn", Nfi.at(1));
  std::cout << "Build core-hole Basis" << std::endl;
  /* Build New Basis */
  Basis nF1(L, N1+1, true);
  nF1.Fermion();
  Basis nF2(L, N2, true);
  nF2.Fermion();
  std::vector<Basis> nBases;
  nBases.push_back(nF1);
  nBases.push_back(nF2);
  /* Update the new Hamiltonian */
  std::cout << "Build core-hole Hamiltonian" << std::endl;
  Hamiltonian<ComplexType> nHam( nBases );
  Vtmp.clear();
  for ( auto val : Vin ){
    Vtmp.push_back(ComplexType(val));
  }
  Vloc.clear();
  Vloc.push_back(Vtmp);
  Vloc.push_back(Vtmp);
  Utmp.clear();
  for ( auto val : Uin ){
    Utmp.push_back(ComplexType(val));
  }
  Uloc.clear();
  Uloc.push_back(Utmp);
  Uloc.push_back(Utmp);
  nHam.BuildLocalHamiltonian(Vloc, Uloc, nBases);
  nHam.BuildHoppingHamiltonian(nBases, lattice);
  nHam.BuildTotalHamiltonian();
  /* Apply the operator */
  std::cout << "Create core-hole state" << std::endl;
  ComplexVectorType VecInit = OperateCdagger( Bases, Vec, nBases, CHloc, 0, ham, nHam);
  VecInit.normalize();
  ComplexVectorType VecT = VecInit;
  Nfi = Ni( nBases, VecT, nHam );
  Niall = Nfi.at(0) + Nfi.at(1);
  file->saveVector("Obs", "NewNup", Nfi.at(0));
  file->saveVector("Obs", "NewNdn", Nfi.at(1));
  delete file;
  if ( dynamics ){
    ComplexType Prefactor = ComplexType(0.0, -1.0e0*dt);/* NOTE: hbar = 1 */
    std::cout << "Begin dynamics......" << std::endl;
    for (size_t cntT = 1; cntT <= Tsteps; cntT++) {
      nHam.expH(Prefactor, VecT);
      if ( cntT % 2 == 0 ){
        HDF5IO file2("XASDYN.h5");
        std::string gname = "Obs-";
        gname.append(std::to_string((unsigned long long)cntT));
        gname.append("/");
        ComplexType Lecho = VecInit.dot(VecT);
        Nfi = Ni( nBases, VecT, nHam );
        file2.saveNumber(gname, "Lecho", Lecho);
        file2.saveVector(gname, "Nup", Nfi.at(0));
        file2.saveVector(gname, "Ndn", Nfi.at(1));
      }
    }
  }else{
    int MaxNumPeak = 20;
    std::vector<double> PeakLocation, PeakWeight;
    RealVectorType EigVals = RealVectorType::Zero(MaxNumPeak);
    ComplexMatrixType EigVecs = ComplexMatrixType::Zero(MaxNumPeak, nHam.GetTotalHilbertSpace());
    EigVecs.row(0) = VecInit;
    nHam.eigh(EigVals, EigVecs, MaxNumPeak, false);
    std::cout << "new " << EigVals[0] << " " << EigVals[1] << std::endl;
    peaks(EG, VecInit, EigVals, EigVecs, PeakLocation, PeakWeight, MaxNumPeak);
    /* Print to check. */
    for ( int i = 0; i < 4; i++ ) std::cout << PeakLocation.at(i) << " " << PeakWeight.at(i) << std::endl;
    HDF5IO file2("XAS2.h5");
    std::string gname = "Obs";
    file2.saveVector(gname, "EigVals", EigVals);
    file2.saveStdVector(gname, "PeakLocation", PeakLocation);
    file2.saveStdVector(gname, "PeakWeight", PeakWeight);
  }
  return 0;
}

