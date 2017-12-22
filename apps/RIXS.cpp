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

#include "mkl.h"

#include <mpi.h>

#ifndef DEBUG
  #define DEBUG 1
#endif

#ifndef NumCores
  #define NumCores 20
#endif

void LoadParameters( const std::string filename, int &L, int &N1, int &N2, RealType &Uinit, RealType &Uloc, RealType &Vloc, int &Tsteps, RealType &dT, int& S2){
    HDF5IO file(filename);
    L = file.loadInt("Parameters", "L");
    N1 = file.loadInt("Parameters", "N1");
    N2 = file.loadInt("Parameters", "N2");
    Uinit = file.loadReal("Parameters", "Uinit");
    Uloc = file.loadReal("Parameters", "U");
    Vloc = file.loadReal("Parameters", "V");
    Tsteps = file.loadInt("Parameters", "Tsteps");
    dT = file.loadReal("Parameters", "dt");
    S2 = file.loadInt("Parameters", "S2");
}

ComplexVectorType OperateCdagger( const std::vector<Basis> OldBases, const ComplexVectorType Vin,
  const std::vector<Basis> NewBases, const int CHloc, const int species,
  Hamiltonian<ComplexType> OldHam, Hamiltonian<ComplexType> NewHam  ){
  size_t NewHilbertSpace = NewBases.at(0).getHilbertSpace() * NewBases.at(1).getHilbertSpace();
  ComplexVectorType Vout = ComplexVectorType::Zero(NewHilbertSpace);
  std::vector<int> OldFup = OldBases.at(0).getFStates();
  std::vector<int> OldFdn = OldBases.at(1).getFStates();
  std::vector<size_t> OldFupTag = OldBases.at(0).getFTags();
  std::vector<size_t> OldFdnTag = OldBases.at(1).getFTags();
  std::vector<size_t> NewFupTag = NewBases.at(0).getFTags();
  std::vector<size_t> NewFdnTag = NewBases.at(1).getFTags();
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
  size_t NewHilbertSpace = NewBases.at(0).getHilbertSpace() * NewBases.at(1).getHilbertSpace();
  ComplexVectorType Vout = ComplexVectorType::Zero(NewHilbertSpace);
  std::vector<int> OldFup = OldBases.at(0).getFStates();
  std::vector<int> OldFdn = OldBases.at(1).getFStates();
  std::vector<size_t> OldFupTag = OldBases.at(0).getFTags();
  std::vector<size_t> OldFdnTag = OldBases.at(1).getFTags();
  std::vector<size_t> NewFupTag = NewBases.at(0).getFTags();
  std::vector<size_t> NewFdnTag = NewBases.at(1).getFTags();
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

std::vector< RealVectorType > Ni( const std::vector<Basis> &Bases, const ComplexVectorType &Vec,
  Hamiltonian<ComplexType> &ham ){
  std::vector< RealVectorType > out;
  RealVectorType tmp1 = RealVectorType::Zero(Bases.at(0).getL());//(Bases.at(0).getL(), 0.0e0);
  RealVectorType tmp2 = RealVectorType::Zero(Bases.at(1).getL());//(Bases.at(1).getL(), 0.0e0);
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

int main(int argc, char const *argv[]) {
  Eigen::setNbThreads(NumCores);
  mkl_set_num_threads(NumCores);
  // Initialize MPI
  MPI_Init(NULL, NULL);
  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  INFO("Eigen3 uses " << Eigen::nbThreads() << " threads.");
  int L, CHloc;
  const int OBC = 0;
  int N1, N2, S2;
  int Tsteps;
  RealType dt;
  RealType Uinit, Vin, Uin;
  LoadParameters( "conf.h5", L, N1, N2, Uinit, Uin, Vin, Tsteps, dt, S2);
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
  for ( auto &lt : lattice ){
    if ( !(lt->VerifySite()) ) RUNTIME_ERROR("Wrong lattice setup!");
  }
  INFO("DONE!");
  INFO("Build Basis - ");
  Basis F1(L, N1, true);
  F1.Fermion();
  std::vector<int> st1 = F1.getFStates();
  std::vector<size_t> tg1 = F1.getFTags();
  Basis F2(L, N2, true);
  F2.Fermion();
  std::vector<int> st2 = F2.getFStates();
  std::vector<size_t> tg2 = F2.getFTags();
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
  std::vector<ComplexType> Utmp(L, ComplexType(Uinit, 0.0e0) );
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
  ham.eigh(Vals, Vecs);
  ComplexVectorType Vec = Vecs.col(0);
  INFO("GS energy = " << Vals[0]);
  INFO("DONE!");


  int CHloc1 = world_rank / L;
  int CHloc2 = world_rank % L;


  std::vector< RealVectorType > Nfi = Ni( Bases, Vec, ham );
  RealVectorType Niall = Nfi.at(0) + Nfi.at(1);
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
  Vtmp.at(CHloc1) = Vin;
  Vloc.clear();
  Vloc.push_back(Vtmp);
  Vloc.push_back(Vtmp);
  Utmp.at(CHloc1) = Uin;
  Uloc.clear();
  Uloc.push_back(Utmp);
  Uloc.push_back(Utmp);
  nHam.BuildLocalHamiltonian(Vloc, Uloc, nBases);
  nHam.BuildHoppingHamiltonian(nBases, lattice);
  nHam.BuildTotalHamiltonian();
  /* Apply the operator */
  std::cout << "Create core-hole state" << std::endl;

  std::string filename1 = "WF";
  filename1.append(std::to_string((unsigned long long)CHloc1));
  filename1.append("-");
  filename1.append(std::to_string((unsigned long long)CHloc2));
  filename1.append(".h5");
  HDF5IO *file = new HDF5IO(filename1);
  file->saveVector("GS", "EVec", Vec);
  file->saveVector("GS", "EVal", Vals);
  file->saveVector("Obs", "Nup", Nfi.at(0));
  file->saveVector("Obs", "Ndn", Nfi.at(1));

  ComplexVectorType VecInit = OperateCdagger( Bases, Vec, nBases, CHloc1, 0, ham, nHam);
  VecInit.normalize();
  ComplexVectorType VecT = VecInit;
  Nfi = Ni( nBases, VecT, nHam );
  Niall = Nfi.at(0) + Nfi.at(1);

  std::string gname = "t-0";
  file->saveVector(gname, "wf", VecT);
  file->saveVector(gname, "Nup", Nfi.at(0));
  file->saveVector(gname, "Ndn", Nfi.at(1));

  ComplexType Prefactor = ComplexType(0.0, -1.0e0*dt);/* NOTE: hbar = 1 */
  for (size_t cntT = 1; cntT <= Tsteps; cntT++) {
    nHam.expH(Prefactor, VecT);
    // if ( cntT % 20 == 0 ){
    //   gname = "t-";
    //   gname.append(std::to_string((unsigned long long)cntT));
    //   Nfi = Ni( nBases, VecT, nHam );
    //   file->saveVector(gname, "wf", VecT);
    //   file->saveVector(gname, "Nup", Nfi.at(0));
    //   file->saveVector(gname, "Ndn", Nfi.at(1));
    // }
  }

  if ( S2 ){
    /* Build New Basis */
    Bases.clear();
    Basis sF1(L, N1+1, true);
    sF1.Fermion();
    Basis sF2(L, N2-1, true);
    sF2.Fermion();
    Bases.push_back(sF1);
    Bases.push_back(sF2);
    /* Update the new Hamiltonian */
    Vtmp = std::vector<ComplexType>(L, ComplexType(0.0e0, 0.0e0) );
    Vloc.clear();
    Vloc.push_back(Vtmp);
    Vloc.push_back(Vtmp);
    Utmp = std::vector<ComplexType>(L, ComplexType(Uinit, 0.0e0) );
    Uloc.clear();
    Uloc.push_back(Utmp);
    Uloc.push_back(Utmp);
    ham.BuildLocalHamiltonian(Vloc, Uloc, Bases);
    ham.BuildHoppingHamiltonian(Bases, lattice);
    ham.BuildTotalHamiltonian();
  }

  ComplexVectorType VecS = OperateC( nBases, VecT, Bases, CHloc2, S2, nHam, ham);
  VecS.normalize();
  Nfi = Ni( nBases, VecS, nHam );
  Niall = Nfi.at(0) + Nfi.at(1);

  gname = "s-0";
  file->saveVector(gname, "wf", VecS);
  file->saveVector(gname, "Nup", Nfi.at(0));
  file->saveVector(gname, "Ndn", Nfi.at(1));

  Prefactor = ComplexType(0.0, -1.0e0*dt);/* NOTE: hbar = 1 */
  for (size_t cntT = 1; cntT <= 3000; cntT++) {
    ham.expH(Prefactor, VecS);
    if ( cntT % 20 == 0 ){
      gname = "s-";
      gname.append(std::to_string((unsigned long long)cntT));
      Nfi = Ni( Bases, Vecs, ham );
      file->saveVector(gname, "wf", VecS);
      file->saveVector(gname, "Nup", Nfi.at(0));
      file->saveVector(gname, "Ndn", Nfi.at(1));
    }
  }

  delete file;
  return 0;
}

