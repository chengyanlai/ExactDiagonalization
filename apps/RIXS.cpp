#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <stdlib.h>// atof
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
#ifdef MPIPARALLEL
  #include <mpi.h>
#endif

#ifndef DEBUG
  #define DEBUG 1
#endif

#ifndef NumCores
  #define NumCores 20
#endif

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

std::vector< RealVectorType > Ni( const std::vector<Basis> &Bases, const ComplexVectorType &Vec,
  Hamiltonian<ComplexType> &ham ){
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
  int L = std::atoi(argv[1]);
  const int OBC = 0;
  int N1 = std::atoi(argv[2]);
  int N2 = N1;
  int S2 = std::atoi(argv[3]);;
  int Tsteps = std::atoi(argv[4]);
  const RealType dt = 0.005;
  RealType Uinit  = std::atof(argv[5]);
  RealType Vin = std::atof(argv[6]);
  RealType Uin = std::atof(argv[7]);
#ifdef MPIPARALLEL
  // Initialize MPI
  MPI_Init(NULL, NULL);
  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
#else
  int world_size = 1;
  // int world_rank = 196;
  for ( size_t world_rank = 0; world_rank < L; world_rank++ ){
#endif
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
  std::vector<int> st1 = F1.GetFStates();
  std::vector<size_t> tg1 = F1.GetFTags();
  Basis F2(L, N2, true);
  F2.Fermion();
  std::vector<int> st2 = F2.GetFStates();
  std::vector<size_t> tg2 = F2.GetFTags();
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
  ComplexVectorType Vec = Vecs.row(0);
  INFO("GS energy = " << Vals[0]);
  INFO("DONE!");
  /* Density calculation */
  std::vector< RealVectorType > Nfi = Ni( Bases, Vec, ham );
  RealVectorType Niall = Nfi.at(0) + Nfi.at(1);

  int CHloc = world_rank % L;
  std::cout << "Core-hole at " << CHloc << std::endl;

  /* Save the initial G.S. wavefunction */
  std::string filename1 = "WF-";
  filename1.append(std::to_string((unsigned long long)CHloc));
  filename1.append(".h5");
  HDF5IO *file = new HDF5IO(filename1);
  file->saveVector("GS", "EVec", Vec);
  file->saveVector("GS", "EVal", Vals);
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
  Vtmp.at(CHloc) = Vin;
  Vloc.clear();
  Vloc.push_back(Vtmp);
  Vloc.push_back(Vtmp);
  Utmp.at(CHloc) = Uin;
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

  std::string gname = "t-0";
  file->saveVector(gname, "wf", VecT);
  file->saveVector(gname, "Nup", Nfi.at(0));
  file->saveVector(gname, "Ndn", Nfi.at(1));
  delete file;

  std::cout << "Time evolution - T " << std::endl;
  std::cout << " H dim = " << nHam.GetTotalHilbertSpace() << " WF dim = " << VecT.rows() << std::endl;
  ComplexType Prefactor = ComplexType(0.0, -1.0e0*dt);/* NOTE: hbar = 1 */
  for (size_t cntT = 1; cntT <= Tsteps; cntT++) {
    nHam.expH(Prefactor, VecT);
      if ( cntT % 20 == 0 ){
        std::string filename3 = "XASDYN-";
        filename3.append( std::to_string((unsigned long long)CHloc) );
        filename3.append(".h5");
        HDF5IO file3(filename3);
        std::string gname = "Obs-";
        gname.append(std::to_string((unsigned long long)cntT));
        gname.append("/");
        ComplexType Lecho = VecInit.dot(VecT);
        Nfi = Ni( nBases, VecT, nHam );
        file3.saveNumber(gname, "Lecho", Lecho);
        file3.saveVector(gname, "Nup", Nfi.at(0));
        file3.saveVector(gname, "Ndn", Nfi.at(1));
      }
  }

  if ( S2 ){
    std::cout << "New basis required due to spin flip" << std::endl;
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

  std::cout << "Destroy core-hole state" << std::endl;
  ComplexVectorType VecS = OperateC( nBases, VecT, Bases, CHloc, S2, nHam, ham);
  VecS.normalize();
  Nfi = Ni( Bases, VecS, ham );
  Niall = Nfi.at(0) + Nfi.at(1);

  std::string filename2 = "SData-";
  filename2.append(std::to_string((unsigned long long)CHloc));
  filename2.append(".h5");
  HDF5IO *file2 = new HDF5IO(filename2);
  gname = "s-0";
  file2->saveVector(gname, "wf", VecS);
  file2->saveVector(gname, "Nup", Nfi.at(0));
  file2->saveVector(gname, "Ndn", Nfi.at(1));
  delete file2;

  std::cout << "Time evolution - S " << std::endl;
  std::cout << " H dim = " << ham.GetTotalHilbertSpace() << " WF dim = " << VecS.rows() << std::endl;
  Prefactor = ComplexType(0.0, -1.0e0*dt);/* NOTE: hbar = 1 */
  for (size_t cntT = 1; cntT <= 3000; cntT++) {
    ham.expH(Prefactor, VecS);
    if ( cntT % 20 == 0 ){
      file2 = new HDF5IO(filename2);
      gname = "s-";
      gname.append(std::to_string((unsigned long long)cntT));
      Nfi = Ni( Bases, VecS, ham );
      file2->saveVector(gname, "wf", VecS);
      file2->saveVector(gname, "Nup", Nfi.at(0));
      file2->saveVector(gname, "Ndn", Nfi.at(1));
      delete file2;
    }
  }
#ifdef MPIPARALLEL
  MPI_Finalize();
#else
  }
#endif
  return 0;
}

