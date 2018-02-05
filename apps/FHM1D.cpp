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
#include "src/Hamiltonian/Hamiltonian.hpp"
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

void LoadEqmParameters( const std::string filename, int &L, int &OBC, int &N1, int &N2, std::vector<RealType> &Jls, std::vector<RealType> &Uls, std::vector<RealType> &Vls){
  HDF5IO file(filename);
  file.LoadNumber("Parameters", "L", L);
  file.LoadNumber("Parameters", "OBC", OBC);
  file.LoadNumber("Parameters", "N1", N1);
  file.LoadNumber("Parameters", "N2", N2);
  file.LoadStdVector("Parameters", "J", Jls);
  file.LoadStdVector("Parameters", "U", Uls);
  file.LoadStdVector("Parameters", "V", Vls);
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
  size_t NewHilbertSpace = NewBases.at(0).getHilbertSpace() * NewBases.at(1).getHilbertSpace();
  ComplexVectorType Vout(NewHilbertSpace, arma::fill::zeros);
  std::vector<int> OldFup = OldBases.at(0).getFStates();
  std::vector<int> OldFdn = OldBases.at(1).getFStates();
  std::vector<size_t> OldFupTag = OldBases.at(0).getFTags();
  std::vector<size_t> OldFdnTag = OldBases.at(1).getFTags();
  std::vector<size_t> NewFupTag = NewBases.at(0).getFTags();
  std::vector<size_t> NewFdnTag = NewBases.at(1).getFTags();
  size_t NewFupIdx, NewFdnIdx, OldFupIdx, OldFdnIdx;
  // int Update = false;
  for ( auto OldFupState : OldFup ){
    if (Species == 1){
      NewFupIdx = NewFupTag.at(OldFupState);// Find their indices
      OldFupIdx = OldFupTag.at(OldFupState);// Find their indices
    }else if ( Species == 0 && Type == 1 && !(btest(OldFupState, CoreHole)) ){// c^\dagger_up
      int NewFupState = ibset(OldFupState, CoreHole);
      NewFupIdx = NewFupTag.at(NewFupState);// Find their indices
      OldFupIdx = OldFupTag.at(OldFupState);// Find their indices
    }else if ( Species == 0 && Type == -1 && btest(OldFupState, CoreHole) ){// c_up
      int NewFupState = ibclr(OldFupState, CoreHole);
      NewFupIdx = NewFupTag.at(NewFupState);// Find their indices
      OldFupIdx = OldFupTag.at(OldFupState);// Find their indices
    }else{
      continue;// next OldFupState
    }
    for ( auto OldFdnState : OldFdn ){
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
      Vout(nid) = Vin(oid);
    }
  }
  return Vout;
}

DTM SingleParticleDensityMatrix( const int species, const std::vector<Basis> &Bases, const DTV &Vec, Hamiltonian<DT> &Ham0 ){
  size_t L = Bases.at(species).getL();
  DTM CM(L, L, arma::fill::zeros);
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

std::vector<DTV> Ni( const std::vector<Basis> &Bases, const DTV &Vec, Hamiltonian<DT> &Ham0 ){
  std::vector<DTV> out;
  DTV tmp1(Bases.at(0).getL(), arma::fill::zeros);//(Bases.at(0).getL(), 0.0e0);
  DTV tmp2(Bases.at(1).getL(), arma::fill::zeros);//(Bases.at(1).getL(), 0.0e0);
  std::vector<int> f1 = Bases.at(0).getFStates();
  std::vector<int> f2 = Bases.at(1).getFStates();
  size_t f1id = 0, f2id = 0;
  for ( int &nf2 : f2 ){
    std::vector<size_t> ids(2,f2id);
    f1id = 0;
    for ( int &nf1 : f1 ){
      ids.at(0) = f1id;
      size_t id = Ham0.DetermineTotalIndex(ids);
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

DTM NiNj( const std::vector<Basis> &Bases, const DTV &Vec, Hamiltonian<DT> &Ham0, const int species1, const int species2 ){
  DTM out(Bases.at(0).getL(), Bases.at(1).getL(), arma::fill::zeros );
  std::vector<int> f1 = Bases.at(0).getFStates();
  std::vector<int> f2 = Bases.at(1).getFStates();
  size_t f1id = 0, f2id = 0;
  for ( int &nf2 : f2 ){
    std::vector<size_t> ids(2,f2id);
    f1id = 0;
    for ( int &nf1 : f1 ){
      ids.at(0) = f1id;
      size_t id = Ham0.DetermineTotalIndex(ids);
      int Ns1 = (species1 == 1)? nf1 : nf2;
      int Ns2 = (species2 == 1)? nf1 : nf2;
      for (size_t cnt1 = 0; cnt1 < Bases.at(0).getL(); cnt1++) {
        for (size_t cnt2 = 0; cnt2 < Bases.at(1).getL(); cnt2++) {
          if ( btest(Ns1, cnt1) && btest(Ns2, cnt2) ) {
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

void Equilibrium(const std::string prefix){
  std::ofstream LogOut;
  LogOut.open(prefix + "fhm.1d.eqm", std::ios::app);
  int L;
  int OBC;
  int N1, N2;
  std::vector<RealType> Jin, Uin, Vin;
  try{
    /* Load parameters from file */
    H5::Exception::dontPrint();
    H5::H5File::isHdf5("conf.h5");
    LoadEqmParameters( "conf.h5", L, OBC, N1, N2, Jin, Uin, Vin);
  }catch(H5::FileIException){
    L = 4;
    OBC = 1;
    N1 = 2;
    N2 = 2;
    Jin = std::vector<RealType>(L-1, 1.0);// OBC
    Uin = std::vector<RealType>(L, 1.0);
    Vin = std::vector<RealType>(L, 0.0);
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
  Hamiltonian<DT> Ham0( Bases );
  // Potential
  std::vector<DT> Vtmp(Vin.begin(), Vin.end());
  std::vector< std::vector<DT> > Vloc = vec(Vtmp, Vtmp);
  // Interaction
  std::vector<DT> Uloc(Uin.begin(), Uin.end());
  Ham0.FermiHubbardModel(Bases, lattice, Vloc, Uloc);
  Ham0.CheckHermitian();
  LogOut << Ham0.GetTotalHilbertSpace() << " DONE!" << std::endl;
  LogOut << "Diagonalize Hamiltonian - " << std::flush;
  RealVectorType Vals;
  DTM Vecs;
  Ham0.eigh(Vals, Vecs, 2);
  // Ham0.diag(Vals, Vecs);// Full spectrum
  LogOut << "DONE!" << std::endl;
  LogOut << "\tGS energy = " << Vals[0] << std::endl;
  LogOut << "\tFES energy = " << Vals[1] << std::endl;
  DTV Vec = Vecs.col(0);
  // ComplexSparseMatrixType H0 = Ham0.GetTotalHamiltonian();
  // std::cout << arma::cdot(Vec, H0 * Vec) << std::endl;
  file->SaveVector("GS", "EVec", Vec);
  file->SaveVector("GS", "EVal", Vals);
  std::vector<DTV> Nfi = Ni( Bases, Vec, Ham0 );
  LogOut << " Up Spin - " << std::endl;
  LogOut << Nfi.at(0) << std::endl;
  DT NupT = arma::sum(Nfi.at(0));
  LogOut << "Total N_up = " << NupT << std::endl;
  LogOut << " Down Spin - " << std::endl;
  LogOut << Nfi.at(1) << std::endl;
  DT NdnT = arma::sum(Nfi.at(1));
  LogOut << "Total N_down = " << NdnT << std::endl;
  LogOut << " N_i - " << std::endl;
  DTV Niall = Nfi.at(0) + Nfi.at(1);
  LogOut << Niall << std::endl;
  DTM NupNdn = NiNj( Bases, Vec, Ham0, 0, 1 );
  LogOut << " Correlation NupNdn" << std::endl;
  LogOut << NupNdn << std::endl;
  DTM NupNup = NiNj( Bases, Vec, Ham0, 0, 0 );
  LogOut << " Correlation NupNup" << std::endl;
  LogOut << NupNup << std::endl;
  DTM NdnNdn = NiNj( Bases, Vec, Ham0, 1, 1 );
  LogOut << " Correlation NdnNdn" << std::endl;
  LogOut << NdnNdn << std::endl;
  LogOut << " N_i^2 - " << std::endl;
  DTM Ni2 = NupNup.diag() + NdnNdn.diag() + 2.0e0 * NupNdn.diag();
  LogOut << Ni2 << std::endl;
  DTM CMUp = SingleParticleDensityMatrix( 0, Bases, Vec, Ham0 );
  LogOut << CMUp << std::endl;
  DTM CMDn = SingleParticleDensityMatrix( 1, Bases, Vec, Ham0 );
  LogOut << CMDn << std::endl;
  file->SaveVector("Obs", "Nup", Nfi.at(0));
  file->SaveNumber("Obs", "NupT", NupT);
  file->SaveVector("Obs", "Ndn", Nfi.at(1));
  file->SaveNumber("Obs", "NdnT", NdnT);
  file->SaveMatrix("Obs", "NupNdn", NupNdn);
  file->SaveMatrix("Obs", "NupNup", NupNup);
  file->SaveMatrix("Obs", "NdnNdn", NdnNdn);
  file->SaveMatrix("Obs", "CMUp", CMUp);
  file->SaveMatrix("Obs", "CMDn", CMDn);
  delete file;
  LogOut.close();
}

void PumpDynamics(const std::string prefix, const int MeasureEvery = 10, const int SaveWFEvery = 1000000 ){
  std::ofstream LogOut;
  LogOut.open(prefix + "pump.fhm.1d", std::ios::app);
  int L;
  int OBC;
  int N1, N2;
  std::vector<RealType> Jeqm, Ueqm, Veqm, At;
  int TSteps;
  RealType dt;
  try{
    /* Load parameters from file */
    H5::Exception::dontPrint();
    H5::H5File::isHdf5("conf.h5");
    LoadEqmParameters( "conf.h5", L, OBC, N1, N2, Jeqm, Ueqm, Veqm);
    LoadPumpParameters( "conf.h5", At, TSteps, dt);
  }catch(H5::FileIException){
    L = 4;
    OBC = 1;
    N1 = 2;
    N2 = 2;
    Jeqm = std::vector<RealType>(L-1, 1.0);// OBC
    Ueqm = std::vector<RealType>(L, 1.0);
    Veqm = std::vector<RealType>(L, 0.0);
    TSteps = 10;
    dt = 0.005;
    At = std::vector<RealType>(TSteps, 0.25 * PI);
  }
  LogOut << "Build Eqm Lattice - " << std::flush;
  std::vector<ComplexType> J(Jeqm.begin(), Jeqm.end());
  std::vector< Node<ComplexType>* > EqmLattice = NN_1D_Chain(L, J, OBC);
  LogOut << "DONE!" << std::endl;
  LogOut << "Build Eqm Basis - " << std::flush;
  Basis F1(L, N1, true);
  F1.Fermion();
  Basis F2(L, N2, true);
  F2.Fermion();
  std::vector<Basis> Bases;
  Bases.push_back(F1);
  Bases.push_back(F2);
  LogOut << "DONE!" << std::endl;
  LogOut << "Build Eqm Hamiltonian - " << std::flush;
  Hamiltonian<ComplexType> Ham0( Bases );
  // Potential
  std::vector<ComplexType> Vtmp(Veqm.begin(), Veqm.end());
  std::vector< std::vector<ComplexType> > Vloc = vec(Vtmp, Vtmp);
  // Interaction
  std::vector<ComplexType> Uloc(Ueqm.begin(), Ueqm.end());
  Ham0.FermiHubbardModel(Bases, EqmLattice, Vloc, Uloc);
  LogOut << "Hermitian = " << Ham0.CheckHermitian() << ", Hilbert space = " << Ham0.GetTotalHilbertSpace() << ", DONE!" << std::endl;
std::cout << Ham0.GetTotalHamiltonian() << std::endl;
  // Load Wavefunction
  LogOut << "Load Eqm Wavefunction - " << std::flush;
  HDF5IO *EqmFile = new HDF5IO("FHMChainData.h5");
  ComplexVectorType VecInit;
  EqmFile->LoadVector("GS", "EVec", VecInit);
  delete EqmFile;
  LogOut << VecInit.n_rows << " DONE!" << std::endl;

  /* Pumping */
  ComplexType Prefactor = ComplexType(0.0, -1.0e0*dt);/* NOTE: hbar = 1 */
  size_t PSteps = At.size();
  ComplexVectorType VecPump = VecInit;
  HDF5IO* file2 = new HDF5IO("Pump.h5");
  std::string gname = "Obs-0/";
  ComplexType Lecho = arma::cdot(VecInit, VecPump);
  std::vector<DTV> Nfi = Ni( Bases, VecPump, Ham0 );
  file2->SaveNumber(gname, "Lecho", Lecho);
  file2->SaveVector(gname, "Nup", Nfi.at(0));
  file2->SaveVector(gname, "Ndn", Nfi.at(1));
  delete file2;
  LogOut << "Begin pumping......" << std::endl;
  for (size_t cntP = 1; cntP <= PSteps; cntP++) {
    if ( OBC ){
      J = std::vector<ComplexType>(L - 1, exp( ComplexType(0.0, 1.0) * At.at(cntP-1)) );
    } else{
      J = std::vector<ComplexType>(L, exp( ComplexType(0.0, 1.0) * At.at(cntP-1)) );
    }
    std::vector< Node<ComplexType>* > PLattice = NN_1D_Chain(L, J, OBC);
    // Potential
    Vtmp = std::vector<ComplexType>(Veqm.begin(), Veqm.end());
    std::vector< std::vector<ComplexType> > Vt = vec(Vtmp, Vtmp);
    // Interaction
    std::vector<ComplexType> Ut(Ueqm.begin(), Ueqm.end());
    Ham0.FermiHubbardModel(Bases, PLattice, Vt, Ut);
std::cout << Ham0.GetTotalHamiltonian() << std::endl;
std::cin.get();
    if ( Ham0.CheckHermitian() ) LogOut << "non-Hermitian Hamiltonian at PStep = " << cntP << std::endl;
    // Evolve the state
    Ham0.expH(Prefactor, VecPump);
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
  Ham0.FermiHubbardModel(Bases, EqmLattice, Vloc, Uloc);
  ComplexSparseMatrixType H0 = Ham0.GetTotalHamiltonian();
  ComplexType Energy = arma::cdot(VecPump, H0 * VecPump);
  HDF5IO* file3 = new HDF5IO("PumpWF.h5");
  gname = "WF";
  file3->SaveVector(gname, "Vec", VecPump);
  file3->SaveNumber(gname, "Energy", Energy);
  delete file3;
  LogOut << "Finished pumping!!" << std::endl;
}

void XASDynamics(const std::string prefix, const int MeasureEvery = 2, const int SaveWFEvery = 1000000 ){
  std::ofstream LogOut;
  LogOut.open(prefix + "xas.fhm.1d", std::ios::app);
  int L;
  int OBC;
  int N1, N2;
  std::vector<RealType> Jeqm, Ueqm, Veqm, Uch, Vch;
  int TSteps;
  RealType dt;
  int CoreHole, Species, Type;
  try{
    /* Load parameters from file */
    H5::Exception::dontPrint();
    H5::H5File::isHdf5(prefix + "conf.h5");
    LoadEqmParameters( prefix + "conf.h5", L, OBC, N1, N2, Jeqm, Ueqm, Veqm);
    LoadXASParameters( prefix + "conf.h5", Uch, Vch, TSteps, dt, CoreHole, Species, Type);
  }catch(H5::FileIException){
    L = 12;
    OBC = 1;
    N1 = 6;
    N2 = 6;
    Jeqm = std::vector<RealType>(L-1, 1.0);// OBC
    Uch = std::vector<RealType>(L, 1.0);
    Vch = std::vector<RealType>(L, 0.0);
    CoreHole = 6;
    Vch.at(CoreHole) = -3.0;
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
  Hamiltonian<ComplexType> EqmHam( EqmBases );
  // Potential
  std::vector<ComplexType> Vw(Veqm.begin(), Veqm.end());
  std::vector< std::vector<ComplexType> > EqmVloc = vec(Vw, Vw);
  // Interaction
  std::vector<ComplexType> EqmUloc(Ueqm.begin(), Ueqm.end());
  EqmHam.FermiHubbardModel(EqmBases, Lattice, EqmVloc, EqmUloc);
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
  Hamiltonian<ComplexType> CoreHoleHam( CoreHoleBases );
  // Potential
  std::vector<ComplexType> CoreHoleVw(Vch.begin(), Vch.end());
  std::vector< std::vector<ComplexType> > CoreHoleVloc = vec(CoreHoleVw, CoreHoleVw);
  // Interaction
  std::vector<ComplexType> CoreHoleUloc(Uch.begin(), Uch.end());
  CoreHoleHam.FermiHubbardModel(CoreHoleBases, Lattice, CoreHoleVloc, CoreHoleUloc);
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
  LogOut << "Operate on Wave function with core hole @ " << CoreHole << ", on species - " << Species << ", type - " << Type << std::endl;
  ComplexVectorType VecInit = Operate( VecInput, CoreHole, Species, Type, EqmBases, CoreHoleBases, EqmHam, CoreHoleHam );
  LogOut << " DONE!" << std::endl;
  LogOut.close();
}

/* main program */
int main(int argc, char *argv[]){
  if ( argc < 2 ) RUNTIME_ERROR(" Need at least one argument to run program. Use 0 to run Equilibrium.");
#if defined(MKL)
  mkl_set_dynamic(0);
  mkl_set_num_threads(NumCores);
#endif
  if ( std::atoi(argv[1]) == 0 ){
    Equilibrium("");
  }else if ( std::atoi(argv[1]) == 1 ){
    PumpDynamics("");
  // }else if ( std::atoi(argv[1]) == 2 ){
    // XASDynamics("");
  }else if ( std::atoi(argv[1]) == 10 ){
    Equilibrium("");
    PumpDynamics("");
  }//else if ( std::atoi(argv[1]) == 20 ){
  //   Equilibrium("");
  //   XASDynamics("");
  // }else if ( std::atoi(argv[1]) == 210 ){
  //   Equilibrium("");
  //   PumpDynamics("");
  //   XASDynamics("");
  // }
  return 0;
}
