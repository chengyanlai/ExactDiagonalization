#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <iterator>
#include "src/EDType.hpp"
#include "src/Node/Node.hpp"
#include "src/Lattice/preset.hpp"
#include "src/Basis/Basis.hpp"
#include "src/Hamiltonian/Holstein/Holstein.hpp"
#include "src/hdf5io/hdf5io.hpp"
#include "src/numeric/lapack.h"

#ifdef MKL
  #include "mkl.h"
#endif

#if not defined(NumCores)
  #define NumCores 16
#endif

#ifdef MPIPARALLEL
  #include <mpi.h>
#endif

void LoadEqmParameters( const std::string filename, int& L, int& N, RealType& G, RealType& W){
  HDF5IO h5f(filename);
  h5f.LoadNumber("Parameters", "L", L);
  h5f.LoadNumber("Parameters", "N", N);
  h5f.LoadNumber("Parameters", "G", G);
  h5f.LoadNumber("Parameters", "W", W);
}

ComplexVectorType NPhonon( const std::vector<Basis> &Bases, const ComplexVectorType &Vec, const Holstein<ComplexType>& Ham){
  int L = Bases.at(0).GetL();
  ComplexVectorType out(L, arma::fill::zeros);
  std::vector< std::vector<int> > b = Bases.at(0).GetBStates();
  assert( b.size() * L == Vec.size() );
  for ( size_t cnt = 0; cnt < L; cnt++ ){// phonon r
    int coff = 0;
    for ( auto &nbi : b ){
      for (size_t cnt2 = 0; cnt2 < L; cnt2++) {// fermion r
        size_t idx = Ham.DetermineTotalIndex( vec<size_t>(cnt2, coff) );
        out.at(cnt) += (RealType)nbi.at(cnt) * std::pow(std::abs(Vec(idx)), 2);
      }
      coff++;
    }
  }
  return out;
}

ComplexMatrixType NFermion( const std::vector<Basis> &Bases, const ComplexVectorType &Vec, const Holstein<ComplexType>& Ham){
  int L = Bases.at(0).GetL();
  ComplexMatrixType out(L, L, arma::fill::zeros);
  assert( Bases.at(0).GetHilbertSpace() * L == Vec.size() );
  for ( size_t s1 = 0; s1 < L; s1++ ){
    for ( size_t s2 = 0; s2 < L; s2++ ){
      for ( size_t b = 0; b < Bases.at(0).GetHilbertSpace(); b++ ){
        size_t idx1 = Ham.DetermineTotalIndex( vec<size_t>(s1, b) );
        size_t idx2 = Ham.DetermineTotalIndex( vec<size_t>(s2, b) );
        out.at(s1, s2) += Vec(idx1) * Conjg(Vec(idx2));
      }
    }
  }
  return out;
}

ComplexVectorType JPsi( const std::vector<Basis> &Bases, const ComplexVectorType &Vec, const Holstein<ComplexType>& Ham){
  int L = Bases.at(0).GetL();
  ComplexVectorType out(Vec.size(), arma::fill::zeros);
  for ( size_t left = 0; left < L; left++ ){
    size_t right = ( left == L - 1 ) ? 0 : left + 1;
    for ( size_t b = 0; b < Bases.at(0).GetHilbertSpace(); b++ ){
      size_t idx1 = Ham.DetermineTotalIndex( vec<size_t>(left, b) );
      size_t idx2 = Ham.DetermineTotalIndex( vec<size_t>(right, b) );
      out(idx1) += ComplexType(0.0, 1.0) * Vec(idx2);
      out(idx2) += ComplexType(0.0,-1.0) * Vec(idx1);
    }
  }
  return out;
}

void Equilibrium(const std::string prefix, int NEV){
  std::ofstream LogOut;
  LogOut.open(prefix + "Holstein.R.eqm", std::ios::app);
  int L = 4;
  int N = 5 * L;
  const RealType Jin = 1.0;
  RealType Win = 0.50;
  RealType Gin = 1.0;

  try{
    H5::Exception::dontPrint();
    H5::H5File::isHdf5(prefix + "conf.h5");
    LoadEqmParameters( prefix + "conf.h5", L, N, Gin, Win);
  }catch(H5::FileIException){
    LogOut << "Use default settings." << std::endl;
  }

  LogOut << "Build Basis - max. phonon quanta = " << N << ", ends up with Phonon Hilbert space = " << std::flush;
  Basis P1(L, N);
  P1.PhononR();
  LogOut << P1.GetHilbertSpace() << std::flush;
  // LogOut << P1 << std::endl;
  std::vector<Basis> Bases;
  Bases.push_back(P1);
  LogOut << " DONE!" << std::endl;

  LogOut << "Build Hamiltonian - " << std::flush;
  Holstein<ComplexType> Ham0( L, Bases );
  Ham0.HolsteinModelR(L, Bases, Win, Gin, Jin);
  LogOut << "Hermitian = " << Ham0.CheckHermitian() << ", Hilbert space = " << Ham0.GetTotalHilbertSpace() << ", DONE!" << std::endl;

  LogOut << "Diagonalize Hamiltonian to find GS " << std::flush;
  RealVectorType Vals;
  ComplexMatrixType Vecs;
  if ( Ham0.GetTotalHilbertSpace() > 1200 ){
    Ham0.eigh(Vals, Vecs, NEV, true);
  }else{
    Ham0.diag(Vals, Vecs);// Full spectrum
    NEV = Ham0.GetTotalHilbertSpace();
  }
  LogOut << "DONE!" << std::endl;
  // Sort eigenvalues
  std::vector<size_t> IndexOrder(Vals.size());
  std::iota( IndexOrder.begin(), IndexOrder.end(), 0 );
  std::sort( IndexOrder.begin(), IndexOrder.end(), [&Vals](size_t i1, size_t i2) {return Vals[i1] < Vals[i2];});
  // End sorting
  LogOut << "\tGS energy = " << std::setprecision(12) << Vals[0] << std::endl;
  LogOut << "\tFES energy = " << std::setprecision(12) << Vals[1] << std::endl;
  HDF5IO* file = new HDF5IO(prefix + "HolsteinR.h5");
  for ( size_t j = 0; j < IndexOrder.size(); j++ ){
    size_t i = IndexOrder.at(j);
    std::string gname = "S-";
    gname.append( std::to_string( (unsigned long)i) );
    file->SaveNumber(gname, "EVal", Vals[i]);
    ComplexVectorType Vec = Vecs.col(i);
    file->SaveVector(gname, "EVec", Vec);
    ComplexVectorType Npi = NPhonon( Bases, Vec, Ham0);
    file->SaveVector(gname, "Phonon", Npi);
    ComplexMatrixType Nfi = NFermion( Bases, Vec, Ham0);
    file->SaveMatrix(gname, "Fermion", Nfi);
    LogOut << i << ": E = " << Vals[i] << ", Np = " << arma::accu(Npi) << ", Nf = " << arma::trace(Nfi) << std::endl;
  }
  delete file;
  LogOut.close();
}

void LoadDynParameters( const std::string filename, int& L, int& N, RealType& G, RealType& W, int& TSteps, RealType& dt){
  HDF5IO h5f(filename);
  h5f.LoadNumber("Parameters", "L", L);
  h5f.LoadNumber("Parameters", "N", N);
  h5f.LoadNumber("Parameters", "G", G);
  h5f.LoadNumber("Parameters", "W", W);
  h5f.LoadNumber("Parameters", "TSteps", TSteps);
  h5f.LoadNumber("Parameters", "dt", dt);
}

void LoadAlphas(const std::string filename, std::vector<ComplexType>& Alphas){
  HDF5IO h5f(filename);
  std::vector<RealType> Ar, Ap;
  h5f.LoadStdVector("Parameters", "AlphaReal", Ar);
  h5f.LoadStdVector("Parameters", "AlphaPhase", Ap);
  assert( Ar.size() == Ap.size() );
  Alphas.clear();
  for ( size_t i = 0; i < Ar.size(); i++ ){
    Alphas.push_back( Ar.at(i) * exp(ComplexType(0.0, 1.0)*Ap.at(i)) );
  }
}

std::vector<RealType> LocalEnergy(ComplexVectorType& VecIn, const std::vector<Basis>& Bases, const Holstein<ComplexType>& Ham0 ){
  std::vector<RealType> Eis;
  int L = Bases.at(0).GetL();
  std::vector< std::vector<int> > b = Bases.at(0).GetBStates();
  for ( size_t i = 0; i < L; i++ ){
    ComplexVectorType tmp(Ham0.GetTotalHilbertSpace(), arma::fill::zeros);
    for ( size_t j = 0; j < b.size(); j++ ){
      size_t idx = Ham0.DetermineTotalIndex( vec<size_t>(i, j) );
      tmp(idx) = VecIn(idx);
    }
    tmp = arma::normalise(tmp);
    RealType E0 = RealPart( arma::cdot(tmp, Ham0.GetTotalHamiltonian() * tmp) );
    Eis.push_back(E0);
  }
  return Eis;
}

ComplexVectorType CoherentState( std::ofstream& LogOut, const std::vector<ComplexType>& alphas, const std::vector<Basis>& Bases, const Holstein<ComplexType>& Ham0){
  std::vector<ComplexVectorType> Vecs;
  std::vector< std::vector<int> > b = Bases.at(0).GetBStates();
  int L = Bases.at(0).GetL();
  for ( size_t i = 0; i < L; i++ ){
    Vecs.push_back( ComplexVectorType(Ham0.GetTotalHilbertSpace(), arma::fill::zeros) );
  }
  for ( size_t i = 0; i < b.size(); i++ ){
    std::vector<int> nb = b.at(i);
    ComplexType Coff(0.0e0, 0.0e0);
    for ( size_t j = 0; j < nb.size(); j++ ){
      int Nbj = nb.at(j);
      ComplexType CoffTmp(1.0e0, 0.0e0);
      for ( size_t k = 1; k <= Nbj; k++ ){
        ComplexType tmp = alphas.at(j) / std::sqrt(RealType(k));
        CoffTmp *= tmp;
      }
      if ( j == 0 ) Coff = CoffTmp;
      else Coff *= CoffTmp;
    }
    for ( size_t k = 0; k < L; k++ ){
      size_t idx = Ham0.DetermineTotalIndex( vec<size_t>(k, i) );
      Vecs.at(k)(idx) = 1.0e-1 * Coff;
    }
  }
  std::vector<RealType> Eis;
  RealMatrixType ME(L, L, arma::fill::zeros);
  for ( size_t i = 0; i < L; i++ ){
    Vecs.at(i) = arma::normalise(Vecs.at(i));
    RealType E0 = RealPart( arma::cdot(Vecs.at(i), Ham0.GetHCouple() * Vecs.at(i)) );
    Eis.push_back(E0);
    ME.at(i,i) = E0;
    if ( i == L - 1 ){
      ME.at(0, i) = -1.0;
      ME.at(i, 0) = -1.0;
    }else{
      ME.at(i+1, i) = -1.0;
      ME.at(i, i+1) = -1.0;
    }
  }
  PrintVector(LogOut, Eis);
  RealVectorType eigval;
  RealMatrixType eigvec;
  arma::eig_sym(eigval, eigvec, ME);
  ComplexVectorType out(Ham0.GetTotalHilbertSpace(), arma::fill::zeros);
  for ( size_t i = 0; i < L; i++ ) out += eigvec.col(0)[i] * Vecs.at(i);
  return arma::normalise(out);
  // int index = std::distance(Eis.begin(), std::min_element(Eis.begin(), Eis.end()) );
  // return Vecs.at(index);
}

void Dynamics(const std::string prefix, const std::string InitialState, const int S1, const int S2, const int MeasureEvery, const int SaveWFEvery ){
  std::ofstream LogOut;
  LogOut.open(prefix + "Holstein.R.dyn", std::ios::app);
  int L = 4;
  int N = 3 * L;
  const RealType Jin = 1.0;
  RealType Win = 0.60;
  RealType Gin = 1.0;
  int TSteps = 20000;
  RealType dt = 0.005;
  try{
    H5::Exception::dontPrint();
    H5::H5File::isHdf5(prefix + "conf.h5");
    LoadDynParameters( prefix + "conf.h5", L, N, Gin, Win, TSteps, dt);
  }catch(H5::FileIException){
    LogOut << "Use default settings." << std::endl;
  }

  LogOut << "Build Basis - max. phonon quanta = " << N << ", ends up with Phonon Hilbert space = " << std::flush;
  Basis P1(L, N);// Get rid og k=0 phonon mode, so using L - 1!!
  P1.PhononR();
  LogOut << P1.GetHilbertSpace() << std::flush;
  std::vector<Basis> Bases;
  Bases.push_back(P1);
  LogOut << " DONE!" << std::endl;

  LogOut << "Build Hamiltonian - " << std::flush;
  Holstein<ComplexType> Ham0( L, Bases );
  Ham0.HolsteinModelR(L, Bases, Win, Gin, Jin);
  LogOut << "Hermitian = " << Ham0.CheckHermitian() << ", Hilbert space = " << Ham0.GetTotalHilbertSpace() << ", DONE!" << std::endl;

  LogOut << "Load Wavefunction - " << InitialState << " " << std::flush;
  ComplexVectorType VecInit(Ham0.GetTotalHilbertSpace(), arma::fill::randn);
  std::string SaveFile = "QuenchState";
  if ( InitialState == "E" && S1 >= 0 && S2 >= 0 ){
    try{
      H5::Exception::dontPrint();
      H5::H5File::isHdf5(prefix + "HolsteinR.h5");
      HDF5IO* file = new HDF5IO(prefix + "HolsteinR.h5");
      ComplexVectorType V1, V2;
      std::string gname = "S-";
      gname.append( std::to_string((unsigned long)S1) );
      SaveFile.append( "-" );
      SaveFile.append( std::to_string((unsigned long)S1) );
      LogOut << " <" << gname << std::flush;
      file->LoadVector(gname, "EVec", V1);
      gname = "S-";
      gname.append( std::to_string((unsigned long)S2) );
      SaveFile.append( "-" );
      SaveFile.append( std::to_string((unsigned long)S2) );
      LogOut << "|" << gname << std::flush;
      file->LoadVector(gname, "EVec", V2);
      delete file;
      VecInit = ( V1 + V2 );
      LogOut << "> = " << arma::cdot(V1, V2) << " DONE." << std::endl;
    }catch(H5::FileIException){
      RUNTIME_ERROR("Can not load eigensate from file - HolsteinR.h5. ");
    }
  }else if ( InitialState == "Z" ){
    SaveFile.append( "-Z-" );
    SaveFile.append( std::to_string(S1) );
    /* Target phonon mode - coherent state */
    std::vector<ComplexType> alphas;
    try{
      H5::Exception::dontPrint();
      H5::H5File::isHdf5(prefix + "conf.h5");
      LoadAlphas( prefix + "conf.h5", alphas );
    }catch(H5::FileIException){
      alphas.push_back( 2.0 * exp(0.00*PI*ComplexType(0.0,1.0)) );
      alphas.push_back( 2.0 * exp(0.25*PI*ComplexType(0.0,1.0)) );
      alphas.push_back( 2.0 * exp(0.50*PI*ComplexType(0.0,1.0)) );
      alphas.push_back( 2.0 * exp(0.75*PI*ComplexType(0.0,1.0)) );
      LogOut << "Use default settings." << std::endl;
    }
    LogOut << "Fermion is at GS, and phonon coherent state: alpha_i = " << std::endl;
    PrintVector(LogOut, alphas);
    LogOut << "Local effective potential energy = " << std::endl;
    VecInit = CoherentState(LogOut, alphas, Bases, Ham0);
    LogOut << " DONE!" << std::endl;
  }else{
    SaveFile.append( "-R" );
    LogOut << "Use random initial - " << VecInit.n_rows << std::endl;
  }
  VecInit = arma::normalise(VecInit);
  LogOut << "Energy of this initial state E = " <<  arma::cdot(VecInit, Ham0.GetTotalHamiltonian() * VecInit) << std::endl;

  ComplexType Prefactor = ComplexType(0.0, -1.0e0*dt);/* NOTE: hbar = 1 */
  ComplexVectorType VecDyn = VecInit;
  HDF5IO* file2 = new HDF5IO(prefix + SaveFile + ".h5");
  std::string gname = "Obs-0/";
  ComplexType Lecho = arma::cdot(VecInit, VecDyn);
  file2->SaveNumber(gname, "F", Lecho);
  ComplexVectorType JVec = JPsi(Bases, VecInit, Ham0);
  ComplexVectorType JVecDyn = JPsi(Bases, VecDyn, Ham0);
  ComplexType JJt = arma::cdot(JVecDyn, JVec);
  file2->SaveNumber(gname, "JJt", JJt);
  ComplexVectorType Npi = NPhonon( Bases, VecDyn, Ham0);
  file2->SaveVector(gname, "Phonon", Npi);
  ComplexMatrixType Nfi = NFermion( Bases, VecDyn, Ham0);
  file2->SaveMatrix(gname, "Fermion", Nfi);
  std::vector<RealType> Ei = LocalEnergy(VecDyn, Bases, Ham0 );
  file2->SaveStdVector(gname, "Ei", Ei);
  delete file2;
  // HDF5IO* file3 = new HDF5IO(prefix + SaveFile + "-WF.h5");
  // gname = "WF";
  // file3->SaveVector(gname, "Vec", VecDyn);
  // delete file3;
  LogOut << "Begin dynamics......" << std::endl;
  for (size_t cntP = 1; cntP <= TSteps; cntP++) {
    // Evolve the state
    Ham0.expH(Prefactor, VecDyn);
    VecDyn = arma::normalise(VecDyn);
    Ham0.expH(Prefactor, JVec);
    JVec = arma::normalise(JVec);
    if ( cntP % MeasureEvery == 0 ){
      file2 = new HDF5IO(prefix + SaveFile + ".h5");
      std::string gname = "Obs-";
      gname.append( std::to_string((unsigned long long)cntP ));
      gname.append("/");
      Lecho = arma::cdot(VecInit, VecDyn);
      file2->SaveNumber(gname, "F", Lecho);
      JVecDyn = JPsi(Bases, VecDyn, Ham0);
      JJt = arma::cdot(JVecDyn, JVec);
      file2->SaveNumber(gname, "JJt", JJt);
      ComplexVectorType Npi = NPhonon( Bases, VecDyn, Ham0);
      file2->SaveVector(gname, "Phonon", Npi);
      ComplexMatrixType Nfi = NFermion( Bases, VecDyn, Ham0);
      file2->SaveMatrix(gname, "Fermion", Nfi);
      Ei = LocalEnergy(VecDyn, Bases, Ham0 );
      file2->SaveStdVector(gname, "Ei", Ei);
      delete file2;
    }
    // if ( cntP % SaveWFEvery == 0 ){
    //   file3 = new HDF5IO(prefix + SaveFile + "-WF.h5");
    //   gname = "WF-";
    //   gname.append( std::to_string((unsigned long long)cntP ));
    //   gname.append("/");
    //   file3->SaveVector(gname, "Vec", VecDyn);
    //   delete file3;
    // }
  }
  // file3 = new HDF5IO(prefix + SaveFile + "-WF.h5");
  // gname = "WF";
  // file3->SaveVector(gname, "Vec", VecDyn);
  // delete file3;
  LogOut << "Finished dynamics!!" << std::endl;

  LogOut.close();
}

/* main program */
int main(int argc, char *argv[]){
  if ( argc < 2 ) RUNTIME_ERROR(" Need at least one argument to run program. Use 0 to run Equilibrium.");
#ifdef MKL
  mkl_set_num_threads(NumCores);
#endif
  int world_size;
  int world_rank;
  std::vector<std::string> MPIFolders(1, "");
#ifdef MPIPARALLEL
  std::vector<std::string> MPIFolders = GetPrefix("MPIFolders");
  // Initialize MPI
  MPI_Init(NULL, NULL);
  // Get the number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  // Get the rank of the process
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  assert ( MPIFolders.size() == world_size );
  if ( std::atoi(argv[1]) == 0 ){
    int NEV = 40;
    if ( argc > 2 ) NEV = std::atoi(argv[2]);
    Equilibrium("", NEV);
  }else if ( std::atoi(argv[1]) == 1 ){
    std::string InitialState = "R";
    int S1 = -1, S2 = -1;
    int SaveWFEvery = 100, MeasureEvery = 20;
    if ( argc > 2 ) InitialState = argv[2];
    if ( argc > 4 ){
      S1 = std::atoi(argv[3]);
      S2 = std::atoi(argv[4]);
    }
    if ( argc > 5 ) SaveWFEvery = std::atoi(argv[5]);
    if ( argc > 6 ) MeasureEvery = std::atoi(argv[6]);
    Dynamics(MPIFolders.at(world_rank), InitialState, S1, S2, MeasureEvery, SaveWFEvery);
  }
  MPI_Finalize();
#else
  world_size = MPIFolders.size();
  world_rank = 0;
  if ( std::atoi(argv[1]) == 0 ){
    int NEV = 40;
    if ( argc > 2 ) NEV = std::atoi(argv[2]);
    Equilibrium("", NEV);
  }else if ( std::atoi(argv[1]) == 1 ){
    std::string InitialState = "R";
    int S1 = 0, S2 = -1;
    int SaveWFEvery = 100, MeasureEvery = 20;
    if ( argc > 2 ) InitialState = argv[2];
    if ( argc > 4 ){
      S1 = std::atoi(argv[3]);
      S2 = std::atoi(argv[4]);
    }
    if ( argc > 5 ) SaveWFEvery = std::atoi(argv[5]);
    if ( argc > 6 ) MeasureEvery = std::atoi(argv[6]);
    Dynamics(MPIFolders.at(world_rank), InitialState, S1, S2, MeasureEvery, SaveWFEvery);
  }
#endif
  return 0;
}
