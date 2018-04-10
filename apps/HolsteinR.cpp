#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
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
  for ( size_t cnt = 0; cnt < L; cnt++ ){// phono r
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
  Basis P1(L, N);// Get rid og k=0 phonon mode, so using L - 1!!
  P1.Phonon();
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
  ComplexVectorType Vec0(Ham0.GetTotalHilbertSpace(), arma::fill::zeros);
  ComplexVectorType Vec1(Ham0.GetTotalHilbertSpace(), arma::fill::zeros);
  ComplexVectorType Vec2(Ham0.GetTotalHilbertSpace(), arma::fill::zeros);
  std::vector< std::vector<int> > b = Bases.at(0).GetBStates();
  for ( size_t i = 0; i < b.size(); i++ ){
    size_t FermionLocation = 0;
    size_t idx = Ham0.DetermineTotalIndex( vec<size_t>(FermionLocation, i) );
    Vec0(idx) = VecIn(idx);
    FermionLocation = 1;
    idx = Ham0.DetermineTotalIndex( vec<size_t>(FermionLocation, i) );
    Vec1(idx) = VecIn(idx);
    FermionLocation = 2;
    idx = Ham0.DetermineTotalIndex( vec<size_t>(FermionLocation, i) );
    Vec2(idx) = VecIn(idx);
  }
  Vec0 = arma::normalise(Vec0);
  Vec1 = arma::normalise(Vec1);
  Vec2 = arma::normalise(Vec2);
  RealType E0 = RealPart( arma::cdot(Vec0, Ham0.GetTotalHamiltonian() * Vec0) );
  RealType E1 = RealPart( arma::cdot(Vec1, Ham0.GetTotalHamiltonian() * Vec1) );
  RealType E2 = RealPart( arma::cdot(Vec2, Ham0.GetTotalHamiltonian() * Vec2) );
  return vec(E0, E1, E2);
}

ComplexVectorType CoherentState( std::ofstream& LogOut, const std::vector<ComplexType>& alphas, const std::vector<Basis>& Bases, const Holstein<ComplexType>& Ham0){
  ComplexVectorType Vec0(Ham0.GetTotalHilbertSpace(), arma::fill::zeros);
  ComplexVectorType Vec1(Ham0.GetTotalHilbertSpace(), arma::fill::zeros);
  ComplexVectorType Vec2(Ham0.GetTotalHilbertSpace(), arma::fill::zeros);
  std::vector< std::vector<int> > b = Bases.at(0).GetBStates();
  for ( size_t i = 0; i < b.size(); i++ ){
    std::vector<int> nb = b.at(i);
    ComplexType Coff(0.0e0, 0.0e0);
    for ( size_t j = 0; j < nb.size(); j++ ){
      int Nbj = nb.at(j);
      ComplexType CoffTmp(1.0e0, 0.0e0);
      for ( size_t k = 1; k <= Nbj; k++ ){
        // Check this
        ComplexType tmp = alphas.at(j) / std::sqrt(RealType(k));
        CoffTmp *= tmp;
        // if ( std::abs(tmp) < 1.0e-10 ) break;
      }
      if ( j == 0 ) Coff = CoffTmp;
      else Coff *= CoffTmp;
    }
    size_t FermionLocation = 0;
    size_t idx = Ham0.DetermineTotalIndex( vec<size_t>(FermionLocation, i) );
    Vec0(idx) = 1.0e-1 * Coff;
    FermionLocation = 1;
    idx = Ham0.DetermineTotalIndex( vec<size_t>(FermionLocation, i) );
    Vec1(idx) = 1.0e-1 * Coff;
    FermionLocation = 2;
    idx = Ham0.DetermineTotalIndex( vec<size_t>(FermionLocation, i) );
    Vec2(idx) = 1.0e-1 * Coff;
  }
  Vec0 = arma::normalise(Vec0);
  Vec1 = arma::normalise(Vec1);
  Vec2 = arma::normalise(Vec2);
  RealType E0 = RealPart( arma::cdot(Vec0, Ham0.GetTotalHamiltonian() * Vec0) );
  RealType E1 = RealPart( arma::cdot(Vec1, Ham0.GetTotalHamiltonian() * Vec1) );
  RealType E2 = RealPart( arma::cdot(Vec2, Ham0.GetTotalHamiltonian() * Vec2) );
  /* At zero */
  LogOut << E0 << " " << E1 << " " << E2 << std::endl;
  if ( E0 < E1 ){
    if ( E0 < E2 ) return Vec0;
    else return Vec2;
  }else{
    if ( E1 < E2 ) return Vec1;
    else return Vec2;
  }
  /* Eigen state */
  // RealMatrixType Mat(3,3,arma::fill::zeros);
  // Mat(0,0) = E0;
  // Mat(1,1) = E1;
  // Mat(2,2) = E2;
  // Mat.diag(1) -= 1.0e0;
  // Mat.diag(-1) -= 1.0e0;
  // Mat.diag(2) -= 1.0e0;
  // Mat.diag(-2) -= 1.0e0;
  // RealType* EigVec = (RealType*)malloc( 3 * 3 * sizeof(RealType) );
  // RealType* Eig = (RealType*)malloc( 3 * sizeof(RealType) );
  // syDiag(Mat.memptr(), 3, Eig, EigVec);
  // LogOut << "Eigenvalues - " << std::endl;
  // LogOut << Eig[0] << " " << Eig[1] << " " << Eig[2] << std::endl;
  // LogOut << EigVec[0] << " " << EigVec[1] << " " << EigVec[2] << std::endl;
  // return EigVec[0] * Vec0 + EigVec[1] * Vec1 + EigVec[2] * Vec2;
}

void Dynamics(const std::string prefix, const std::string InitialState, const int S1, const int S2, const int MeasureEvery, const int SaveWFEvery ){
  std::ofstream LogOut;
  LogOut.open(prefix + "Holstein.R.dyn", std::ios::app);
  int L = 3;
  int N = 15 * L;
  const RealType Jin = 1.0;
  RealType Win = 0.30;
  RealType Gin = 10.0;
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
  P1.Phonon();
  LogOut << P1.GetHilbertSpace() << std::flush;
  std::vector<Basis> Bases;
  Bases.push_back(P1);
  LogOut << " DONE!" << std::endl;

  LogOut << "Build Hamiltonian - " << std::flush;
  Holstein<ComplexType> Ham0( L, Bases );
  Ham0.HolsteinModelR(L, Bases, Win, Gin, Jin);
  LogOut << "Hermitian = " << Ham0.CheckHermitian() << ", Hilbert space = " << Ham0.GetTotalHilbertSpace() << ", DONE!" << std::endl;

  LogOut << "Load Wavefunction - " << InitialState << " " << std::flush;
  // ComplexVectorType VecInit(Ham0.GetTotalHilbertSpace(), arma::fill::zeros);
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
    SaveFile.append( "-Z" );
    /* Target phonon mode - coherent state */
    std::vector<ComplexType> alphas;
    try{
      H5::Exception::dontPrint();
      H5::H5File::isHdf5(prefix + "conf.h5");
      LoadAlphas( prefix + "conf.h5", alphas );
    }catch(H5::FileIException){
      alphas.push_back( 5.2 * exp(0.3*PI*ComplexType(0.0,1.0)) );
      alphas.push_back( 7.9 * exp(0.7*PI*ComplexType(0.0,1.0)) );
      alphas.push_back(10.1 * exp(-0.2*PI*ComplexType(0.0,1.0)) );
      alphas.push_back( 3.2 * exp(0.4*PI*ComplexType(0.0,1.0)) );
      LogOut << "Use default settings." << std::endl;
    }
    LogOut << "Fermion is at GS, and phonon coherent state: alpha_i = " << std::flush;
    for ( auto val : alphas ){
      LogOut << val << " " << std::flush;
    }
    VecInit = CoherentState(LogOut, alphas, Bases, Ham0);
    LogOut << " DONE!" << std::endl;
  }else{
    SaveFile.append( "-R" );
    // VecInit = ComplexVectorType(Ham0.GetTotalHilbertSpace(), arma::fill::randn);
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
  ComplexVectorType Npi = NPhonon( Bases, VecDyn, Ham0);
  file2->SaveVector(gname, "Phonon", Npi);
  ComplexMatrixType Nfi = NFermion( Bases, VecDyn, Ham0);
  file2->SaveMatrix(gname, "Fermion", Nfi);
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
    if ( cntP % MeasureEvery == 0 ){
      file2 = new HDF5IO(prefix + SaveFile + ".h5");
      std::string gname = "Obs-";
      gname.append( std::to_string((unsigned long long)cntP ));
      gname.append("/");
      ComplexType Lecho = arma::cdot(VecInit, VecDyn);
      file2->SaveNumber(gname, "F", Lecho);
      ComplexVectorType Npi = NPhonon( Bases, VecDyn, Ham0);
      file2->SaveVector(gname, "Phonon", Npi);
      ComplexMatrixType Nfi = NFermion( Bases, VecDyn, Ham0);
      file2->SaveMatrix(gname, "Fermion", Nfi);
      std::vector<RealType> Ei = LocalEnergy(VecDyn, Bases, Ham0 );
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
    Dynamics("", InitialState, S1, S2, MeasureEvery, SaveWFEvery);
  }
  return 0;
}
