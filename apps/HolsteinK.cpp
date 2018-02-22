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

#ifdef MKL
  #include "mkl.h"
#endif

void LoadEqmParameters( const std::string filename, int& L, int& N, RealType& G, RealType& W, RealType& EShift, std::string& Method){
  HDF5IO h5f(filename);
  h5f.LoadNumber("Parameters", "L", L);
  h5f.LoadNumber("Parameters", "N", N);
  h5f.LoadNumber("Parameters", "G", G);
  h5f.LoadNumber("Parameters", "W", W);
  h5f.LoadNumber("Parameters", "EShift", EShift);
  int Mid = 0;
  h5f.LoadNumber("Parameters", "Method", Mid);
  if ( Mid == 0 ) Method = "SR";
  else if ( Mid == 2 ) Method = "SM";
  else if ( Mid == 1 ) Method = "LR";
}

ComplexVectorType NPhonon( const std::vector<Basis> &Bases, const ComplexVectorType &Vec, const Holstein<ComplexType>& Ham){
  int L = Bases.at(0).GetL() + 1;
  // int N = Bases.at(0).GetN();
  ComplexVectorType out(L-1, arma::fill::zeros);
  std::vector< std::vector<int> > b = Bases.at(0).GetBStates();
  assert( b.size() * L == Vec.size() );
  for ( size_t cnt = 0; cnt < L; cnt++ ){
    int coff = 0;
    for ( auto &nbi : b ){
      for (size_t cnt = 0; cnt < L-1; cnt++) {
        size_t idx = Ham.DetermineTotalIndex( vec<size_t>(cnt, coff) );
        out.at(cnt) += (RealType)nbi.at(cnt) * std::pow(std::abs(Vec(idx)), 2);
      }
      coff++;
    }
  }
  return out;
}

ComplexMatrixType NFermion( const std::vector<Basis> &Bases, const ComplexVectorType &Vec, const Holstein<ComplexType>& Ham){
  int L = Bases.at(0).GetL() + 1;
  // int N = Bases.at(0).GetN();
  ComplexMatrixType out(L, L, arma::fill::zeros);
  assert( Bases.at(0).GetHilbertSpace() * L == Vec.size() );
  for ( size_t s1 = 0; s1 < L; s1++ ){
    for ( size_t s2 = 0; s2 < L; s2++ ){
      int coff = 0;
      for ( size_t b = 0; b < Bases.at(0).GetHilbertSpace(); b++ ){
        size_t idx1 = Ham.DetermineTotalIndex( vec<size_t>(s1, b) );
        size_t idx2 = Ham.DetermineTotalIndex( vec<size_t>(s2, b) );
        out.at(s1, s2) += Vec(idx2) * Conjg(Vec(idx1));
      }
    }
  }
  return out;
}

void Equilibrium(const std::string prefix, int NEV){
  std::ofstream LogOut;
  LogOut.open(prefix + "Holstein.K.eqm", std::ios::app);
  int L = 4;
  int N = 10 * L;
  const RealType Jin = 1.0;
  RealType Win = 10.0;
  RealType Gin = 10.0;
  RealType EShift = 0;
  std::string Target;
  if ( std::abs(EShift) > 1.0e-5 ) Target = "SM";
  else Target = "SR";

  try{
    H5::Exception::dontPrint();
    H5::H5File::isHdf5(prefix + "conf.h5");
    LoadEqmParameters( prefix + "conf.h5", L, N, Gin, Win, EShift, Target);
  }catch(H5::FileIException){
    LogOut << "Use default settings." << std::endl;
  }
  /* k points */
  LogOut << "Build k-point - " << std::flush;
  std::vector<int> Kn;
  Kn.clear();
  Kn.push_back(0);
  for ( size_t n = 1; n <= L/2 - 1; n++ ){
    Kn.push_back(+n);
    Kn.push_back(-n);
  }
  Kn.push_back(L/2);
  // PrintVector(Kn, 3, " ");
  assert( Kn.size() == L );
  LogOut << Kn.size() << " k-points DONE!" << std::endl;

  LogOut << "Build Basis - max. phonon quanta = " << N << ", ends up with Phonon Hilbert space = " << std::flush;
  Basis P1(L-1, N);// Get rid og k=0 phonon mode, so using L - 1!!
  P1.PhononK();
  LogOut << P1.GetHilbertSpace() << std::flush;
  // LogOut << P1 << std::endl;
  std::vector<Basis> Bases;
  Bases.push_back(P1);
  LogOut << " DONE!" << std::endl;

  LogOut << "Build Hamiltonian - " << std::flush;
  Holstein<ComplexType> Ham0( Kn.size(), Bases );
  Ham0.HolsteinModelK(Kn, Bases, Win, Gin, Jin, N);
  LogOut << "Hermitian = " << Ham0.CheckHermitian() << ", Hilbert space = " << Ham0.GetTotalHilbertSpace() << ", DONE!" << std::endl;

  LogOut << "Shift energy by - " << EShift << std::flush;
  Ham0.ShiftEnergy(EShift);
  LogOut << ", Hermitian = " << Ham0.CheckHermitian() << ", Hilbert space = " << Ham0.GetTotalHilbertSpace() << ", DONE!" << std::endl;
  LogOut << "Diagonalize Hamiltonian to find - " << Target << " - " << std::flush;
  RealVectorType Vals;
  ComplexMatrixType Vecs;
  if ( Ham0.GetTotalHilbertSpace() > 1200 ){
    Ham0.eigh(Vals, Vecs, NEV, true, Target);
  }else{
    Ham0.diag(Vals, Vecs);// Full spectrum
    NEV = Ham0.GetTotalHilbertSpace();
  }
  LogOut << "DONE!" << std::endl;
  LogOut << "\tGS energy = " << std::setprecision(12) << Vals[0] << std::endl;
  LogOut << "\tFES energy = " << std::setprecision(12) << Vals[1] << std::endl;
  HDF5IO* file = new HDF5IO(prefix + "HolsteinK.h5");
  for ( size_t i = 0; i < NEV; i++ ){
    std::string gname = "S-";
    file->SaveNumber(gname, "EVal", Vals[i]);
    gname.append( std::to_string( (unsigned long)i) );
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

void Dynamics(const std::string prefix, const int S1, const int S2, const int MeasureEvery = 20, const int SaveWFEvery = 1000000, const std::string InitialState = "R" ){
  std::ofstream LogOut;
  LogOut.open(prefix + "Holstein.K.dyn", std::ios::app);
  int L = 2;
  int N = 10 * L;
  const RealType Jin = 1.0;
  RealType Win = 10.0;
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
  /* k points */
  LogOut << "Build k-point - " << std::flush;
  std::vector<int> Kn;
  Kn.clear();
  Kn.push_back(0);
  for ( size_t n = 1; n <= L/2 - 1; n++ ){
    Kn.push_back(+n);
    Kn.push_back(-n);
  }
  Kn.push_back(L/2);
  // PrintVector(Kn, 3, " ");
  assert( Kn.size() == L );
  LogOut << Kn.size() << " k-points DONE!" << std::endl;

  LogOut << "Build Basis - max. phonon quanta = " << N << ", ends up with Phonon Hilbert space = " << std::flush;
  Basis P1(L-1, N);// Get rid og k=0 phonon mode, so using L - 1!!
  P1.PhononK();
  LogOut << P1.GetHilbertSpace() << std::flush;
  std::vector<Basis> Bases;
  Bases.push_back(P1);
  LogOut << " DONE!" << std::endl;

  LogOut << "Build Hamiltonian - " << std::flush;
  Holstein<ComplexType> Ham0( Kn.size(), Bases );
  Ham0.HolsteinModelK(Kn, Bases, Win, Gin, Jin, N);
  LogOut << "Hermitian = " << Ham0.CheckHermitian() << ", Hilbert space = " << Ham0.GetTotalHilbertSpace() << ", DONE!" << std::endl;

  LogOut << "Load Wavefunction - " << std::flush;
  ComplexVectorType VecInit(Ham0.GetTotalHilbertSpace(), arma::fill::randn);
  try{
    H5::Exception::dontPrint();
    H5::H5File::isHdf5(prefix + "HolsteinK.h5");
    HDF5IO* file = new HDF5IO(prefix + "HolsteinK.h5");
    ComplexVectorType V1, V2;
    std::string gname = "S-";
    gname.append( std::to_string((unsigned long)S1) );
    LogOut << " <" << gname << std::flush;
    file->LoadVector(gname, "EVec", V1);
    gname = "S-";
    gname.append( std::to_string((unsigned long)S2) );
    LogOut << "|" << gname << std::flush;
    file->LoadVector(gname, "EVec", V2);
    delete file;
    VecInit = ( V1 + V2 );
    LogOut << "> = " << arma::cdot(V1, V2) << " DONE." << std::endl;
  }catch(H5::FileIException){
    if ( InitialState == "Zero" ){
      for ( size_t f = 0; f < L; f++){
        for ( size_t b = 0; Bases.at(0).GetHilbertSpace(); b++ ){
          size_t idx = Ham0.DetermineTotalIndex( vec<size_t>(f, b) );
          if ( b == 0 ) VecInit[idx] = ComplexType(1.0e0, 0.0e0);
          else VecInit[idx] = ComplexType(0.0e0, 0.0e0);
        }
      }
    }else{
      LogOut << "Use random initial." << std::endl;
    }
  }
  VecInit = arma::normalise(VecInit);

  ComplexType Prefactor = ComplexType(0.0, -1.0e0*dt);/* NOTE: hbar = 1 */
  ComplexVectorType VecDyn = VecInit;
  HDF5IO* file2 = new HDF5IO("QuenchState.h5");
  std::string gname = "Obs-0/";
  file2->SaveNumber(gname, "S1", S1);
  file2->SaveNumber(gname, "S2", S2);
  ComplexType Lecho = arma::cdot(VecInit, VecDyn);
  file2->SaveNumber(gname, "F", Lecho);
  ComplexVectorType Npi = NPhonon( Bases, VecDyn, Ham0);
  file2->SaveVector(gname, "Phonon", Npi);
  ComplexMatrixType Nfi = NFermion( Bases, VecDyn, Ham0);
  file2->SaveMatrix(gname, "Fermion", Nfi);
  delete file2;
  LogOut << "Begin dynamics......" << std::endl;
  for (size_t cntP = 1; cntP <= TSteps; cntP++) {
    // Evolve the state
    Ham0.expH(Prefactor, VecDyn);
    VecDyn = arma::normalise(VecDyn);
    if ( cntP % MeasureEvery == 0 ){
      file2 = new HDF5IO("QuenchState.h5");
      std::string gname = "Obs-";
      gname.append( std::to_string((unsigned long long)cntP ));
      gname.append("/");
      ComplexType Lecho = arma::cdot(VecInit, VecDyn);
      file2->SaveNumber(gname, "F", Lecho);
      ComplexVectorType Npi = NPhonon( Bases, VecDyn, Ham0);
      file2->SaveVector(gname, "Phonon", Npi);
      ComplexMatrixType Nfi = NFermion( Bases, VecDyn, Ham0);
      file2->SaveMatrix(gname, "Fermion", Nfi);
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
  HDF5IO* file3 = new HDF5IO("QuenchStateWF.h5");
  gname = "WF";
  file3->SaveVector(gname, "Vec", VecDyn);
  delete file3;
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
    int SaveWFEvery = 1000000, MeasureEvery = 20;
    int S1 = 8, S2 = 9;
    if ( argc > 3 ){
      S1 = std::atoi(argv[2]);
      S2 = std::atoi(argv[3]);
    }
    if ( argc > 4 ) SaveWFEvery = std::atoi(argv[4]);
    if ( argc > 5 ) MeasureEvery = std::atoi(argv[5]);
    Dynamics("", S1, S2, MeasureEvery, SaveWFEvery);
  }
  return 0;
}
