#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include "src/EDType.hpp"
#include "src/Node/Node.hpp"
#include "src/Lattice/preset.hpp"
#include "src/Basis/Basis.hpp"
#include "src/Hamiltonian/Holstein/Holstein.hpp"
#include "src/hdf5io/hdf5io.hpp"

#ifdef MKL
  #include "mkl.h"
#endif

#if not defined(NumCores)
  #define NumCores 16
#endif

//? cSpell:words eigenstate phonon Diagonalize

const int WithoutK0Phonon = 1;
int L = 4;
int N = 5 * L;
const RealType Jin = 1.0;
RealType Win = 0.50;
RealType Gin = 1.00;

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

RealVectorType NPhonon( const std::vector<Basis> &Bases, const RealVectorType &Vec){
  int L = Bases.at(0).GetL();
  RealVectorType out(L-1, arma::fill::zeros);
  std::vector< std::vector<int> > b = Bases.at(0).GetBStates();
  assert( b.size() == Vec.size() );
  int cnt = 0;
  for ( auto &nbi : b ){
    for (size_t cnt2 = 0; cnt2 < L-1; cnt2++) {//* Phonon k
      out.at(cnt2) += (RealType)nbi.at(cnt2+WithoutK0Phonon) * std::pow(std::abs(Vec(cnt)), 2);
    }
    cnt++;
  }
  return out;
}

RealVectorType NFermion( const std::vector<Basis> &Bases, const RealVectorType &Vec){
  int L = Bases.at(0).GetL();
  RealVectorType out(L, arma::fill::zeros);
  std::vector< std::vector<int> > b = Bases.at(0).GetBStates();
  assert( b.size() == Vec.size() );
  int cnt = 0;
  for ( auto &nbi : b ){
    out.at(nbi.at(0)) += std::pow(std::abs(Vec(cnt)), 2);
    cnt++;
  }
  return out;
}

std::vector<int> CkCq( const std::vector<int> &Kn, const std::map<int, std::vector<Basis> > BasisMap, const int ki, const int qi, const int Kti, int &Kfi){//! <c^\dagger_k c_q>
  std::vector<int> out;
  assert( ki != qi );
  Basis Bk = BasisMap.at(Kti).at(0);
  int Kf = Kn.at(Kti) + Kn.at(ki) - Kn.at(qi);
  Kfi = DeltaKIndex( Kf, Kn );
  Basis Bq = BasisMap.at(Kfi).at(0);
  std::vector< std::vector<int> > b = Bk.GetBStates();
  int cnt = 0;
  for ( auto &nbi : b ){
    if ( nbi.at(0) == qi ){
      std::vector<int> wb = nbi;
      wb.at(0) = ki;//* c^\dagger_k c_q |KetState>
      double wb_tag = BosonBasisTag(wb);
      size_t index = Bq.GetIndexFromTag( wb_tag );//* This is the index for BraState
      if ( !(Bq.DummyCheckState(wb_tag, wb)) ) throw std::runtime_error("States mismatch in CkCq");
      out.push_back(index);
    }else{
      out.push_back(-1);
    }
    cnt++;
  }
  return out;
}

RealVectorType Jq0State(const std::vector<Basis> &Bases, const std::vector<int> &Kn, const RealVectorType &Vec){
  int L = Bases.at(0).GetL();
  std::vector< std::vector<int> > b = Bases.at(0).GetBStates();
  RealVectorType out = Vec;
  assert( b.size() == Vec.size() );
  assert( b.size() == out.size() );
  int cnt = 0;
  for ( auto &nbi : b ){
    int Kfi = Kn.at(nbi.at(0));
    RealType Kf = Kfi * PI / RealType(L/2);
    out(cnt) = sin(Kf) * Vec(cnt);//! No normalizaztion
    cnt++;
  }
  return out;
}

RealMatrixType Jq0Vecs(const std::vector<Basis> &Bases, const std::vector<int> &Kn, const RealMatrixType &Vecs){
  int L = Bases.at(0).GetL();
  std::vector< std::vector<int> > b = Bases.at(0).GetBStates();
  RealMatrixType out(Vecs.n_rows, Vecs.n_cols, arma::fill::zeros);
  for ( size_t i = 0; i < Vecs.n_cols; i++ ){
    RealVectorType vec = Vecs.col(i);
    RealVectorType J0vec = Jq0State(Bases, Kn, vec);
    out.col(i) = J0vec;
  }
  return out;
}

void SpectralPeaks( const double Eg, const RealVectorType AS, const RealVectorType &EigVal, const RealMatrixType &EigVec,
  std::vector<double> &PeakLocation, std::vector<double> &PeakWeight, const int MaxNumPeak){
  PeakLocation.clear();
  PeakWeight.clear();
  for ( size_t i = 0; i < EigVec.n_cols; i++){
    PeakLocation.push_back(EigVal(i) - Eg);
    RealType val = arma::cdot(EigVec.col(i), AS);
    PeakWeight.push_back( RealPart(val * Conjg(val)) );
    if ( i > MaxNumPeak ) break;
  }
  assert( PeakLocation.size() == PeakWeight.size() );
}

void Equilibrium(const std::string prefix, int NEV, int NumPeaks){
  std::ofstream LogOut;
  LogOut.open(prefix + "Holstein.K.eqm", std::ios::app);
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
  //* k points */
  LogOut << "Build k-point - " << std::flush;
  std::vector<int> Kn;
  Kn.clear();
  Kn.push_back(0);
  for ( size_t n = 1; n <= L/2 - 1; n++ ){
    Kn.push_back(+n);
    Kn.push_back(-n);
  }
  Kn.push_back(L/2);
  // //PrintVector(Kn, 3, " ");
  assert( Kn.size() == L );
  LogOut << Kn.size() << " k-points DONE!" << std::endl;

  std::map<int, std::vector<Basis> > BasisSets;
  for ( size_t k = 0; k < Kn.size(); k++ ){
    bool FullSpectrum = false;
    int TargetK = Kn.at(k);
    LogOut << "Build Basis - max. phonon quanta = " << N << " with target total momentum K = " << double(TargetK) << ", ends up with Hilbert space = " << std::flush;
    Basis P1(L, N);
    P1.PhononK(Kn, TargetK, WithoutK0Phonon);
    LogOut << P1.GetHilbertSpace() << std::flush;
    // //std::cout << TargetK << std::endl;
    // //std::cout << P1 << std::endl;
    std::vector<Basis> Bases;
    Bases.push_back(P1);
    BasisSets[k] = Bases;
    LogOut << " DONE!" << std::endl;

    LogOut << "Build Hamiltonian - " << std::flush;
    Holstein<RealType> Ham0( Bases );
    Ham0.HolsteinModelK(Kn, Bases, Win, Gin, Jin, N, WithoutK0Phonon);
    LogOut << "Hermitian = " << Ham0.CheckHermitian() << ", Hilbert space = " << Ham0.GetTotalHilbertSpace() << ", DONE!" << std::endl;
    LogOut << "Diagonalize Hamiltonian to find - " << Target << " - " << std::flush;
    RealVectorType Vals;
    RealMatrixType Vecs;
    if ( Ham0.GetTotalHilbertSpace() > 2100 ){
      Ham0.eigh(Vals, Vecs, NEV, true, Target);//* ARPACK
    }else{
      Ham0.diag(Vals, Vecs);//* Full spectrum
      NEV = Ham0.GetTotalHilbertSpace();
      FullSpectrum = true;
    }
    LogOut << "DONE!" << std::endl;
    LogOut << "\tEnergy 0 : " << std::setprecision(12) << Vals[0] << std::endl;
    LogOut << "\tEnergy 1 : " << std::setprecision(12) << Vals[1] << std::endl;
    LogOut << "\tEnergy 2 : " << std::setprecision(12) << Vals[2] << std::endl;
    LogOut << "\tEnergy 3 : " << std::setprecision(12) << Vals[3] << std::endl;
    LogOut << "\tEnergy 4 : " << std::setprecision(12) << Vals[4] << std::endl;
    LogOut << "\tEnergy 5 : " << std::setprecision(12) << Vals[5] << std::endl;
    HDF5IO* file = new HDF5IO(prefix + "Holstein.K.h5");
    std::string gnameK = std::to_string( (unsigned long) k );
    file->SaveNumber(gnameK, "K", TargetK);
    RealMatrixType JVecs = Jq0Vecs(Bases, Kn, Vecs);
    RealMatrixType Jmn = Vecs.t() * JVecs;//* <m| j |n>
    file->SaveMatrix(gnameK, "Jmn", Jmn);
    for ( size_t i = 0; i < NEV; i++ ){
      std::string gname = gnameK;
      gname.append("/S-");
      gname.append( std::to_string( (unsigned long)i) );
      file->SaveNumber(gname, "EVal", Vals[i]);
      RealVectorType Vec = Vecs.col(i);
      file->SaveVector(gname, "EVec", Vec);
      RealVectorType Npi = NPhonon( Bases, Vec);
      file->SaveVector(gname, "Phonon", Npi);
      RealVectorType Nfi = NFermion( Bases, Vec);
      file->SaveVector(gname, "Fermion", Nfi);
      LogOut << std::setw(3) << i << ": E = " << std::setprecision(12) << std::setw(15) << Vals[i] << ", Np = " << std::setw(14) << RealPart(arma::accu(Npi)) << ", Nf = [" << std::flush;
      for ( size_t j = 0; j < L; j++) LogOut << std::setw(16) << RealPart(Nfi.at(j)) << ", " << std::flush;
      LogOut << "], Total = " << RealPart(arma::accu(Nfi)) << std::endl;
      if ( NumPeaks ){
        std::vector<double> PeakLocation;
        PeakLocation.clear();
        std::vector<double> PeakWeight;
        PeakWeight.clear();
        RealVectorType JVec = Jq0State(Bases, Kn, Vec);
        RealVectorType wVals(NumPeaks, arma::fill::zeros);
        RealMatrixType wVecs(Ham0.GetTotalHilbertSpace(), NumPeaks, arma::fill::zeros);
        if ( !(FullSpectrum) ){
          wVecs.col(0) = JVec;//* It got normalized inside SpectralH.
          Ham0.SpectralH( wVals, wVecs, JVec, NumPeaks );
        }else{
          NumPeaks = Vecs.n_cols;
          wVals = Vals;
          wVecs = Vecs;
          // //wVecs.col(0) = JVec;
          // //Ham0.SpectralH( wVals, wVecs, JVec, NumPeaks );
        }
        // JVec = arma::normalise(JVec);//! Check this to be 1
        SpectralPeaks(Vals[i], JVec, wVals, wVecs, PeakLocation, PeakWeight, NumPeaks);
        file->SaveStdVector(gname, "JJw", PeakLocation);
        file->SaveStdVector(gname, "JJh", PeakWeight);
        double Wt = std::accumulate(PeakWeight.begin(), PeakWeight.end(), 0.0);
        file->SaveNumber(gname, "Wt", Wt);
      }
    }
    delete file;
  }
  for ( int ki = 0; ki < L; ki++){
    for ( int qi = 0; qi < L; qi++){
      if ( ki == qi ) continue;//* No momentum change - diagonal!
      for ( int PsiK = 0; PsiK < L; PsiK++){
        int PsiB;
        std::vector<int> CkCqIndex = CkCq( Kn, BasisSets, ki, qi, PsiK, PsiB);
        HDF5IO* file = new HDF5IO(prefix + "Holstein.K.h5");
        std::string gname = "S";
        gname.append( std::to_string( (unsigned long)PsiB) );
        gname.append( "-Cd" );
        gname.append( std::to_string( (unsigned long)ki) );
        gname.append( "C" );
        gname.append( std::to_string( (unsigned long)qi) );
        gname.append( "-S" );
        gname.append( std::to_string( (unsigned long)PsiK) );
        file->SaveStdVector("CkCqMap", gname, CkCqIndex);
        delete file;
      }
    }
  }
  LogOut.close();
}

//* main program
int main(int argc, char *argv[]){
  if ( argc < 2 ) RUNTIME_ERROR(" Need at least one argument to run program. Use 0 to run Equilibrium.");
  #ifdef MKL
    mkl_set_num_threads(NumCores);
  #endif
  if ( std::atoi(argv[1]) == 0 ){
    int NEV = 40;
    if ( argc > 2 ) NEV = std::atoi(argv[2]);
    int NumPeaks = 0;
    if ( argc > 3 ) NumPeaks = std::atoi(argv[3]);
    Equilibrium("", NEV, NumPeaks);
  }
  return 0;
}
