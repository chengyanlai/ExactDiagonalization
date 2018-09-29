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

const int WithoutK0Phonon = 1;
int L = 2;
int N = 10 * L;
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

ComplexVectorType NPhonon( const std::vector<Basis> &Bases, const ComplexVectorType &Vec, const Holstein<ComplexType>& Ham){
  int L = Bases.at(0).GetL();
  ComplexVectorType out(L-1, arma::fill::zeros);
  std::vector< std::vector<int> > b = Bases.at(0).GetBStates();
  assert( b.size() == Vec.size() );
  int coff = 0;
  for ( auto &nbi : b ){
    for (size_t cnt2 = 0; cnt2 < L-1; cnt2++) {//* Phonon k
      out.at(cnt2) += (RealType)nbi.at(cnt2+WithoutK0Phonon) * std::pow(std::abs(Vec(coff)), 2);
    }
    coff++;
  }
  return out;
}

ComplexVectorType NFermion( const std::vector<Basis> &Bases, const ComplexVectorType &Vec, const Holstein<ComplexType>& Ham){
  int L = Bases.at(0).GetL();
  ComplexVectorType out(L, arma::fill::zeros);
  std::vector< std::vector<int> > b = Bases.at(0).GetBStates();
  assert( b.size() == Vec.size() );
  int coff = 0;
  for ( auto &nbi : b ){
    out.at(nbi.at(0)) += std::pow(std::abs(Vec(coff)), 2);
    coff++;
  }
  return out;
}

void Equilibrium(const std::string prefix, int NEV){
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

  for ( size_t k = 0; k < Kn.size(); k++ ){
    int TargetK = Kn.at(k);
    LogOut << "Build Basis - max. phonon quanta = " << N << " with target total momentum K = " << double(TargetK) << ", ends up with Hilbert space = " << std::flush;
    Basis P1(L, N);//! Get rid og k=0 phonon mode, so using L - 1!!
    P1.PhononK(Kn, TargetK, WithoutK0Phonon);
    LogOut << P1.GetHilbertSpace() << std::flush;
    // //std::cout << TargetK << std::endl;
    // //std::cout << P1 << std::endl;
    std::vector<Basis> Bases;
    Bases.push_back(P1);
    LogOut << " DONE!" << std::endl;

    LogOut << "Build Hamiltonian - " << std::flush;
    Holstein<ComplexType> Ham0( Bases );
    Ham0.HolsteinModelK(Kn, Bases, Win, Gin, Jin, N, WithoutK0Phonon);
    LogOut << "Hermitian = " << Ham0.CheckHermitian() << ", Hilbert space = " << Ham0.GetTotalHilbertSpace() << ", DONE!" << std::endl;
    LogOut << "Diagonalize Hamiltonian to find - " << Target << " - " << std::flush;
    RealVectorType Vals;
    ComplexMatrixType Vecs;
    if ( Ham0.GetTotalHilbertSpace() > 1200 ){
      Ham0.eigh(Vals, Vecs, NEV, true, Target);//* ARPACK
    }else{
      Ham0.diag(Vals, Vecs);//* Full spectrum
      NEV = Ham0.GetTotalHilbertSpace();
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
    for ( size_t i = 0; i < NEV; i++ ){
      std::string gname = gnameK;
      gname.append("/S-");
      gname.append( std::to_string( (unsigned long)i) );
      file->SaveNumber(gname, "EVal", Vals[i]);
      ComplexVectorType Vec = Vecs.col(i);
      file->SaveVector(gname, "EVec", Vec);
      ComplexVectorType Npi = NPhonon( Bases, Vec, Ham0);
      file->SaveVector(gname, "Phonon", Npi);
      ComplexVectorType Nfi = NFermion( Bases, Vec, Ham0);
      file->SaveVector(gname, "Fermion", Nfi);
      LogOut << std::setw(3) << i << ": E = " << std::setprecision(12) << std::setw(15) << Vals[i] << ", Np = " << std::setw(14) << arma::accu(Npi).real() << ", Nf = [" << std::flush;
      for ( size_t j = 0; j < L; j++) LogOut << std::setw(16) << Nfi.at(j).real() << ", " << std::flush;
      LogOut << "], Total = " << arma::accu(Nfi).real() << std::endl;
    }
    delete file;
  }
  LogOut.close();
}

void Spectral();

void LoadDynParameters( const std::string filename, int& L, int& N, RealType& G, RealType& W, int& TSteps, RealType& dt){
  HDF5IO h5f(filename);
  h5f.LoadNumber("Parameters", "L", L);
  h5f.LoadNumber("Parameters", "N", N);
  h5f.LoadNumber("Parameters", "G", G);
  h5f.LoadNumber("Parameters", "W", W);
  h5f.LoadNumber("Parameters", "TSteps", TSteps);
  h5f.LoadNumber("Parameters", "dt", dt);
}

void Dynamics(const std::string prefix, const int TargetKIdx, const std::string S1, const std::string S2, const int MeasureEvery, const int SaveWFEvery ){
  std::ofstream LogOut;
  LogOut.open(prefix + "Holstein.K.dyn", std::ios::app);
  int TSteps = 20000;
  RealType dt = 0.005;
  try{
    H5::Exception::dontPrint();
    H5::H5File::isHdf5(prefix + "conf.h5");
    LoadDynParameters( prefix + "conf.h5", L, N, Gin, Win, TSteps, dt);
  }catch(H5::FileIException){
    LogOut << "Use default settings." << std::endl;
  }
  //* k points
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

  int TargetK = Kn.at(TargetKIdx);
  LogOut << "Build Basis - max. phonon quanta = " << N << " with target total momentum K = " << double(TargetK) << ", ends up with Hilbert space = " << std::flush;
  Basis P1(L, N);//! Get rid og k=0 phonon mode, so using L - 1!!
  P1.PhononK(Kn, TargetK, WithoutK0Phonon);
  LogOut << P1.GetHilbertSpace() << std::flush;

  std::vector<Basis> Bases;
  Bases.push_back(P1);
  LogOut << " DONE!" << std::endl;

  LogOut << "Build Hamiltonian - " << std::flush;
  Holstein<ComplexType> Ham0( Bases );
  Ham0.HolsteinModelK(Kn, Bases, Win, Gin, Jin, N, WithoutK0Phonon);
  LogOut << "Hermitian = " << Ham0.CheckHermitian() << ", Hilbert space = " << Ham0.GetTotalHilbertSpace() << ", DONE!" << std::endl;

  LogOut << "Load Wavefunction ... " << std::flush;
  ComplexVectorType VecInit(Ham0.GetTotalHilbertSpace(), arma::fill::randn);
  std::string gnameK = std::to_string( (unsigned long) TargetKIdx );
  std::string SaveFile = "QuenchState-" + gnameK + "-" + S1 + "-" + S2;
  try{
    H5::Exception::dontPrint();
    H5::H5File::isHdf5(prefix + "Holstein.K.h5");
    HDF5IO* file = new HDF5IO(prefix + "Holstein.K.h5");
    ComplexVectorType V1, V2;
    file->LoadVector(gnameK + "/S-" + S1 +"/", "EVec", V1);
    file->LoadVector(gnameK + "/S-" + S2 +"/", "EVec", V2);
    delete file;
    VecInit = ( V1 + V2 );
    LogOut << "<" << S1 << "|" << S2 << "> = " << arma::cdot(V1, V2) << " DONE." << std::endl;
  }catch(H5::FileIException){
    RUNTIME_ERROR("Can not load eigensate from file - Holstein.K.h5. ");
  }
  VecInit = arma::normalise(VecInit);

  ComplexType Prefactor = ComplexType(0.0, -1.0e0*dt);//* NOTE: hbar = 1
  ComplexVectorType VecDyn = VecInit;
  HDF5IO* file2 = new HDF5IO(prefix + SaveFile + ".h5");
  std::string gname = "Obs-0/";
  ComplexType Lecho = arma::cdot(VecInit, VecDyn);
  file2->SaveNumber(gname, "F", Lecho);
  ComplexVectorType Npi = NPhonon( Bases, VecDyn, Ham0);
  file2->SaveVector(gname, "Phonon", Npi);
  ComplexVectorType Nfi = NFermion( Bases, VecDyn, Ham0);
  file2->SaveVector(gname, "Fermion", Nfi);
  delete file2;
  HDF5IO* file3 = new HDF5IO(prefix + SaveFile + "-WF.h5");
  gname = "WF";
  file3->SaveVector(gname, "Vec", VecDyn);
  delete file3;
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
      ComplexVectorType Nfi = NFermion( Bases, VecDyn, Ham0);
      file2->SaveVector(gname, "Fermion", Nfi);
      delete file2;
    }
    if ( cntP % SaveWFEvery == 0 ){
      file3 = new HDF5IO(prefix + SaveFile + "-WF.h5");
      gname = "WF-";
      gname.append( std::to_string((unsigned long long)cntP ));
      gname.append("/");
      file3->SaveVector(gname, "Vec", VecDyn);
      delete file3;
    }
  }
  file3 = new HDF5IO(prefix + SaveFile + "-WF.h5");
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
    int TargetKIdx = std::atoi(argv[2]);
    std::string S1 = argv[3];
    std::string S2 = argv[4];
    int SaveWFEvery = 100, MeasureEvery = 20;
    if ( argc > 5 ) SaveWFEvery = std::atoi(argv[5]);
    if ( argc > 6 ) MeasureEvery = std::atoi(argv[6]);
    Dynamics("", TargetKIdx, S1, S2, MeasureEvery, SaveWFEvery);
  }
  return 0;
}
