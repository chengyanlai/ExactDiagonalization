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

void LoadEqmParameters( const std::string filename, int& L, int& N, std::vector<RealType>& Momentum, std::vector<RealType>&G, std::vector<RealType>& W){
  HDF5IO h5f(filename);
  h5f.LoadNumber("Parameters", "L", L);
  h5f.LoadNumber("Parameters", "N", N);
  h5f.LoadStdVector("Parameters", "Momentum", Momentum);
  h5f.LoadStdVector("Parameters", "G", G);
  h5f.LoadStdVector("Parameters", "W", W);
}

void Equilibrium(const std::string prefix, const bool SaveBasis = false ){
  /* Holstein model in equilibrium
      This function solve the spectrum of Holstein model and reproduce the result from PRB 60, 1633 (1999).
  */
  std::ofstream LogOut;
  LogOut.open(prefix + "Holstein.1d.eqm.log", std::ios::app);
  const int OBC = 0;
  int N = 9;
  int L = 2 * N;
  std::vector<RealType> Momentum, Win, Gin;
  std::vector<RealType> Jin(L, 1.0);// t0 = 1

  try{
    H5::Exception::dontPrint();
    H5::H5File::isHdf5(prefix + "conf.h5");
    LoadEqmParameters( prefix + "conf.h5", L, N, Momentum, Gin, Win);
  }catch(H5::FileIException){
    Momentum.push_back(0.0);
    // for ( int i = 0; i < 100; i ++ ) Momentum.push_back( 0.01 * RealType(i) );
    // Win = std::vector<RealType>(L,10.00);// omega_0 = 10
    // Gin = std::vector<RealType>(L,20.00);// g = 20; \lambda = g^2 / (2 t0 \omega) = 20
    Win = std::vector<RealType>(L, 1.0);// testing
    Gin = std::vector<RealType>(L, 1.0);// testing - 1
    // Gin = std::vector<RealType>(L, sqrt(2.0));// testing - 1
  }

  // LogOut << "Build Lattice - " << std::endl;
  // std::vector<DT> JWork(Jin.begin(), Jin.end());
  // const std::vector< Node<DT>* > lattice = NN_1D_Chain(L, JWork, OBC);
  // for ( auto &lt : lattice ){
  //   if ( !(lt->VerifySite(LogOut)) ) RUNTIME_ERROR("Wrong lattice setup!");
  // }
  // LogOut << "DONE!" << std::endl;
  LogOut << "Build Basis - " << std::flush;
  Basis B1(L, N);
  std::string BasisFile = "LFS-L";
  BasisFile.append( std::to_string( (unsigned long)L) );
  BasisFile.append( "N" );
  BasisFile.append( std::to_string( (unsigned long)N) );
  BasisFile.append( ".h5" );
  try{
    H5::Exception::dontPrint();
    H5::H5File::isHdf5(prefix + BasisFile);
    B1.Load(prefix + BasisFile, "LFS");
    LogOut << B1.GetHilbertSpace() << " loaded from " << BasisFile << std::flush;
  }catch(H5::FileIException){
    B1.Phonon();
    LogOut << B1.GetHilbertSpace() << std::flush;
    if ( SaveBasis ) B1.Save(prefix + BasisFile, "LFS" );
  }
  std::vector<Basis> Bases;
  Bases.push_back(B1);
  LogOut << " DONE!" << std::endl;
  LogOut << "Build Hamiltonian - " << std::flush;
  Holstein<ComplexType> Ham0( Bases );
  std::vector<ComplexType> Wloc(Win.begin(), Win.end());
  std::vector<ComplexType> Gloc(Gin.begin(), Gin.end());

  for( int i = 0; i < Momentum.size(); i++ ){
    // Ham0.HolsteinModel( Bases, Momentum, lattice,  Wloc,  Gloc );
    Ham0.HolsteinModel( Bases, Momentum.at(i), Jin.at(0),  Wloc,  Gloc );
    // ComplexMatrixType H(Ham0.GetTotalHamiltonian());
    // std::cout << H << std::endl;
    LogOut << "Hermitian = " << Ham0.CheckHermitian() << ", Hilbert space = " << Ham0.GetTotalHilbertSpace() << ", DONE!" << std::endl;
    LogOut << "Diagonalize Hamiltonian - " << std::flush;
    RealVectorType Vals;
    ComplexMatrixType Vecs;
    if ( Ham0.GetTotalHilbertSpace() > 1200 ){
      Ham0.eigh(Vals, Vecs, 10);
      ComplexVectorType Vec = Vecs.col(0);
      HDF5IO* file = new HDF5IO(prefix + "Holstein.h5");
      std::string gname = "M";
      gname.append( std::to_string((unsigned long long)i) );
      file->SaveNumber(gname, "Momentum", Momentum.at(i));
      file->SaveMatrix(gname, "EVEcs", Vecs);
      file->SaveVector(gname, "EVals", Vals);
      file->SaveVector(gname, "GS", Vec);
      delete file;
    }else{
      Ham0.diag(Vals, Vecs);// Full spectrum
      ComplexVectorType Vec = Vecs.col(0);
      std::ofstream SOut;
      SOut.open(prefix + "Holstein.1d.eqm.spectrum.L" + std::to_string((unsigned long)L) + "N" + std::to_string( (unsigned long)N), std::ios::app);
      SOut << std::setprecision(5) << Momentum.at(i) << "\t";
      for ( int j = 0; j < 30 ; j++){
        SOut << std::setprecision(12) << Vals[j] << "\t";
      }
      SOut<< std::endl;
      SOut.close();
    }
    LogOut << "DONE!" << std::endl;
    LogOut << "\tGS energy = " << std::setprecision(12) << Vals[0] << std::endl;
    LogOut << "\tFES energy = " << std::setprecision(12) << Vals[1] << std::endl;
  }
  LogOut << "Finished" << std::endl;
  LogOut.close();
}

void LoadiQDynParameters( const std::string filename, int& L, int& N, std::vector<RealType>& Momentum, std::vector<RealType>&G, std::vector<RealType>& W, int& TSteps, RealType& dt){
  HDF5IO h5f(filename);
  h5f.LoadNumber("Parameters", "L", L);
  h5f.LoadNumber("Parameters", "N", N);
  h5f.LoadStdVector("Parameters", "Momentum", Momentum);
  h5f.LoadStdVector("Parameters", "G", G);
  h5f.LoadStdVector("Parameters", "W", W);
  h5f.LoadNumber("Parameters", "TSteps", TSteps);
  h5f.LoadNumber("Parameters", "dt", dt);
}

void iQDynamics(const std::string prefix, const int MeasureEvery = 20, const int SaveWFEvery = 1000000, const bool SaveBasis = false ){
  /* Interaction quench
      The initial state is k = \pi without phonons.
      At t = 0, the \lambda is turn on and we monitor the dynamics.
      This function reporduce the result in PRB 94 014304 (2016).
  */
  std::ofstream LogOut;
  LogOut.open(prefix + "Holstein.1d.iQdyn.log", std::ios::app);
  const int OBC = 0;
  int N = 16;
  int L = 2 * N;
  RealType Momentum = 1.0;
  std::vector<RealType> WDyn(L, 2.0);// Used in PRB 94 014304 (2016)
  std::vector<RealType> GDyn(L, sqrt(2.0));// Used in PRB 94 014304 (2016)
  std::vector<RealType> Jin(L, 1.0);// t0 = 1
  int TSteps = 10000;
  RealType dt = 0.005;

  try{
    H5::Exception::dontPrint();
    H5::H5File::isHdf5(prefix + "conf.h5");
    std::vector<RealType> Min;
    LoadiQDynParameters( prefix + "conf.h5", L, N, Min, GDyn, WDyn, TSteps, dt);
    Momentum = Min.at(0);
  }catch(H5::FileIException){
    LogOut << " Use default parameters" << std::endl;;
  }

  LogOut << "Build Basis - " << std::flush;
  Basis B1(L, N);
  std::string BasisFile = "LFS-L";
  BasisFile.append( std::to_string( (unsigned long)L) );
  BasisFile.append( "N" );
  BasisFile.append( std::to_string( (unsigned long)N) );
  BasisFile.append( ".h5" );
  try{
    H5::Exception::dontPrint();
    H5::H5File::isHdf5(prefix + BasisFile);
    B1.Load(prefix + BasisFile, "LFS");
    LogOut << B1.GetHilbertSpace() << " loaded from " << BasisFile << std::flush;
  }catch(H5::FileIException){
    B1.Phonon();
    LogOut << B1.GetHilbertSpace() << std::flush;
    if ( SaveBasis ) B1.Save(prefix + BasisFile, "LFS" );
  }
  std::vector<Basis> Bases;
  Bases.push_back(B1);
  LogOut << " DONE!" << std::endl;
  LogOut << "Build Hamiltonian - " << std::flush;
  Holstein<ComplexType> HamDyn( Bases );
  std::vector<ComplexType> Wloc(WDyn.begin(), WDyn.end());
  std::vector<ComplexType> Gloc(GDyn.begin(), GDyn.end());
  HamDyn.HolsteinModel( Bases, Momentum, Jin.at(0),  Wloc,  Gloc );
  LogOut << "Hermitian = " << HamDyn.CheckHermitian() << ", Hilbert space = " << HamDyn.GetTotalHilbertSpace() << ", DONE!" << std::endl;
  LogOut << "Initialize the k = " << Momentum << " state - " << std::flush;
  ComplexVectorType VecT(HamDyn.GetTotalHilbertSpace(), arma::fill::zeros);
  VecT.at(0) = 1.0e0;// k = Momentum state
  LogOut << " DONE!" << std::endl;

  ComplexSparseMatrixType Hk = HamDyn.GetHKinetic();

  LogOut << "Begin dynamics ...... " << std::endl;
  ComplexType Prefactor(0.0, 1.0 * dt);
  ComplexType Ek = arma::cdot(VecT, Hk * VecT);
  HDF5IO* file1 = new HDF5IO(prefix + "Holstein-iQ.h5");
  std::string gname = "obs-0";
  file1->SaveNumber(gname, "Ek", Ek);
  delete file1;
  HDF5IO* file2 = new HDF5IO("Holstein-iQWF.h5");
  gname = "WF-0";
  file2->SaveVector(gname, "Vec", VecT);
  delete file2;
  for ( size_t Ts = 1; Ts <= TSteps; Ts++ ){
    HamDyn.expH(Prefactor, VecT);
    if ( Ts % MeasureEvery == 0 ){
      Ek = arma::cdot(VecT, Hk * VecT);
      file1 = new HDF5IO(prefix + "Holstein-iQ.h5");
      std::string gname = "obs-";
      gname.append( std::to_string((unsigned long long)Ts ));
      gname.append("/");
      file1->SaveNumber(gname, "Ek", Ek);
      delete file1;
    }
    if ( Ts % SaveWFEvery == 0 ){
      file2 = new HDF5IO("Holstein-iQWF.h5");
      gname = "WF-";
      gname.append( std::to_string((unsigned long long)Ts ));
      gname.append("/");
      file2->SaveVector(gname, "Vec", VecT);
      delete file2;
    }
  }
  file2 = new HDF5IO("Holstein-iQWF.h5");
  gname = "WF";
  file2->SaveVector(gname, "Vec", VecT);
  delete file2;
  LogOut << "Finished dynamics!" << std::endl;
  LogOut.close();
}

/* main program */
int main(int argc, char *argv[]){
  if ( argc < 2 ) RUNTIME_ERROR(" Need at least one argument to run program. Use 0 to run Equilibrium.");
#if defined(MKL)
  mkl_set_num_threads(NumCores);
#endif
  if ( std::atoi(argv[1]) == 0 ){
    Equilibrium("");
  }else if ( std::atoi(argv[1]) == 1 ){
    iQDynamics("");
  }
  return 0;
}
