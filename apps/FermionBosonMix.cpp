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

#ifndef NumCores
#define NumCores 10
#endif

#ifdef MPIPARALLEL
  #include <mpi.h>
#endif

std::vector<std::string> GetPrefix(const std::string FileName){
  std::ifstream file(FileName);
  std::vector<std::string> out;
  std::string s;
  while (std::getline(file, s)){
    out.push_back(s + "/");
  }
  return out;
}

void LoadPulse(const std::string filename, int& Tf, double& dT, std::vector<double>& Ect, int& SaveWFEvery){
  HDF5IO file(filename);
  std::string gname = "Input";
  Tf = file.loadInt(gname, "Tf");
  dT = file.loadReal(gname, "dT");
  Ect.clear();
  file.loadStdVector(gname, "Ect", Ect);
  // SaveWFEvery = file.loadInt(gname, "SaveWFEvery");
  SaveWFEvery = 40;
}

void LoadParameters( const std::string filename, int &BL, int &FL, int &maxLocalB,
  double &Jd, double &Jc, double &Ed, double &Ec, double &Vc, std::vector<std::vector<double> > &DeltaDC){
  HDF5IO file(filename);
  std::string gname = "Input";
  BL = file.loadInt(gname, "BL");
  FL = file.loadInt(gname, "FL");
  maxLocalB = file.loadInt(gname, "maxLocalB");
  Jd = file.loadReal(gname, "Jd");
  Jc = file.loadReal(gname, "Jc");
  Ed = file.loadReal(gname, "Ed");
  Ec = file.loadReal(gname, "Ec");
  Vc = file.loadReal(gname, "Vc");
  DeltaDC.clear();
  std::string setName = "DeltaDC-";
  for ( size_t i = 0; i < FL; i++){
    std::string tmp = setName;
    tmp.append(std::to_string((unsigned long long)i));
    std::vector<double> work;
    work.clear();
    file.loadStdVector(gname, tmp, work);
    DeltaDC.push_back( work );
  }
}

void density( const std::vector<Basis> &bs, const RealVectorType GS, Hamiltonian<double> Ham,
  std::vector<double> &nb, std::vector<double> &nf){
  nb.clear();
  nb.assign(bs.at(0).getL(), 0.0e0);
  nf.clear();
  nf.assign(bs.at(1).getL(), 0.0e0);
  size_t bid = 0;
  for ( std::vector<int> b : bs.at(0).getBStates() ){
    size_t fid = 0;
    for ( int f : bs.at(1).getFStates() ){
      for ( size_t fl = 0; fl < bs.at(1).getL(); fl++ ){
        if ( btest(f, fl) ){
          size_t id = Ham.DetermineTotalIndex(vec<size_t>(bid, fid));
          nf.at(fl) += GS(id) * GS(id);
        }
      }
      for ( size_t bl = 0; bl < bs.at(0).getL(); bl++ ){
        if ( b.at(bl) > 0 ){
          size_t id = Ham.DetermineTotalIndex(vec<size_t>(bid, fid));
          nb.at(bl) += GS(id) * GS(id);
        }
      }
      fid++;
    }
    bid++;
  }
}

RealVectorType AdaggerState( const int site, const std::vector<Basis> &bs,
  Hamiltonian<double> Ham, const RealVectorType State, const int maxLocalB ){
  assert( !(bs.at(0).getType()) );
  RealVectorType psi = RealVectorType::Zero(State.rows());
  size_t bid1 = 0;
  for ( std::vector<int> b : bs.at(0).getBStates() ){
    std::vector<int> nb = b;
    if( (b.at(site) < maxLocalB && maxLocalB) || (b.at(site) < bs.at(0).getN() && !(maxLocalB)) ){
      nb.at(site) += 1;
      size_t bid2 = bs.at(0).getIndexFromTag( BosonBasisTag(nb) );
      for ( size_t j = 0; j < bs.at(1).getHilbertSpace(); j++ ){
        size_t id1 = Ham.DetermineTotalIndex(vec<size_t>(bid1, j));
        size_t id2 = Ham.DetermineTotalIndex(vec<size_t>(bid2, j));
        psi(id2) = State(id1);
      }
    }
    bid1++;
  }
  return psi;
}

RealVectorType AState( const int site, const std::vector<Basis> &bs,
  Hamiltonian<double> Ham, const RealVectorType State, const int maxLocalB ){
  assert( !(bs.at(0).getType()) );
  RealVectorType psi = RealVectorType::Zero(State.rows());
  size_t bid1 = 0;
  for ( std::vector<int> b : bs.at(0).getBStates() ){
    std::vector<int> nb = b;
    if( b.at(site) > 0 && maxLocalB ){
      nb.at(site) -= 1;
      size_t bid2 = bs.at(0).getIndexFromTag( BosonBasisTag(nb) );
      for ( size_t j = 0; j < bs.at(1).getHilbertSpace(); j++ ){
        size_t id1 = Ham.DetermineTotalIndex(vec<size_t>(bid1, j));
        size_t id2 = Ham.DetermineTotalIndex(vec<size_t>(bid2, j));
        psi(id2) = State(id1);
      }
    }
    bid1++;
  }
  return psi;
}

void peaks( const double EG, const RealVectorType AS, const RealVectorType &EigVal, const RealMatrixType &EigVec,
  std::vector<double> &PeakLocation, std::vector<double> &PeakWeight, const int MaxNumPeak){
  PeakLocation.clear();
  PeakWeight.clear();
  for ( size_t i = 0; i < EigVec.rows(); i++){
    PeakLocation.push_back(EigVal(i) - EG);
    RealVectorType An = EigVec.row(i);
    double val = An.dot(AS);
    PeakWeight.push_back( val * val );
    if ( i > MaxNumPeak ) break;
  }
  assert( PeakLocation.size() == PeakWeight.size() );
}

std::vector<Basis> BuildBasis(const int& BL, const int& FL, const int& maxLocalB){
  std::vector<Basis> Bs;
  // For bosons
  int BN = BL;
  Basis Bosons(BL, BN);
  Bosons.Boson(maxLocalB);
  Bs.push_back(Bosons);
  /* Print to check. Only for small matrix */
  // std::vector< std::vector<int> > BStates = Bosons.getBStates();
  // std::vector<RealType> BTags = Bosons.getBTags();
  // for( size_t i = 0; i < BTags.size(); ++i ){
  //   OASOut << BTags.at(i) << ": " << std::flush;
  //   Bosons.printBosonBasis( BStates.at(i) );
  // }

  // For fermions
  int FN = FL;
  Basis Fermions(FL, FN, true);
  Fermions.Fermion(0);
  Bs.push_back(Fermions);
  /* Print to check. Only for small matrix */
  // std::vector<int> FStates = Fermions.getFStates();
  // std::vector<size_t> FTags = Fermions.getFTags();
  // for( size_t i = 0; i < FTags.size(); ++i ){
  //   OASOut << FTags.at(i) << ": " << std::flush;
  //   Fermions.printFermionBasis( FStates.at(i) );
  //   OASOut << std::endl;
  // }
  return Bs;
}

template<typename T>
std::vector< std::vector< Node<T>* > > BuildLattice( const int& BL, const T& Jd, const int& FL, const T& Jc){
  // For bosonic chain
  bool BOBC = false;
  std::vector<T> BJ;
  for (size_t cnt = 0; cnt < BL; cnt++) {
    BJ.push_back(Jd);
  }
  if ( BL < 3 ) {
    BOBC = true;
    BJ.pop_back();
  }
  const std::vector< Node<T>* > BLattice = NN_1D_Chain(BL, BJ, BOBC);
  // For fermionic chain
  const bool FOBC = true;
  std::vector<T> FJ;
  for (size_t cnt = 0; cnt < FL-1; cnt++) {
    FJ.push_back(Jc);
  }
  const std::vector< Node<T>* > FLattice = NN_1D_Chain(FL, FJ, FOBC);
  std::vector< std::vector< Node<T>* > > LT;
  LT.push_back(BLattice);
  LT.push_back(FLattice);
  return LT;
}

template<typename T>
Hamiltonian<T> BuildHamiltonian( const int& maxLocalB, const std::vector<Basis>& Bs, const double& Ed, const double& Ec, const double& Vc, std::vector< std::vector< Node<T>* > > LT, const std::vector<std::vector<double> >& DeltaDC ){
  Hamiltonian<T> Ham(Bs);
  // Local potential
  std::vector<std::vector<T> > Vls;
  std::vector<T> Vloc(Bs.at(0).getL(), T(Ed));
  Vls.push_back(Vloc);
  Vloc.assign(Bs.at(1).getL(), T(Ec));
  Vls.push_back(Vloc);
  // Local interaction
  //    Here, U for fermion is NN density density terms.
  std::vector<std::vector<T> > Uls;
  std::vector<T> Uloc(Bs.at(0).getL(), T(0.0e0));
  Uls.push_back(Uloc);
  Uloc.assign(Bs.at(1).getL(), T(Vc));
  Uls.push_back(Uloc);
  Ham.BuildLocalHamiltonian(Vls, Uls, Bs);
  Ham.BuildHoppingHamiltonian(Bs, LT);
  std::vector< std::tuple<int, int, T> > DeltaTerm;
  for ( size_t i = 0; i < Bs.at(1).getL(); i++){
    for (size_t j = 0; j < Bs.at(0).getL(); j++){
      DeltaTerm.push_back(std::make_tuple(j, i, T(DeltaDC.at(i).at(j))));
    }
  }
  Ham.BuildHybridHamiltonian( 0, 1, DeltaTerm, Bs, maxLocalB);
  Ham.BuildTotalHamiltonian();
  /* Print to check. Only for small matrix */
  // RealSparseMatrixType w2 = Ham.getTotalHamiltonian();
  // RealMatrixType h2(w2);
  // OASOut << h2 << std::endl;
  return Ham;
}

void OAS(const std::string prefix, const int dynamics=0){
  std::ofstream OASOut;
  OASOut.open(prefix + "trace.oas", std::ios::app);
  /* Parameters */
  int BL = 2;
  int FL = 2;
  int maxLocalB = 1;
  double Jd = 0.020e0;
  double Jc = 0.010e0;
  double Ed = 3.20e0;
  double Ec = 3.50e0;
  double Vc = 0.00e0;
  std::vector<std::vector<double> > DeltaDC;

  std::string filename = prefix;
  filename.append("confs.h5");
  LoadParameters( filename, BL, FL, maxLocalB, Jd, Jc, Ed, Ec, Vc, DeltaDC);
  /* Basis */
  std::vector<Basis> Bs = BuildBasis(BL, FL, maxLocalB);

  /* Hamiltonian */
  std::vector< std::vector< Node<RealType>* > > LT = BuildLattice( BL, Jd, FL, Jc);
  Hamiltonian<double> Ham = BuildHamiltonian( maxLocalB, Bs, Ed, Ec, Vc, LT, DeltaDC );

  /* Get ground state and excited states */
  int GetFullSpectrum = 1;
  if ( Ham.getTotalHilbertSpace() > 5000 ){
    GetFullSpectrum = 0;
  }
  RealVectorType EigVals, GS;
  RealMatrixType EigVecs;
  if ( GetFullSpectrum ){
    Ham.diag(EigVals, EigVecs);
  }else{
    Ham.eigh(EigVals, EigVecs, 2);
  }
  GS = EigVecs.row(0);
  double EG = EigVals[0];
  OASOut << EigVals[0] << " " << EigVals[1] << std::endl;
  std::vector<double> Nf, Nb;
  density( Bs, GS, Ham, Nb, Nf);

  /* Get equilibrium spectrum */
  size_t MaxNumPeak;
  if ( GetFullSpectrum ) MaxNumPeak =  Ham.getTotalHilbertSpace();
  else MaxNumPeak = 20;
  std::vector<std::vector<double> > PeakLocations, PeakWeights;
  for ( size_t cnt = 0; cnt < BL; cnt ++ ){
    std::vector<double> PeakLocation, PeakWeight;
    RealVectorType AS = AdaggerState( cnt, Bs, Ham, GS, maxLocalB);
    if ( !(GetFullSpectrum) ){
      EigVecs.row(0) = AS;
      Ham.eigh(EigVals, EigVecs, MaxNumPeak, false);
      OASOut << "new " << EigVals[0] << " " << EigVals[1] << std::endl;
    }
    peaks(EG, AS, EigVals, EigVecs, PeakLocation, PeakWeight, MaxNumPeak);
    /* Print to check. */
    for ( int i = 0; i < 4; i++ ) OASOut << PeakLocation.at(i) << " " << PeakWeight.at(i) << std::endl;
    PeakLocations.push_back(PeakLocation);
    PeakWeights.push_back(PeakWeight);
  }

  /* save results */
  filename = prefix;
  filename.append("plex.h5");
  HDF5IO *file = new HDF5IO(filename);
  file->saveNumber("Input", "BL", BL);
  file->saveNumber("Input", "FL", FL);
  file->saveNumber("Input", "maxLocalB", maxLocalB);
  file->saveNumber("Input", "Jd", Jd);
  file->saveNumber("Input", "Jc", Jc);
  file->saveNumber("Input", "Ed", Ed);
  file->saveNumber("Input", "Ec", Ec);
  file->saveNumber("Input", "Vc", Vc);
  std::string setName = "DeltaDC-";
  for ( size_t i = 0; i < FL; i++){
    std::string tmp = setName;
    tmp.append(std::to_string((unsigned long long)i));
    file->saveStdVector("Input", tmp, DeltaDC.at(i));
  }
  std::string gname = "obs";
  file->saveNumber(gname, "FullSpectrum", GetFullSpectrum);
  file->saveNumber(gname, "EG", EG);
  file->saveVector(gname, "EigVals", EigVals);
  file->saveMatrix(gname, "EigVecs", EigVecs);
  file->saveStdVector(gname, "Nf", Nf);
  file->saveStdVector(gname, "Nb", Nb);
  std::string gname2 = "peak";
  for ( size_t cnt = 0; cnt < BL; cnt ++ ){
    std::string dname1 = "Ld-";
    dname1.append(std::to_string((unsigned long long)cnt));
    file->saveStdVector(gname2, dname1, PeakLocations.at(cnt));
    std::string dname2 = "Wd-";
    dname2.append(std::to_string((unsigned long long)cnt));
    file->saveStdVector(gname2, dname2, PeakWeights.at(cnt));
  }
  delete file;

  /* Run dynamics */
  if ( dynamics ){
    int Tf, SaveWFEvery;
    int Kmax = 20;
    double dT;
    std::vector<double> Ect;
    LoadPulse(prefix+"pulse.h5", Tf, dT, Ect, SaveWFEvery);
    ComplexType Prefactor = ComplexType(0.0, -1.0 * dT );
    std::vector< std::vector< Node<ComplexType>* > > LTc = BuildLattice( BL, ComplexType(Jd, 0.0e0), FL, ComplexType(Jc, 0.0e0) );

    /* save inputs */
    filename = prefix;
    filename.append("trplex.h5");
    HDF5IO *file1 = new HDF5IO(filename);
    file1->saveNumber("Input", "BL", BL);
    file1->saveNumber("Input", "FL", FL);
    file1->saveNumber("Input", "maxLocalB", maxLocalB);
    file1->saveNumber("Input", "Jd", Jd);
    file1->saveNumber("Input", "Jc", Jc);
    file1->saveNumber("Input", "Ed", Ed);
    file1->saveStdVector("Input", "Ect", Ect);
    file1->saveNumber("Input", "Vc", Vc);
    file1->saveNumber("Input", "dT", dT);
    file1->saveNumber("Input", "Tf", Tf);
    for ( size_t i = 0; i < FL; i++){
      std::string tmp1 = "DeltaDC-";
      tmp1.append(std::to_string((unsigned long long)i));
      file1->saveStdVector("Input", tmp1, DeltaDC.at(i));
    }
    delete file1;

    for ( size_t cnt = 0; cnt < BL; cnt ++ ){
      RealVectorType AS = AdaggerState( cnt, Bs, Ham, GS, maxLocalB );
      ComplexVectorType VecT = AS.cast<ComplexType>();
      size_t TStep = 0;
      Hamiltonian<ComplexType> HamT;
      std::vector<ComplexType> At;
      At.clear();
      while ( TStep < Tf ){
        if ( TStep < Ect.size() ){// update Hamiltonian
          HamT = BuildHamiltonian( maxLocalB, Bs, Ed, Ect.at(TStep), Vc, LTc, DeltaDC );
        }
        HamT.expH( Prefactor, VecT, Kmax );
        At.push_back( AS.dot(VecT) );
        TStep++;
        if ( TStep % SaveWFEvery ){
          filename = prefix;
          filename.append("trplex.h5");
          HDF5IO *file3 = new HDF5IO(filename);
          std::string gname3 = "wf";
          std::string sname3 = "vec-";
          sname3.append(std::to_string((unsigned long long)TStep));
          file3->saveVector(gname3, sname3, VecT);
          delete file3;
        }
      }
      filename = prefix;
      filename.append("trplex.h5");
      HDF5IO *file2 = new HDF5IO(filename);
      std::string gname2 = "obs";
      std::string sname2 = "at-";
      sname2.append(std::to_string((unsigned long long)cnt));
      file2->saveStdVector(gname2, sname2, At);
      delete file2;
    }
  }
}

/* main program */
int main(int argc, char *argv[]){
  if ( argc < 2 ) RUNTIME_ERROR(" Need at least one argument to run program. ");
#ifdef MKL
  mkl_set_dynamic(0);
  mkl_set_num_threads(NumCores);
#endif
  int world_size;
  int world_rank;
  std::vector<std::string> MPIFolders = GetPrefix("MPIFolders");
#ifdef MPIPARALLEL
  // Initialize MPI
  MPI_Init(NULL, NULL);
  // Get the number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  // Get the rank of the process
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  assert ( MPIFolders.size() == world_size );
  OAS( MPIFolders.at(world_rank), std::atoi(argv[1]) );
  MPI_Finalize();
#else
  world_size = MPIFolders.size();
  world_rank = 0;
  for(;world_rank < world_size; world_rank++){
    OAS( MPIFolders.at(world_rank), std::atoi(argv[1]) );
  }
#endif
  return 0;
}
