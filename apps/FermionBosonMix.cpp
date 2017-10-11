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

#ifdef MKL
  #include "mkl.h"
#endif

#ifndef NumCores
#define NumCores 10
#endif

#ifdef MPIPARALLEL
  #include <mpi.h>
#endif

void LoadParameters( const std::string filename,
  int &BL, int &FL, int &maxLocalB,
  double &Jbb, double &Jff,
  double &Vbb, double &Vff,
  double &Uff, std::vector<std::vector<double> > &DeltaDC){
  HDF5IO file(filename);
  std::string gname = "Input";
  BL = file.loadInt(gname, "BL");
  FL = file.loadInt(gname, "FL");
  maxLocalB = file.loadInt(gname, "maxLocalB");
  Jbb = file.loadReal(gname, "Jbb");
  Jff = file.loadReal(gname, "Jff");
  Vbb = file.loadReal(gname, "Vbb");
  Vff = file.loadReal(gname, "Vff");
  Uff = file.loadReal(gname, "Uff");
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

void peaks( const RealVectorType AS, const RealVectorType &EigVal, const RealMatrixType &EigVec,
  std::vector<double> &PeakLocation, std::vector<double> &PeakWeight){
  PeakLocation.clear();
  PeakWeight.clear();
  for ( size_t i = 0; i < EigVec.rows(); i++){
    PeakLocation.push_back(EigVal(i) - EigVal(0));
    RealVectorType An = EigVec.row(i);
    double val = An.dot(AS);
    PeakWeight.push_back( val * val );
  }
  assert( PeakLocation.size() == PeakWeight.size() );
}

int main(int argc, char const *argv[]) {
  Eigen::setNbThreads(NumCores);
#ifdef MKL
  mkl_set_num_threads(NumCores);
#endif
  /* Parameters */
  int BL = 2;
  int FL = 2;
  int maxLocalB = 1;
  double Jbb = 0.020e0;
  double Jff = 0.010e0;
  double Vbb = 3.20e0;
  double Vff = 3.50e0;
  double Uff = 0.00e0;
  std::vector<std::vector<double> > DeltaDC;
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
  assert( argc > 1 );
  size_t RunSets = atoi(argv[1]);
  for ( size_t world_rank = 0; world_rank < RunSets; world_rank++){
#endif
  std::string prefix = "Input-";
  prefix.append(std::to_string((unsigned long long)world_rank));
  prefix.append("/confs.h5");
  LoadParameters( prefix, BL, FL, maxLocalB, Jbb, Jff, Vbb, Vff, Uff, DeltaDC);
  /* Basis */
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
  //   std::cout << BTags.at(i) << ": " << std::flush;
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
  //   std::cout << FTags.at(i) << ": " << std::flush;
  //   Fermions.printFermionBasis( FStates.at(i) );
  //   std::cout << std::endl;
  // }

  /* Hamiltonian */
  // Local potential
  std::vector<std::vector<double> > Vls;
  std::vector<double> Vloc(BL, Vbb);
  Vls.push_back(Vloc);
  Vloc.assign(FL, Vff);
  Vls.push_back(Vloc);
  // Local interaction
  //    Here, U for fermion is NN density density terms.
  std::vector<std::vector<double> > Uls;
  std::vector<double> Uloc(BL, 0.0e0);
  Uls.push_back(Uloc);
  Uloc.assign(FL, Uff);
  Uls.push_back(Uloc);
  Hamiltonian<double> Ham(Bs);
  Ham.BuildLocalHamiltonian(Vls, Uls, Bs);
  // For bosonic chain
  bool BOBC = false;
  std::vector<double> BJ;
  for (size_t cnt = 0; cnt < BL; cnt++) {
    BJ.push_back(Jbb);
  }
  if ( BL < 3 ) {
    BOBC = true;
    BJ.pop_back();
  }
  const std::vector< Node<RealType>* > BLattice = NN_1D_Chain(BL, BJ, BOBC);
  // For fermionic chain
  const bool FOBC = true;
  std::vector<double> FJ;
  for (size_t cnt = 0; cnt < FL-1; cnt++) {
    FJ.push_back(Jff);
  }
  const std::vector< Node<RealType>* > FLattice = NN_1D_Chain(FL, FJ, FOBC);
  std::vector< std::vector< Node<RealType>* > > LT;
  LT.push_back(BLattice);
  LT.push_back(FLattice);
  Ham.BuildHoppingHamiltonian(Bs, LT);
  Ham.BuildTotalHamiltonian();
  std::vector< std::tuple<int, int, double> > DeltaTerm;
  for ( size_t i = 0; i < FL; i++){
    for (size_t j = 0; j < BL; j++){
      DeltaTerm.push_back(std::make_tuple(j, i, DeltaDC.at(i).at(j)));
    }
  }
  Ham.AddHybridHamiltonian( 0, 1, DeltaTerm, Bs, maxLocalB);
  /* Print to check. Only for small matrix */
  // RealSparseMatrixType w2 = Ham.getTotalHamiltonian();
  // RealMatrixType h2(w2);
  // std::cout << h2 << std::endl;

  int GetFullSpectrum = 1;
  if ( Ham.getTotalHilbertSpace() > 5000 ){
    GetFullSpectrum = 0;
  }
  RealVectorType EigVals, GS;
  RealMatrixType EigVecs;
  size_t MaxNumPeak = 20;
  if ( GetFullSpectrum ){
    Ham.diag(EigVals, EigVecs);
  }else{
    Ham.eigh(EigVals, EigVecs, MaxNumPeak);
  }
  GS = EigVecs.row(0);
  // std::cout << EigVals[0] << " " << EigVals[1] << std::endl;

  std::vector<double> Nf, Nb;
  density( Bs, GS, Ham, Nb, Nf);
  // std::cout << "Nb - " << std::endl;
  // for( auto val : Nb){
  //   std::cout << val << " " << std::flush;
  // }
  // std::cout << "\nNf - " << std::endl;
  // for( auto val : Nf){
  //   std::cout << val << " " << std::flush;
  // }
  // std::cout << std::endl;

  std::vector<std::vector<double> > PeakLocations, PeakWeights;
  for ( size_t cnt = 0; cnt < BL; cnt ++ ){
    std::vector<double> PeakLocation, PeakWeight;
    RealVectorType AS = AdaggerState( cnt, Bs, Ham, GS, maxLocalB);
    peaks(AS, EigVals, EigVecs, PeakLocation, PeakWeight);
    for ( int i = 0; i < 4; i++ ) std::cout << PeakLocation.at(i) << " " << PeakWeight.at(i) << std::endl;
    PeakLocations.push_back(PeakLocation);
    PeakWeights.push_back(PeakWeight);
    /* NOTE: Why this krylovExpansion is not working? */
    // double threshNorm = 1.0e-12;
    // RealVectorType wVals;
    // RealMatrixType wVecs;
    // Ham.krylovExpansion( AS, wVals, wVecs, MaxNumPeak, threshNorm );
    // std::cout << wVals[0] << " " << wVals[1] << std::endl;
    // std::cout << wVecs.rows() << " " << wVecs.cols() << std::endl;
    // wVecs.transposeInPlace();
    // peaks(AS, wVals, wVecs, PeakLocation, PeakWeight);
    // std::cout << wVecs.rows() << " " << wVecs.cols() << std::endl;
    // std::cout << "cp" << std::endl;
    // for ( int i = 0; i < 4; i++ ) std::cout << PeakLocation.at(i) << " " << PeakWeight.at(i) << std::endl;
  }

  // save results
  std::string filename = "Input-";
  filename.append(std::to_string((unsigned long long)world_rank));
  filename.append("/plex.h5");
  HDF5IO *file = new HDF5IO(filename);
  file->saveNumber("Input", "BL", BL);
  file->saveNumber("Input", "FL", FL);
  file->saveNumber("Input", "maxLocalB", maxLocalB);
  file->saveNumber("Input", "Jbb", Jbb);
  file->saveNumber("Input", "Jff", Jff);
  file->saveNumber("Input", "Vbb", Vbb);
  file->saveNumber("Input", "Vff", Vff);
  file->saveNumber("Input", "Uff", Uff);
  std::string setName = "DeltaDC-";
  for ( size_t i = 0; i < FL; i++){
    std::string tmp = setName;
    tmp.append(std::to_string((unsigned long long)i));
    file->saveStdVector("Input", tmp, DeltaDC.at(i));
  }
  std::string gname = "obs";
  file->saveNumber(gname, "FullSpectrum", GetFullSpectrum);
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
#ifdef MPIPARALLEL
  MPI_Finalize();
#else
  }
#endif
}
