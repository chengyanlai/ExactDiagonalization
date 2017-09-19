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

void LoadParameters( const std::string filename, const int set,
  int &BL, int &FL, int &maxLocalB,
  double &Jbb, double &Jff,
  double &Vbbs, double &Vff,
  double &Uff, double &DeltaDCs){
  HDF5IO file(filename);
  std::string gname = "Input-";
  gname.append(std::to_string((unsigned long long)set));
  BL = file.loadInt(gname, "BL");
  FL = file.loadInt(gname, "FL");
  maxLocalB = file.loadInt(gname, "maxLocalB");
  Jbb = file.loadReal(gname, "Jbb");
  Jff = file.loadReal(gname, "Jff");
  Vbb = file.loadReal(gname, "Vbb");
  Vff = file.loadReal(gname, "Vff");
  Uff = file.loadReal(gname, "Uff");
  DeltaDC = file.loadReal(gname, "DeltaDC");
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
  double DeltaDC = -0.050e0;
#ifdef MPIPARALLEL
  // Initialize MPI
  MPI_Init(NULL, NULL);
  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  LoadParameters( "confs.h5", world_rank, BL, FL, maxLocalB, Jbb, Jff, Vbbs, Vff, Uff, DeltaDC);
#else
  LoadParameters( "confs.h5", 0, BL, FL, maxLocalB, Jbb, Jff, Vbb, Vff, Uff, DeltaDC);
#endif
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
  if ( BL == 2 ) {
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
      DeltaTerm.push_back(std::make_tuple(j, i, DeltaDC));
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
  if ( GetFullSpectrum ){
    Ham.diag(EigVals, EigVecs);
    GS = EigVecs.row(0);
  }else{
    std::vector<RealType> Val;
    Ham.eigh(Val, GS);
    EigVals = dMapVector(&Val[0], 2);
  }
  // std::cout << "Eigen energies = " << EigVals(0) << ", " << EigVals(1) << std::endl;

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
  // std::cout << EigVec.row(0) << std::endl;
  std::vector<std::vector<double> > PeakLocations, PeakWeights;
  for ( size_t cnt = 0; cnt < BL; cnt ++ ){
    std::vector<double> PeakLocation, PeakWeight;
    RealVectorType AS = AdaggerState( cnt, Bs, Ham, GS, maxLocalB);
    if ( GetFullSpectrum ){
      peaks(AS, EigVals, EigVecs, PeakLocation, PeakWeight);
    }else{
      size_t Kmax = 20;
      double threshNorm = 1.0e-12;
      RealVectorType wVals;
      RealMatrixType wVecs;
      Ham.krylovExpansion( AS, wVals, wVecs, Kmax, threshNorm );
      wVecs.transposeInPlace();
      peaks(AS, wVals, wVecs, PeakLocation, PeakWeight);
      for ( auto &val : PeakWeight){
        val *= AS.norm();
      }
    }
    PeakLocations.push_back(PeakLocation);
    PeakWeights.push_back(PeakWeight);
  }

  // save results
  std::string filename = "plex";
#ifdef MPIPARALLEL
  filename.append(std::to_string((unsigned long long)world_rank));
#endif
  filename.append(".h5");
  HDF5IO *file = new HDF5IO(filename);
  file->saveNumber("Input", "BL", BL);
  file->saveNumber("Input", "FL", FL);
  file->saveNumber("Input", "maxLocalB", maxLocalB);
  file->saveNumber("Input", "Jbb", Jbb);
  file->saveNumber("Input", "Jff", Jff);
  file->saveNumber("Input", "Vbb", Vbb);
  file->saveNumber("Input", "Vff", Vff);
  file->saveNumber("Input", "Uff", Uff);
  file->saveNumber("Input", "DeltaDC", DeltaDC);
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
#endif
}
