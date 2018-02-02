#include "src/Hamiltonian/Hamiltonian.hpp"
#include "src/numeric/lapack.h"

#ifndef DEBUG
#define DEBUG 0
#endif

template<typename Tnum>
Hamiltonian<Tnum>::Hamiltonian( const std::vector<Basis> &bs )
{
  /* Get each hilbert space and calculate total Hilbert space.
     This has nothing to do with Fermion / Boson modeling.
  */
  HilbertSpaces.clear();
  for ( auto &b : bs ){
    HilbertSpaces.push_back(b.getHilbertSpace());
  }
  size_t TotalDim = getTotalHilbertSpace();
  if (DEBUG) std::cout << "Total Hilbert Space = " << TotalDim << std::endl;
  // H_total.resize(TotalDim, TotalDim);// Eigen3
  // H_total.reserve(3*TotalDim);
  // H_local.resize(TotalDim, TotalDim);
  // H_local.reserve(TotalDim);
  // H_hop.resize(TotalDim, TotalDim);
  // H_hop.reserve(2*TotalDim);
  // H_hybridization.resize(TotalDim,TotalDim);
  // H_hybridization.reserve(2*TotalDim);
}

template<typename Tnum>
void Hamiltonian<Tnum>::BuildLocalHamiltonian(const std::vector< std::vector<Tnum> > &Vloc,const std::vector< std::vector<Tnum> > &Uloc,const std::vector<Basis> &bs ){
  ULongMatrixType Locations;
  VectorType Values;
  assert( Vloc.size() == Uloc.size() );
  assert( Vloc.size() == bs.size() );
  int cnt = 0;
  /* For intra-species local terms*/
  for ( auto &b : bs ){
    assert( b.getL() == Vloc.at(cnt).size() );
    assert( b.getL() == Uloc.at(cnt).size() );
    if( !(b.getType()) ){//boson
      if ( bs.size() == 1 ){
        // BosonIntraLocalPart( Vloc.at(cnt), Uloc.at(cnt), b, hloc );// Eigen3
      }else{
        // BosonIntraLocalPart( cnt, Vloc.at(cnt), Uloc.at(cnt), b, hloc );// Eigen3
      }
    }else{//fermion only has potential
      assert( bs.size() < 3 );//NOTE: Only support this right now.
      // FermionIntraLocalPart( cnt, Vloc.at(cnt), b, hloc );// Eigen3
    }
    cnt++;
  }
  /* For inter-species local terms
     NOTE: Only support two species fermion
           due to FermionInterLocalPart
  */
  std::vector<int> sid;
  sid.push_back(0);
  sid.push_back(1);
  if ( bs.size() == 2 && bs.at(0).getType() && bs.at(1).getType() ){
    // FermionInterLocalPart(sid, Uloc.at(0), bs, hloc);// Eigen3
  }else if ( bs.size() == 2 && !(bs.at(0).getType()) && bs.at(1).getType() ){
    /* NOTE: it is hard coded for only OBC now! */
    std::vector<std::tuple<int, int, Tnum> > betweenSitesVals;
    for ( size_t i = 0; i < bs.at(1).getL()-1; i++ ){
      betweenSitesVals.push_back( std::make_tuple(i, i+1, Uloc.at(1).at(i)) );
    }
    // FermionIntraNN(1, betweenSitesVals, bs.at(1), hloc);// Eigen3
  }else if ( bs.size() == 2 && bs.at(0).getType() && !(bs.at(1).getType()) ){
    /* NOTE: Same as above, it is hard coded for only OBC now! */
    std::vector<std::tuple<int, int, Tnum> > betweenSitesVals;
    for ( size_t i = 0; i < bs.at(1).getL()-1; i++ ){
      betweenSitesVals.push_back( std::make_tuple(i, i+1, Uloc.at(0).at(i)) );
    }
    // FermionIntraNN(0, betweenSitesVals, bs.at(0), hloc);// Eigen3
  }
  // H_local.setFromTriplets(hloc.begin(), hloc.end());// Eigen3
  // if (DEBUG) std::cout << "Non-zero matrix elements = " << hloc.size() << std::endl;// Eigen3
}

template<typename Tnum>
void Hamiltonian<Tnum>::BuildHoppingHamiltonian( const std::vector<Basis> &bs, const std::vector< Node<Tnum>* > &lt ){
  /* NOTE: This functiuon assume all bases live in the same lattice
  and have the same hopping amplitude. */
  // std::vector<Triplet> hhop;// Eigen3
  // hhop.clear();
  ULongMatrixType Locations;
  VectorType Values;
  int cnt = 0;
  for ( auto &b : bs ){
    assert( b.getL() == lt.size() );
    if( !(b.getType()) ){//boson
      if ( bs.size() == 1){
        // BosonIntraHoppingPart( lt, b, hhop );// Eigen3
      }else{
        // BosonIntraHoppingPart( cnt, lt, b, hhop );// Eigen3
      }
    }else{//fermion
      assert( bs.size() == 2);//NOTE: Only support this right now.
      // FermionIntraHoppingPart( cnt, lt, b, hhop );// Eigen3
    }
    cnt++;
  }
  // H_hop.setFromTriplets(hhop.begin(), hhop.end());// Eigen3
  // if (DEBUG) std::cout << "Non-zero matrix elements = " << hhop.size() << std::endl;// Eigen3
}

template<typename Tnum>
void Hamiltonian<Tnum>::BuildHoppingHamiltonian( const std::vector<Basis> &bs, const std::vector< std::vector< Node<Tnum>* > > &lt ){
  /* NOTE: This functiuon assume bases live in its own lattice */
  assert( bs.size() == lt.size() );
  assert( bs.size() == 2);//NOTE: Only support this right now due to FermionIntraHoppingPart.
  // std::vector<Triplet> hhop;// Eigen3
  // hhop.clear();
  ULongMatrixType Locations;
  VectorType Values;
  int cnt = 0;
  for ( size_t i = 0; i < bs.size(); i++ ){
    assert( bs.at(i).getL() == lt.at(i).size() );
    if( bs.at(i).getType() ){//fermion
      // FermionIntraHoppingPart( i, lt.at(i), bs.at(i), hhop );// Eigen3
    }else{//boson
      // BosonIntraHoppingPart( i, lt.at(i), bs.at(i), hhop );// Eigen3
    }
    cnt++;
  }
  // H_hop.setFromTriplets(hhop.begin(), hhop.end());// Eigen3
  // if (DEBUG) std::cout << "Non-zero matrix elements = " << hhop.size() << std::endl;// Eigen3
}

// template<typename Tnum>
// void Hamiltonian<Tnum>::BuildXXZHamiltonian(const Tnum Delta,
//   const std::vector<Basis> &bs, const std::vector< Node<Tnum>* > &lt )
// {
//   // std::vector<Triplet> hhop;// Eigen3
//   // hhop.clear();// Eigen3
//   for ( auto &b : bs ){
//     // SpinOneHalfXXZ( Delta, lt, b, hhop );// Eigen3
//   }
//   // H_hop.setFromTriplets(hhop.begin(), hhop.end());// Eigen3
//   // if (DEBUG) std::cout << "Non-zero matrix elements = " << hhop.size() << std::endl;// Eigen3
// }

// template<typename Tnum>
// void Hamiltonian<Tnum>::BuildTIsingHamiltonian(const Tnum hz,
//   const std::vector<Basis> &bs, const std::vector< Node<Tnum>* > &lt )
// {
//   // std::vector<Triplet> hhop;// Eigen3
//   // hhop.clear();// Eigen3
//   for ( auto &b : bs ){
//     // TIsing( hz, lt, b, hhop );// Eigen3
//   }
//   // H_hop.setFromTriplets(hhop.begin(), hhop.end());// Eigen3
//   // if (DEBUG) std::cout << "Non-zero matrix elements = " << hhop.size() << std::endl;// Eigen3
// }

// template<typename Tnum>
// void Hamiltonian<Tnum>::BuildHybridHamiltonian( const int species1, const int species2, const std::vector< std::tuple<int, int, Tnum> > &hybVals, const std::vector<Basis> &bs, const int maxLocalB ){
//   // std::vector<Triplet> hhyb;// Eigen3
//   // Hybridization( species1, species2, hybVals, bs, hhyb, maxLocalB );// Eigen3
//   // H_hybridization.setFromTriplets(hhyb.begin(), hhyb.end());// Eigen3
//   // if (DEBUG) std::cout << "Non-zero matrix elements = " << hhyb.size() << " from H_hyb!" << std::endl;// Eigen3
// }

template<typename Tnum>
void Hamiltonian<Tnum>::eigh( RealVectorType &Vals, MatrixType &Vecs, const int nev, const bool randomInitial){
  size_t dim = getTotalHilbertSpace();
  if ( randomInitial ) Vecs = MatrixType(nev, dim, arma::fill::randu);
  Tnum* input_ptr = Vecs.memptr();
  std::vector<RealType> Val;
  arpackDiagonalize(dim, input_ptr, Val, nev, /*tol*/0.0e0);
  Vecs = MatrixType(input_ptr, nev, dim);
  Vals = RealVectorType(Val);
}

template<typename Tnum>
void Hamiltonian<Tnum>::diag( RealVectorType &Vals, MatrixType &Vecs){
  size_t dim = getTotalHilbertSpace();
  // convert H_total to dense matrix
  MatrixType Mat(H_total);
  // working space
  RealType* EigVec = (T*)malloc( dim * dim * sizeof(T) );
  RealType* Eig = (RealType*)malloc( dim * sizeof(RealType) );
  syDiag(dMat.memptr(), dim, Eig, EigVec);
  Vecs = MatrixType(EigVec, dim, dim);
  Vals = RealVectorType(Eig, dim);
}

template<>
void Hamiltonian<ComplexType>::expH( const ComplexType Prefactor, ComplexVectorType& Vec, const size_t Kmax ){
  RealVectorType KVals(Kmax);// armadillo
  ComplexMatrixType KVecs(Kmax, Vec.n_rows);// armadillo
  KVecs.row(0) = Vec;
  eigh( KVals, KVecs, Kmax, false );
  ComplexMatrixType Dmat(Kmax, Kmax);// armadillo
  for (size_t cnt = 0; cnt < Kmax; cnt++) {
    Dmat(cnt,cnt) = exp( Prefactor * KVals(cnt) );
  }
  ComplexMatrixType KVecsT = KVecs.st();// copy of transpose without conj
  Vec = (KVecsT * Dmat) * (conj(KVecsT) * Vec);// armadillo
}

/* Matrix vector product with MomHamiltonian: y = H_total * x + alpha * y
 * @param x the input vector
 * @param y the output vector
 * @param alpha the scaling value
 */
template<>
void Hamiltonian<RealType>::mvprod(RealType* x, RealType* y, RealType alpha)const {
  size_t dim = getTotalHilbertSpace();
  RealVectorType Vin(x, dim, true, false);//, copy_aux_mem = true, strict = false)// armadillo
  RealVectorType Vout(y, dim, true, false);//, copy_aux_mem = true, strict = false)// armadillo
  Vout = H_total * Vin + alpha * Vout;
  memcpy(y, Vout.memptr(), dim * sizeof(RealType) );// armadillo
}
template<>
void Hamiltonian<ComplexType>::mvprod(ComplexType* x, ComplexType* y, RealType alpha)const {
  size_t dim = getTotalHilbertSpace();
  ComplexVectorType Vin(x, dim, true, false);//, copy_aux_mem = true, strict = false)// armadillo
  ComplexVectorType Vout(y, dim, true, false);//, copy_aux_mem = true, strict = false)// armadillo
  Vout = H_total * Vin + alpha * Vout;
  memcpy(y, Vout.memptr(), dim * sizeof(ComplexType) );// armadillo
}

template class Hamiltonian<RealType>;
template class Hamiltonian<ComplexType>;
