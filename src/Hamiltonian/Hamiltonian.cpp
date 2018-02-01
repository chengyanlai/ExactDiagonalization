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
  H_total.resize(TotalDim, TotalDim);
  H_total.reserve(3*TotalDim);
  H_local.resize(TotalDim, TotalDim);
  H_local.reserve(TotalDim);
  H_hop.resize(TotalDim, TotalDim);
  H_hop.reserve(2*TotalDim);
  H_hybridization.resize(TotalDim,TotalDim);
  H_hybridization.reserve(2*TotalDim);
}

template<typename Tnum>
void Hamiltonian<Tnum>::BuildLocalHamiltonian(const std::vector< std::vector<Tnum> > &Vloc,const std::vector< std::vector<Tnum> > &Uloc,const std::vector<Basis> &bs ){
  std::vector<Triplet> hloc;
  hloc.clear();
  assert( Vloc.size() == Uloc.size() );
  assert( Vloc.size() == bs.size() );
  int cnt = 0;
  /* For intra-species local terms*/
  for ( auto &b : bs ){
    assert( b.getL() == Vloc.at(cnt).size() );
    assert( b.getL() == Uloc.at(cnt).size() );
    if( !(b.getType()) ){//boson
      if ( bs.size() == 1 ){
        BosonIntraLocalPart( Vloc.at(cnt), Uloc.at(cnt), b, hloc );
      }else{
        BosonIntraLocalPart( cnt, Vloc.at(cnt), Uloc.at(cnt), b, hloc );
      }
    }else{//fermion only has potential
      assert( bs.size() < 3 );//NOTE: Only support this right now.
      FermionIntraLocalPart( cnt, Vloc.at(cnt), b, hloc );
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
    FermionInterLocalPart(sid, Uloc.at(0), bs, hloc);
    // FermionInterLocalPart(sid, Uloc.at(0).at(0), bs, hloc);
  }else if ( bs.size() == 2 && !(bs.at(0).getType()) && bs.at(1).getType() ){
    /* NOTE: it is hard coded for only OBC now! */
    std::vector<std::tuple<int, int, Tnum> > betweenSitesVals;
    for ( size_t i = 0; i < bs.at(1).getL()-1; i++ ){
      betweenSitesVals.push_back( std::make_tuple(i, i+1, Uloc.at(1).at(i)) );
    }
    FermionIntraNN(1, betweenSitesVals, bs.at(1), hloc);
  }else if ( bs.size() == 2 && bs.at(0).getType() && !(bs.at(1).getType()) ){
    /* NOTE: Same as above, it is hard coded for only OBC now! */
    std::vector<std::tuple<int, int, Tnum> > betweenSitesVals;
    for ( size_t i = 0; i < bs.at(1).getL()-1; i++ ){
      betweenSitesVals.push_back( std::make_tuple(i, i+1, Uloc.at(0).at(i)) );
    }
    FermionIntraNN(0, betweenSitesVals, bs.at(0), hloc);
  }
  H_local.setFromTriplets(hloc.begin(), hloc.end());
  if (DEBUG) std::cout << "Non-zero matrix elements = " << hloc.size() << std::endl;
}

template<typename Tnum>
void Hamiltonian<Tnum>::BuildHoppingHamiltonian( const std::vector<Basis> &bs, const std::vector< Node<Tnum>* > &lt ){
  /* NOTE: This functiuon assume all bases live in the same lattice
  and have the same hopping amplitude. */
  std::vector<Triplet> hhop;
  hhop.clear();
  int cnt = 0;
  for ( auto &b : bs ){
    assert( b.getL() == lt.size() );
    if( !(b.getType()) ){//boson
      if ( bs.size() == 1){
        BosonIntraHoppingPart( lt, b, hhop );
      }else{
        BosonIntraHoppingPart( cnt, lt, b, hhop );
      }
    }else{//fermion
      assert( bs.size() == 2);//NOTE: Only support this right now.
      FermionIntraHoppingPart( cnt, lt, b, hhop );
    }
    cnt++;
  }
  H_hop.setFromTriplets(hhop.begin(), hhop.end());
  if (DEBUG) std::cout << "Non-zero matrix elements = " << hhop.size() << std::endl;
}

template<typename Tnum>
void Hamiltonian<Tnum>::BuildHoppingHamiltonian( const std::vector<Basis> &bs, const std::vector< std::vector< Node<Tnum>* > > &lt ){
  /* NOTE: This functiuon assume bases live in its own lattice */
  assert( bs.size() == lt.size() );
  assert( bs.size() == 2);//NOTE: Only support this right now due to FermionIntraHoppingPart.
  std::vector<Triplet> hhop;
  hhop.clear();
  int cnt = 0;
  for ( size_t i = 0; i < bs.size(); i++ ){
    assert( bs.at(i).getL() == lt.at(i).size() );
    if( bs.at(i).getType() ){//fermion
      FermionIntraHoppingPart( i, lt.at(i), bs.at(i), hhop );
    }else{//boson
      BosonIntraHoppingPart( i, lt.at(i), bs.at(i), hhop );
    }
    cnt++;
  }
  H_hop.setFromTriplets(hhop.begin(), hhop.end());
  if (DEBUG) std::cout << "Non-zero matrix elements = " << hhop.size() << std::endl;
}

template<typename Tnum>
void Hamiltonian<Tnum>::BuildXXZHamiltonian(const Tnum Delta,
  const std::vector<Basis> &bs, const std::vector< Node<Tnum>* > &lt )
{
  std::vector<Triplet> hhop;
  hhop.clear();
  for ( auto &b : bs ){
    SpinOneHalfXXZ( Delta, lt, b, hhop );
  }
  H_hop.setFromTriplets(hhop.begin(), hhop.end());
  if (DEBUG) std::cout << "Non-zero matrix elements = " << hhop.size() << std::endl;
}

template<typename Tnum>
void Hamiltonian<Tnum>::BuildTIsingHamiltonian(const Tnum hz,
  const std::vector<Basis> &bs, const std::vector< Node<Tnum>* > &lt )
{
  std::vector<Triplet> hhop;
  hhop.clear();
  for ( auto &b : bs ){
    TIsing( hz, lt, b, hhop );
  }
  H_hop.setFromTriplets(hhop.begin(), hhop.end());
  if (DEBUG) std::cout << "Non-zero matrix elements = " << hhop.size() << std::endl;
}

template<typename Tnum>
void Hamiltonian<Tnum>::BuildHybridHamiltonian( const int species1, const int species2, const std::vector< std::tuple<int, int, Tnum> > &hybVals, const std::vector<Basis> &bs, const int maxLocalB ){
  std::vector<Triplet> hhyb;
  Hybridization( species1, species2, hybVals, bs, hhyb, maxLocalB );
  H_hybridization.setFromTriplets(hhyb.begin(), hhyb.end());
  if (DEBUG) std::cout << "Non-zero matrix elements = " << hhyb.size() << " from H_hyb!" << std::endl;
}

template<>
void Hamiltonian<RealType>::eigh( RealVectorType &Vals, RealMatrixType &Vecs, const int nev, const bool randomInitial){
  size_t dim = getTotalHilbertSpace();
  if ( randomInitial ) Vecs = MatrixType::Random(nev, dim);
  RealType* input_ptr = Vecs.data();
  std::vector<RealType> Val;
  arpackDiagonalize(dim, input_ptr, Val, nev, /*tol*/0.0e0);
  Vals = dMapVector(&Val[0], nev);
  Vecs = MapMatrix(input_ptr, nev, dim);
}

template<>
void Hamiltonian<ComplexType>::eigh( RealVectorType &Vals, ComplexMatrixType &Vecs, const int nev, const bool randomInitial){
  size_t dim = getTotalHilbertSpace();
  if ( randomInitial ) Vecs = MatrixType::Random(nev, dim);
  ComplexType* input_ptr = Vecs.data();
  std::vector<RealType> Val;
  arpackDiagonalize(dim, input_ptr, Val, nev, /*tol*/0.0e0);
  Vals = dMapVector(&Val[0], nev);
  Vecs = MapMatrix(input_ptr, nev, dim);
}

template<>
void Hamiltonian<RealType>::diag( RealVectorType &Vals, RealMatrixType &Vecs){
  size_t dim = getTotalHilbertSpace();
  // convert H_total to dense matrix
  RealMatrixType dMat = MatrixType(H_total);
  // working space
  RealType* EigVec = (RealType*)malloc( dim * dim * sizeof(RealType) );
  RealType* Eig = (RealType*)malloc( dim * sizeof(RealType) );
  syDiag(dMat.data(), dim, Eig, EigVec);
  Vecs = MapMatrix(EigVec, dim, dim);
  Vals = dMapVector(Eig, dim);
}

template<>
void Hamiltonian<ComplexType>::diag( RealVectorType &Vals, ComplexMatrixType &Vecs){
  size_t dim = getTotalHilbertSpace();
  // convert H_total to dense matrix
  ComplexMatrixType dMat = MatrixType(H_total);
  // working space
  ComplexType* EigVec = (ComplexType*)malloc( dim * dim * sizeof(ComplexType) );
  RealType* Eig = (RealType*)malloc( dim * sizeof(RealType) );
  syDiag(dMat.data(), dim, Eig, EigVec);
  Vecs = MapMatrix(EigVec, dim, dim);
  Vals = dMapVector(Eig, dim);
}

template<>
void Hamiltonian<ComplexType>::expH( const ComplexType Prefactor, ComplexVectorType &Vec, const size_t Kmax ){
  RealVectorType KVals = RealVectorType::Zero(Kmax);
  ComplexMatrixType KVecs = ComplexMatrixType::Zero(Kmax, Vec.rows());
  KVecs.row(0) = Vec;
  std::cout << Vec << std::endl;
  std::cout << KVecs.data()[0] << " " << KVecs.data()[1] << " " << KVecs.data()[2] << " " << std::endl;
  std::cin.get();
  eigh( KVals, KVecs, Kmax, false);
  ComplexMatrixType Dmat = ComplexMatrixType::Zero(Kmax, Kmax);
  for (size_t cnt = 0; cnt < Kmax; cnt++) {
    Dmat(cnt,cnt) = exp( Prefactor * KVals(cnt) );
  }
  ComplexVectorType work = Dmat * (KVecs * Vec);
  Vec = KVecs.adjoint() * work;
}

template<>
RealVectorType Hamiltonian<RealType>::expVals( const RealType Prefactor, const RealVectorType Vec){
  RealVectorType out = Prefactor * Vec;
  Eigen::ArrayXd work = out.array().exp();
  return work.matrix();
}

/* Matrix vector product with MomHamiltonian: y = H_total * x + alpha * y
 * @param x the input vector
 * @param y the output vector
 * @param alpha the scaling value
 */
template<>
void Hamiltonian<RealType>::mvprod(RealType* x, RealType* y, RealType alpha)const {
  size_t dim = getTotalHilbertSpace();
  Eigen::Map<RealVectorType> Vin(x, dim);
  Eigen::Map<RealVectorType> Vout(y, dim);
  Vout = H_total * Vin + alpha * Vout;
  memcpy(y, Vout.data(), dim * sizeof(RealType) );
}
template<>
void Hamiltonian<ComplexType>::mvprod(ComplexType* x, ComplexType* y, RealType alpha)const {
  size_t dim = getTotalHilbertSpace();
  Eigen::Map<ComplexVectorType> Vin(x, dim);
  Eigen::Map<ComplexVectorType> Vout(y, dim);
  Vout = H_total * Vin + alpha * Vout;
  memcpy(y, Vout.data(), dim * sizeof(ComplexType) );
}

template class Hamiltonian<RealType>;
template class Hamiltonian<ComplexType>;
