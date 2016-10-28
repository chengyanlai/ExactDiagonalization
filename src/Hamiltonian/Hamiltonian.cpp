#include "src/Hamiltonian/Hamiltonian.hpp"
#include "src/Lanczos/krylov.hpp"

template<typename Tnum>
Hamiltonian<Tnum>::Hamiltonian( const std::vector<Basis> &bs )
{
  // get each hilbert space
  HilbertSpaces.clear();
  for ( auto &b : bs ){
    HilbertSpaces.push_back(b.getHilbertSpace());
  }
  size_t TotalDim = getTotalHilbertSpace();
  INFO("Total Hilbert Space = " << TotalDim);
  H_total.resize(TotalDim, TotalDim);
  H_total.reserve(3*TotalDim);
  H_local.resize(TotalDim, TotalDim);
  H_local.reserve(TotalDim);
  H_hop.resize(TotalDim, TotalDim);
  H_hop.reserve(2*TotalDim);
}

template<typename Tnum>
Hamiltonian<Tnum>::~Hamiltonian(){}

template<typename Tnum>
void Hamiltonian<Tnum>::BuildLocalHamiltonian(
  const std::vector< std::vector<Tnum> > &Vloc,
  const std::vector< std::vector<Tnum> > &Uloc,
  const std::vector<Basis> &bs )
{
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
      assert( bs.size() == 1);//NOTE: Only support this right now.
      BosonIntraLocalPart( cnt, Vloc.at(cnt), Uloc.at(cnt), b, hloc );
    }else{//fermion only has potential
      assert( bs.size() < 3 );//NOTE: Only support this right now.
      FermionIntraLocalPart( cnt, Vloc.at(cnt), b, hloc );
    }
    cnt++;
  }
  /* For inter-species local terms*/
  //NOTE: Only support two species fermion
  std::vector<int> sid;
  sid.push_back(0);
  sid.push_back(1);
  if ( bs.size() == 2 && bs.at(0).getType() && bs.at(1).getType() ){
    FermionInterLocalPart(sid, Uloc.at(0).at(0), bs, hloc);
  }
  H_local.setFromTriplets(hloc.begin(), hloc.end());
  std::cout << "Non-zero matrix elements = " << hloc.size() << std::endl;
}

template<typename Tnum>
void Hamiltonian<Tnum>::BuildHoppingHamiltonian(
  const std::vector<Basis> &bs, const std::vector< Node<Tnum>* > &lt )
{
  /* NOTE: This functiuon assume all bases live in the same lattice
  and have the same hopping amplitude. */
  std::vector<Triplet> hhop;
  hhop.clear();
  int cnt = 0;
  for ( auto &b : bs ){
    assert( b.getL() == lt.size() );
    if( !(b.getType()) ){//boson
      assert( bs.size() == 1);//NOTE: Only support this right now.
      BosonIntraHoppingPart( cnt, lt, b, hhop );
    }else{//fermion
      assert( bs.size() == 2);//NOTE: Only support this right now.
      FermionIntraHoppingPart( cnt, lt, b, hhop );
    }
    cnt++;
  }
  H_hop.setFromTriplets(hhop.begin(), hhop.end());
  std::cout << "Non-zero matrix elements = " << hhop.size() << std::endl;
}

template<typename Tnum>
void Hamiltonian<Tnum>::eigh( std::vector<RealType> &Val, VectorType &Vec,
  const bool FullDiagonalization)
{
  size_t dim = getTotalHilbertSpace();
  Vec = VectorType::Random(dim);
  if( !(FullDiagonalization) ){
    Tnum* input_ptr = Vec.data();
    arpackDiagonalize(dim, input_ptr, Val, /*nev*/2, /*tol*/0.0e0);
    Vec = Eigen::Map<VectorType>(input_ptr, dim);
  }
}

template<>
void Hamiltonian<ComplexType>::expH( const ComplexType Prefactor,
  ComplexVectorType &Vec, const size_t Kmax )
{
  krylov(H_total, Vec, Prefactor, Kmax);
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
