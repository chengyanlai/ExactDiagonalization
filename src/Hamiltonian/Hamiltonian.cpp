#include "src/Hamiltonian/Hamiltonian.hpp"
#include "src/Lanczos/lanczos.hpp"
#include "src/Lanczos/krylov.hpp"

template<typename Tnum, typename Tlabel>
Hamiltonian<Tnum, Tlabel>::Hamiltonian( const std::vector<Basis> &bs )
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

template<typename Tnum, typename Tlabel>
Hamiltonian<Tnum, Tlabel>::~Hamiltonian(){}

template<typename Tnum, typename Tlabel>
void Hamiltonian<Tnum, Tlabel>::BuildLocalHamiltonian(
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
  /*test**/
  // hloc.push_back(hloc.at(0));
  // hloc.push_back(hloc.at(0));
  // hloc.push_back(hloc.at(0));
  // INFO(hloc.size());
  H_local.setFromTriplets(hloc.begin(), hloc.end());
  // INFO(H_local);
}

template<typename Tnum, typename Tlabel>
void Hamiltonian<Tnum, Tlabel>::BuildHoppingHamiltonian(
  const std::vector<Basis> &bs, const std::vector< Node<Tnum, Tlabel>* > &lt )
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
  // INFO(H_hop);
}

template<typename Tnum, typename Tlabel>
void Hamiltonian<Tnum, Tlabel>::eigh( RealType &Val, VectorType &Vec,
  const bool FullDiagonalization)
{
  Vec = VectorType::Random(getTotalHilbertSpace());
  size_t max_iter = 200;
  RealType err_tol = 1.0e-9;
  if( !(FullDiagonalization) ){
    bool converge = LanczosEV(H_total, Vec, Val, max_iter, err_tol);
    if ( !(converge) ) INFO("Lanczos is not converged!");
  }
}

template<>
void Hamiltonian<ComplexType, int>::expH( const ComplexType Prefactor,
  ComplexVectorType &Vec, const size_t Kmax )
{
  krylov(H_total, Vec, Prefactor, Kmax);
}

template class Hamiltonian<RealType, int>;
template class Hamiltonian<ComplexType, int>;
