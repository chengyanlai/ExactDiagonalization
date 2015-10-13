#include "src/Hamiltonian/Hamiltonian.hpp"
#include "src/Lanczos/lanczos.hpp"

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
void Hamiltonian<Tnum, Tlabel>::BuildIntraLocalHamiltonian(
  const std::vector< std::vector<Tnum> > &Vloc,
  const std::vector< std::vector<Tnum> > &Uloc,
  const std::vector<Basis> &bs )
{
  std::vector<Triplet> hloc;
  hloc.clear();
  assert( Vloc.size() == Uloc.size() );
  assert( Vloc.size() == bs.size() );
  int cnt = 0;
  for ( auto &b : bs ){
    assert( b.getL() == Vloc.at(cnt).size() );
    assert( b.getL() == Uloc.at(cnt).size() );
    if( !(b.getType()) ){//boson
      BosonIntraLocalPart( cnt, Vloc.at(cnt), Uloc.at(cnt), b, hloc );
    }
    cnt++;
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
void Hamiltonian<Tnum, Tlabel>::BuildIntraHoppingHamiltonian(
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
      BosonIntraHoppingPart( cnt, lt, b, hhop );
    }
    cnt++;
  }
  H_hop.setFromTriplets(hhop.begin(), hhop.end());
  // INFO(H_hop);
}

template<typename Tnum, typename Tlabel>
void Hamiltonian<Tnum, Tlabel>::eigh( const bool FullDiagonalization )
{
  RealType Val;
  size_t max_iter = 200;
  RealType err_tol = 1.0e-9;
  VectorType Vec = RealVectorType::Random(getTotalHilbertSpace());
  if( !(FullDiagonalization) ){
    bool converge = LanczosEV(getTotalHilbertSpace(), H_total,
      Vec, Val, max_iter, err_tol);
    if ( !(converge) ) INFO("Lanczos is not converged!");
  }
  INFO(Val);
  // INFO(Vec);
}

template class Hamiltonian<RealType, int>;
template class Hamiltonian<RealType, char>;
template class Hamiltonian<ComplexType, int>;
template class Hamiltonian<ComplexType, char>;
