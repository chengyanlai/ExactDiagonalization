#include <cassert>
#include <cmath>//sqrt
#include <map>
#include <tuple>
#include "src/Hamiltonian/Hamiltonian.hpp"

template<typename Tnum>
void Hamiltonian<Tnum>::LocalTerms( const std::vector<Tnum> &Vloc, const std::vector<Tnum> &Uloc, const Basis& bs, std::vector<std::tuple<int, int, Tnum> > &MatElemts ){
  int state_id = 0;
  for ( std::vector<int> b : bs.BStates ){
    int loc = 0;
    Tnum val = (Tnum)0.0;
    for ( int &p : b ){
      val += Vloc.at(loc) * (Tnum)p + (Tnum)0.50e0 * Uloc.at(loc) * (Tnum)(p * (p - 1));
      loc++;
    }
    size_t idx = DetermineTotalIndex( std::vector<size_t>(1, state_id) );
    // hloc.push_back(Triplet(idx, idx, val));
    MatElemts.push_back( std::make_tuple(idx, idx, val) );
    state_id++;
  }
}

template<typename Tnum>
void Hamiltonian<Tnum>::NNHopping( const std::vector< Node<Tnum>* > &lt, const Basis& bs, std::vector<std::tuple<int, int, Tnum> > &MatElemts ){
  int state_id1 = 0;
  for ( std::vector<int> b : bs.BStates ){
    for ( Node<Tnum>* l : lt ) {
      size_t site_i = l->data;
      std::vector< Node<Tnum>* > nn = l->GetNeighbors();
      std::vector< Tnum > nnJ = l->GetJval();
      for (size_t cnt = 0; cnt < l->NumNeighbors(); cnt++) {
        std::vector<int> b_copy = b;
        size_t site_j = nn.at(cnt)->data;
        if ( b_copy.at(site_i) > 0 ){
          Tnum val = (Tnum)(-1.0e0) * nnJ.at(cnt) * sqrt( (Tnum)((b.at(site_j) + 1) * b.at(site_i)) );
          b_copy.at(site_j) = b.at(site_j) + 1;
          b_copy.at(site_i) = b.at(site_i) - 1;
          size_t state_id2 = bs.GetIndexFromTag( BosonBasisTag(b_copy) );
          if ( state_id2 < bs.GetHilbertSpace() ){
            size_t idx1 = DetermineTotalIndex(std::vector<size_t>(1, state_id1));
            size_t idx2 = DetermineTotalIndex(std::vector<size_t>(1, state_id2));
            MatElemts.push_back( std::make_tuple(idx1, idx2, val) );
          }
        }
      }
    }
    state_id1++;
  }
}

template<typename Tnum>
void Hamiltonian<Tnum>::BoseHubbardModel( const std::vector<Basis>& bs, const std::vector< Node<Tnum>* >& lattice, const std::vector<Tnum>& Vloc, const std::vector<Tnum>& Uloc ){
  // Both species live in the same lattice and the same happing amplitudes
  assert( Vloc.size() == Uloc.size() );
  assert( Vloc.size() == bs.at(0).GetL() );
  std::vector<std::tuple<int, int, Tnum> > MatElemts;
  MatElemts.clear();
  /* Local terms contain potential and interaction */
  LocalTerms( Vloc, Uloc, bs.at(0), MatElemts );
  /* For N-N hopping */
  assert( bs.at(0).GetL() == lattice.size() );
  NNHopping( lattice, bs.at(0), MatElemts );
  /* Build H_total */
  BuildTotalHamiltonian( MatElemts );
}

template class Hamiltonian<RealType>;
template class Hamiltonian<ComplexType>;
