#include <cassert>
#include <cmath>//sqrt
#include <map>
#include <tuple>
#include "src/Hamiltonian/BHM/BoseHubbard.hpp"

template<typename Tnum>
void BHM<Tnum>::LocalTerms( const std::vector<Tnum> &Vloc, const std::vector<Tnum> &Uloc, const Basis& bs){
  std::vector<std::tuple<int, int, Tnum> > MatElemtsV, MatElemtsU;
  MatElemtsV.clear();
  MatElemtsU.clear();
  int state_id = 0;
  for ( std::vector<int> b : bs.BStates ){
    int loc = 0;
    Tnum Vval = (Tnum)0.0;
    Tnum Uval = (Tnum)0.0;
    for ( int &p : b ){
      Vval += Vloc.at(loc) * (Tnum)p;
      Uval += (Tnum)0.50e0 * Uloc.at(loc) * (Tnum)(p * (p - 1));
      loc++;
    }
    size_t idx = this->DetermineTotalIndex( std::vector<size_t>(1, state_id) );
    // hloc.push_back(Triplet(idx, idx, val));
    MatElemtsV.push_back( std::make_tuple(idx, idx, Vval) );
    MatElemtsU.push_back( std::make_tuple(idx, idx, Uval) );
    state_id++;
  }
  H_V = BuildSparseHamiltonian( this->GetTotalHilbertSpace(), MatElemtsV );
  H_U = BuildSparseHamiltonian( this->GetTotalHilbertSpace(), MatElemtsU );
}

template<typename Tnum>
void BHM<Tnum>::NNHopping( const std::vector< Node<Tnum>* > &lt, const Basis& bs){
  std::vector<std::tuple<int, int, Tnum> > MatElemts;
  MatElemts.clear();
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
            size_t idx1 = this->DetermineTotalIndex(std::vector<size_t>(1, state_id1));
            size_t idx2 = this->DetermineTotalIndex(std::vector<size_t>(1, state_id2));
            MatElemts.push_back( std::make_tuple(idx1, idx2, val) );
          }
        }
      }
    }
    state_id1++;
  }
  H_J = BuildSparseHamiltonian( this->GetTotalHilbertSpace(), MatElemts );
}

template<typename Tnum>
void BHM<Tnum>::BoseHubbardModel( const std::vector<Basis>& bs, const std::vector< Node<Tnum>* >& lattice, const std::vector<Tnum>& Vloc, const std::vector<Tnum>& Uloc ){
  // Both species live in the same lattice and the same happing amplitudes
  assert( Vloc.size() == Uloc.size() );
  assert( Vloc.size() == bs.at(0).GetL() );
  /* Local terms contain potential and interaction */
  LocalTerms( Vloc, Uloc, bs.at(0) );
  /* For N-N hopping */
  assert( bs.at(0).GetL() == lattice.size() );
  NNHopping( lattice, bs.at(0) );
  /* Build H_total */
  this->H_total = H_J = H_V + H_U;
}

template class BHM<RealType>;
template class BHM<ComplexType>;
