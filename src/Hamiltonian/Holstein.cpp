#include <cassert>
#include <cmath>//exp
#include <tuple>
#include "src/Hamiltonian/Hamiltonian.hpp"

template<typename Tnum>
void Hamiltonian<Tnum>::LocalPhonon( const std::vector<Tnum>& Wloc, const Basis& bs, std::vector<std::tuple<int, int, Tnum> > &MatElemts ){
  int state_id = 0;
  for ( std::vector<int> b : bs.BStates ){
    int loc = 0;
    Tnum val = (Tnum)0.0;
    for ( int &p : b ){
      val += Wloc.at(loc) * (Tnum)p;
      loc++;
    }
    MatElemts.push_back( std::make_tuple(state_id, state_id, val) );
    state_id++;
  }
}

template<typename Tnum>
void Hamiltonian<Tnum>::FermionPhononCoupling( const std::vector<Tnum>& Gloc, const Basis& bs, std::vector<std::tuple<int, int, Tnum> > &MatElemts ){
  std::vector< std::vector<int> > States = bs.BStates;
  size_t rid = 0;
  typename std::vector< std::vector<int> >::const_iterator it = States.begin();
  for (; it != States.end(); ++it ){
    std::vector<int> state = *it;
    RealType tg = bs.CreatePhonon(state);
    if ( tg > -1.0e-5 ){
      size_t cid = bs.GetIndexFromTag(tg);
      if ( state == bs.BStates.at(cid) ) MatElemts.push_back( std::make_tuple(rid, cid, -1.0e0 * Gloc.at(0)) );
    }
    state = *it;
    tg = bs.DestroyPhonon(state);
    if ( tg > -1.0e-5 ){
      size_t cid = bs.GetIndexFromTag(tg);
      if ( state == bs.BStates.at(cid) ) MatElemts.push_back( std::make_tuple(rid, cid, -1.0e0 * Gloc.at(0)) );
    }
    rid++;
  }
}

template<typename Tnum>
void Hamiltonian<Tnum>::FermionNNHopping( const RealType& k, const std::vector< Node<Tnum>* > &lt, const Basis& bs, std::vector<std::tuple<int, int, Tnum> > &MatElemts ){
  int L = bs.GetL();
  int state_id1 = 0;
  for ( std::vector<int> b : bs.BStates ){
    for ( Node<Tnum>* l : lt ) {
      size_t site_i = l->data;
      std::vector< Node<Tnum>* > nn = l->GetNeighbors();
      std::vector< Tnum > nnJ = l->GetJval();
      for (size_t cnt = 0; cnt < l->NumNeighbors(); cnt++) {
        size_t site_j = nn.at(cnt)->data;
        std::vector<int> state1 = b;
        RealType tg1;
        size_t id1 = state_id1;
        if ( site_i ){
          tg1 = bs.FermionJumpRight(state1, site_i);
          id1 = bs.GetIndexFromTag(tg1);
        }
        std::vector<int> state2 = b;
        RealType tg2;
        size_t id2 = state_id1;
        if ( site_j ){
          tg2 = bs.FermionJumpRight(state2, site_j);
          id2 = bs.GetIndexFromTag(tg2);
        }
        if ( state1 == bs.BStates.at(id1) && state2 == bs.BStates.at(id2) ){
          RealType r = RealType(site_i) - RealType(site_j);
          if ( site_i == 0 && site_j == L - 1 ) r = -1.0;
          else if ( site_i == L - 1 && site_j == 0 ) r = 1.0;
          MatElemts.push_back( std::make_tuple(id1, id2, -1.0e0 * nnJ.at(cnt) * exp( ComplexType(0.0, 1.0) * k * PI * r ) ) );// Only Nearest-Neighbor !!
// std::cout << id1 << " " << id2 << " " << site_i << " " << site_j << " " << r << " " << -1.0e0 * nnJ.at(cnt) * exp( ComplexType(0.0, 1.0) * k * PI * r ) << std::endl;
        }
      }
    }
    state_id1++;
  }
}

template<typename Tnum>
void Hamiltonian<Tnum>::HolsteinModel( const std::vector<Basis>& bs, const RealType& k, const std::vector< Node<Tnum>* >& lattice, const std::vector<Tnum>& Wloc, const std::vector<Tnum>& Gloc ){
  // Both species live in the same lattice and the same happing amplitudes
  assert( Wloc.size() == bs.at(0).GetL() );
  assert( Gloc.size() == bs.at(0).GetL() );
  std::vector<std::tuple<int, int, Tnum> > MatElemts;
  MatElemts.clear();
  /* Build H_total */
  LocalPhonon( Wloc, bs.at(0), MatElemts );
  FermionPhononCoupling( Gloc, bs.at(0), MatElemts );
  FermionNNHopping( k, lattice, bs.at(0), MatElemts );
  BuildTotalHamiltonian( MatElemts );
}

template class Hamiltonian<ComplexType>;
