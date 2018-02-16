#include <cassert>
#include <cmath>//exp
#include <tuple>
#include "src/Hamiltonian/Holstein/Holstein.hpp"

template<typename Tnum>
void Holstein<Tnum>::PhononLocal( const std::vector<Tnum>& Wloc, const Basis& bs ){
  std::vector<std::tuple<int, int, Tnum> > MatElemts;
  MatElemts.clear();
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
  H_Phonon = this->BuildSparseHamiltonian( MatElemts );
}

template<typename Tnum>
void Holstein<Tnum>::FermionPhononCoupling( const std::vector<Tnum>& Gloc, const Basis& bs){
  std::vector<std::tuple<int, int, Tnum> > MatElemts;
  std::vector< std::vector<int> > States = bs.BStates;
  size_t rid = 0;
  typename std::vector< std::vector<int> >::const_iterator it = States.begin();
  for (; it != States.end(); ++it ){
    std::vector<int> state = *it;
    RealType tg = bs.CreatePhonon(state);
    RealType Npf = state.at(0);
    if ( tg > -1.0e-5 ){
      size_t cid = bs.GetIndexFromTag(tg);
      if ( state == bs.BStates.at(cid) ) MatElemts.push_back( std::make_tuple(cid, rid, -1.0e0 * Gloc.at(0) * sqrt(Npf) ) );
    }
    state = *it;
    if ( state.at(0) ){
      RealType Npi = state.at(0);
      tg = bs.DestroyPhonon(state);
      size_t cid = bs.GetIndexFromTag(tg);
      if ( state == bs.BStates.at(cid) ) MatElemts.push_back( std::make_tuple(cid, rid, -1.0e0 * Gloc.at(0) * sqrt(Npi) ) );
    }
    rid++;
  }
  H_Couple = this->BuildSparseHamiltonian( MatElemts );
}

template<typename Tnum>
void Holstein<Tnum>::FermionNNHoppingInfinite( const RealType& k, const RealType& J, const Basis& bs, const RealType& Phi){
  std::vector<std::tuple<int, int, Tnum> > MatElemts;
  int id1 = 0;
  for ( std::vector<int> b : bs.BStates ){
    std::vector<int> state = b;
    RealType tg2;
    size_t id2;
    // Jump left
    tg2 = bs.FermionJumpLeft(state, 1);
    id2 = bs.GetIndexFromTag(tg2);
    if ( state == bs.BStates.at(id2) ){
      MatElemts.push_back( std::make_tuple(id1, id2, -1.0e0 * J * exp( ComplexType(0.0, 1.0) * (k + Phi) * PI ) ) );
    }
    state = b;
    tg2 = bs.FermionJumpRight(state, 1);
    id2 = bs.GetIndexFromTag(tg2);
    if ( state == bs.BStates.at(id2) ){
      MatElemts.push_back( std::make_tuple(id1, id2, -1.0e0 * J * exp( ComplexType(0.0,-1.0) * (k + Phi) * PI ) ) );
    }
    id1++;
  }
  H_Kinetic = this->BuildSparseHamiltonian( MatElemts );
}

template<typename Tnum>
void Holstein<Tnum>::HolsteinModel( const std::vector<Basis>& bs, const RealType& k, const RealType& J, const std::vector<Tnum>& Wloc, const std::vector<Tnum>& Gloc ){
  // Both species live in the same lattice and the same happing amplitudes
  assert( Wloc.size() == bs.at(0).GetL() );
  assert( Gloc.size() == bs.at(0).GetL() );
  PhononLocal( Wloc, bs.at(0) );
  FermionPhononCoupling( Gloc, bs.at(0) );
  FermionNNHoppingInfinite( k, J, bs.at(0) );
  /* Build H_total */
  this->H_total = H_Phonon + H_Kinetic + H_Couple;
}



template<typename Tnum>
void Holstein<Tnum>::PhononLocalK( const int& KPoints, const RealType& W, const Basis& bs ){
  std::vector<std::tuple<int, int, Tnum> > MatElemts;
  MatElemts.clear();
  int state_id = 0;
  for ( std::vector<int> b : bs.BStates ){
    int Np = std::accumulate(b.begin(), b.end(), 0);
    Tnum val = (Tnum)W * (Tnum)Np;
    for ( size_t j = 0; j< KPoints; j++ ){
      size_t idx = this->DetermineTotalIndex( vec<int>(j, state_id) );
      MatElemts.push_back( std::make_tuple(idx, idx, val) );
    }
    state_id++;
  }
  H_Phonon = this->BuildSparseHamiltonian( MatElemts );
}

template<typename Tnum>
void Holstein<Tnum>::FermionNNHoppingK( const int& KPoints, const RealType& J, const Basis& bs, const RealType& Phi){
  std::vector<std::tuple<int, int, Tnum> > MatElemts;
  for ( size_t j = 0; j < KPoints; j++ ){
    RealType Kf = RealType(j) * PI / RealType(KPoints);
    Tnum val = -2.0e0 * J * cos( (Kf + Phi) * PI );
    for ( size_t cnt = 0; cnt < bs.GetHilbertSpace(); cnt++ ){
      size_t idx = this->DetermineTotalIndex( vec<int>(j, cnt) );
      MatElemts.push_back( std::make_tuple(idx, idx, val) );
    }
  }
  H_Kinetic = this->BuildSparseHamiltonian( MatElemts );
}

template<typename Tnum>
void Holstein<Tnum>::FermionPhononCouplingK( const std::vector<int>& KPoints, const RealType& G, const Basis& bs){
  std::vector<std::tuple<int, int, Tnum> > MatElemts;
  std::vector< std::vector<int> > States = bs.BStates;
  for ( size_t jF = 1; jF < KPoints.size(); jF++ ){// there are no k=0 phonon mode
    int KF = KPoints.at(jF);
    size_t rid = 0;
    typename std::vector< std::vector<int> >::const_iterator it = States.begin();
    for (; it != States.end(); ++it ){
      std::vector<int> state = *it;
      RealType tg = bs.CreatePhonon(state);
      RealType Npf = state.at(0);
      if ( tg > -1.0e-5 ){
        size_t cid = bs.GetIndexFromTag(tg);
        if ( state == bs.BStates.at(cid) ) MatElemts.push_back( std::make_tuple(cid, rid, -1.0e0 * Gloc.at(0) * sqrt(Npf) ) );
      }
      state = *it;
      if ( state.at(0) ){
        RealType Npi = state.at(0);
        tg = bs.DestroyPhonon(state);
        size_t cid = bs.GetIndexFromTag(tg);
        if ( state == bs.BStates.at(cid) ) MatElemts.push_back( std::make_tuple(cid, rid, -1.0e0 * Gloc.at(0) * sqrt(Npi) ) );
      }
      rid++;
    }
  }
  H_Couple = this->BuildSparseHamiltonian( MatElemts );
}

template<typename Tnum>
void Holstein<Tnum>::HolsteinModelK( const std::vector<int>& KPoints, const std::vector<Basis>& bs, const RealType& W, const RealType& G, const RealType& J ){
  PhononLocalK( KPoints.size(), W, bs.at(0) );
  FermionNNHoppingK( KPoints.size(), J, bs.at(0) );
  FermionPhononCouplingK( KPoints, G, bs.at(0) );
  /* Build H_total */
  this->H_total = H_Phonon + H_Kinetic + H_Couple;
}

template class Holstein<ComplexType>;
