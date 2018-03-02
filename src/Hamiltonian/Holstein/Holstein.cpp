#include <cassert>
#include <cmath>//exp
#include <tuple>
#include <numeric>
#include "src/Hamiltonian/Holstein/Holstein.hpp"

template<typename Tnum>
void Holstein<Tnum>::PhononLocal( const RealType& Wloc, const Basis& bs ){
  std::vector<std::tuple<int, int, Tnum> > MatElemts;
  MatElemts.clear();
  int state_id = 0;
  for ( std::vector<int> b : bs.BStates ){
    int loc = 0;
    Tnum val = (Tnum)Wloc * (Tnum)std::accumulate(b.begin(), b.end(), 0);
    MatElemts.push_back( std::make_tuple(state_id, state_id, val) );
    state_id++;
  }
  H_Phonon = BuildSparseHamiltonian( this->GetTotalHilbertSpace(), MatElemts );
}

template<typename Tnum>
void Holstein<Tnum>::FermionPhononCoupling( const RealType& Gloc, const Basis& bs){
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
      if ( state == bs.BStates.at(cid) ) MatElemts.push_back( std::make_tuple(cid, rid, -1.0e0 * Gloc * sqrt(Npf) ) );
    }
    state = *it;
    if ( state.at(0) ){
      RealType Npi = state.at(0);
      tg = bs.DestroyPhonon(state);
      size_t cid = bs.GetIndexFromTag(tg);
      if ( state == bs.BStates.at(cid) ) MatElemts.push_back( std::make_tuple(cid, rid, -1.0e0 * Gloc * sqrt(Npi) ) );
    }
    rid++;
  }
  H_Couple = BuildSparseHamiltonian( this->GetTotalHilbertSpace(), MatElemts );
}

template<typename Tnum>
void Holstein<Tnum>::FermionNNHopping( const RealType& k, const RealType& J, const Basis& bs, const RealType& Phi){
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
  H_Kinetic = BuildSparseHamiltonian( this->GetTotalHilbertSpace(), MatElemts );
}

template<typename Tnum>
void Holstein<Tnum>::HolsteinModel( const std::vector<Basis>& bs, const RealType& k, const RealType& J, const RealType& Wloc, const RealType& Gloc ){
  // Both species live in the same lattice and the same happing amplitudes
  PhononLocal( Wloc, bs.at(0) );
  FermionPhononCoupling( Gloc, bs.at(0) );
  FermionNNHopping( k, J, bs.at(0) );
  /* Build H_total */
  this->H_total = H_Phonon + H_Kinetic + H_Couple;
}

/* The following function designed to work with Holstein model in momentum space. */

template<typename Tnum>
void Holstein<Tnum>::PhononK( const int& KPoints, const RealType& W, const Basis& bs ){
  std::vector<std::tuple<int, int, Tnum> > MatElemts;
  MatElemts.clear();
  int state_id = 0;
  for ( std::vector<int> b : bs.BStates ){
    int Np = std::accumulate(b.begin(), b.end(), 0);
    Tnum val = (Tnum)W * (Tnum)Np;
    for ( size_t j = 0; j< KPoints; j++ ){
      size_t idx = this->DetermineTotalIndex( vec<size_t>(j, state_id) );
      MatElemts.push_back( std::make_tuple(idx, idx, val) );
    }
    state_id++;
  }
  H_Phonon = BuildSparseHamiltonian( this->GetTotalHilbertSpace(), MatElemts );
}

template<typename Tnum>
void Holstein<Tnum>::FermionK( const int& KPoints, const RealType& J, const Basis& bs, const RealType& Phi){
  std::vector<std::tuple<int, int, Tnum> > MatElemts;
  for ( size_t j = 0; j < KPoints; j++ ){
    RealType Kf = RealType(j) * PI / RealType(KPoints);
    Tnum val = -2.0e0 * J * cos( (Kf + Phi) * PI );
    for ( size_t cnt = 0; cnt < bs.GetHilbertSpace(); cnt++ ){
      size_t idx = this->DetermineTotalIndex( vec<size_t>(j, cnt) );
      MatElemts.push_back( std::make_tuple(idx, idx, val) );
    }
  }
  H_Kinetic = BuildSparseHamiltonian( this->GetTotalHilbertSpace(), MatElemts );
}

template<typename Tnum>
void Holstein<Tnum>::FermionPhononK( const std::vector<int>& KPoints, const RealType& G, const Basis& bs, const int& Nmax){
  std::vector<std::tuple<int, int, Tnum> > MatElemts;
  int L = KPoints.size();
  for ( int Ki = 0; Ki < L; Ki++ ){
    int Kn = KPoints.at(Ki);
    for ( int Qi = 0; Qi < L-1; Qi++ ){// from 1 due to the k=0 phonon mode
      int Qn = KPoints.at( Qi + 1 );// +1 due to the k=0 phonon mode
      int Pn = Kn + Qn;
      FirstBZ(L, Pn);
      int Pi = KnToIndex(Pn);// this is index for KPoints not the phonon basis!!
      std::vector< std::vector<int> > States = bs.BStates;
      size_t rPid = 0;
      typename std::vector< std::vector<int> >::const_iterator it = States.begin();
      for (; it != States.end(); ++it ){
        // Create -q phonon
        // The actual index for phonon basis is Qi
        std::vector<int> state = *it;
        int cQi;
        if ( Qi == L-2 ) cQi = L-2;// Handle pi <-> -pi
        else if ( Qi % 2 == 0 ) cQi = Qi + 1;// plus to minus
        else cQi = Qi - 1;// minus to plus
        state.at(cQi) += 1;
        if ( std::accumulate(state.begin(), state.end(), 0) <= Nmax ){// state exist in basis
          RealType tg = BosonBasisTag(state);
          size_t cPid = bs.GetIndexFromTag(tg);
          size_t rid = this->DetermineTotalIndex( vec<size_t>(Ki, rPid) );
          size_t cid = this->DetermineTotalIndex( vec<size_t>(Pi, cPid) );
          RealType Np = state.at(cQi);
          MatElemts.push_back( std::make_tuple(rid, cid, -1.0e0 * (Tnum)G * sqrt(Np) ) );
        }
        // Destory q phonon
        state = *it;
        int dQi= Qi;
        if ( state.at(dQi) ){
          RealType Np = state.at(dQi);
          state.at(dQi) -= 1;
          RealType tg = BosonBasisTag(state);
          size_t cPid = bs.GetIndexFromTag(tg);
          size_t rid = this->DetermineTotalIndex( vec<size_t>(Ki, rPid) );
          size_t cid = this->DetermineTotalIndex( vec<size_t>(Pi, cPid) );
          MatElemts.push_back( std::make_tuple(rid, cid, -1.0e0 * (Tnum)G * sqrt(Np) ) );
        }
        rPid++;
      }
    }
  }
  H_Couple = BuildSparseHamiltonian( this->GetTotalHilbertSpace(), MatElemts );
}

template<typename Tnum>
void Holstein<Tnum>::HolsteinModelK( const std::vector<int>& KPoints, const std::vector<Basis>& bs, const RealType& W, const RealType& G, const RealType& J, const int& Nmax ){
  PhononK( KPoints.size(), W, bs.at(0) );
  FermionK( KPoints.size(), J, bs.at(0) );
  FermionPhononK( KPoints, G, bs.at(0), Nmax );
  /* Build H_total */
  this->H_total = H_Phonon + H_Kinetic + H_Couple;
}

template class Holstein<ComplexType>;