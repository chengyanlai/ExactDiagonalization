#include <cassert>
#include <cmath>//exp
#include <tuple>
#include <numeric>
#include "src/Hamiltonian/Holstein/Holstein.hpp"

//* The following function designed to work with Holstein model in limited functional space - thermodynamic limit.
template<typename Tnum>
void Holstein<Tnum>::PhononLocal( const RealType& Wloc, const Basis& bs ){
  H_Phonon.zeros();
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
  H_Couple.zeros();
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
      if ( state == bs.BStates.at(cid) ) MatElemts.push_back( std::make_tuple(rid, cid, -1.0e0 * Gloc * sqrt(Npf) ) );
    }
    state = *it;
    if ( state.at(0) ){
      RealType Npi = state.at(0);
      tg = bs.DestroyPhonon(state);
      size_t cid = bs.GetIndexFromTag(tg);
      if ( state == bs.BStates.at(cid) ) MatElemts.push_back( std::make_tuple(rid, cid, -1.0e0 * Gloc * sqrt(Npi) ) );
    }
    rid++;
  }
  H_Couple = BuildSparseHamiltonian( this->GetTotalHilbertSpace(), MatElemts );
}

template<>
void Holstein<ComplexType>::FermionNNHopping( const RealType& k, const RealType& J, const Basis& bs, const RealType& Phi){
  H_Kinetic.zeros();
  std::vector<std::tuple<int, int, ComplexType> > MatElemts;
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

template<>
void Holstein<ComplexType>::HolsteinModel( const std::vector<Basis>& bs, const RealType& k, const RealType& J, const RealType& Wloc, const RealType& Gloc ){
  PhononLocal( Wloc, bs.at(0) );
  FermionPhononCoupling( Gloc, bs.at(0) );
  FermionNNHopping( k, J, bs.at(0) );
  //* Build H_total
  this->H_total = H_Phonon + H_Kinetic + H_Couple;
}


//* The following function designed to work with Holstein model in real space.
template<typename Tnum>
void Holstein<Tnum>::PhononR( const int& RPoints, const RealType& W, const Basis& bs ){
  H_Phonon.zeros();
  std::vector<std::tuple<int, int, Tnum> > MatElemts;
  MatElemts.clear();
  int state_id = 0;
  for ( std::vector<int> b : bs.BStates ){
    int Np = std::accumulate(b.begin(), b.end(), 0);
    Tnum val = (Tnum)W * (Tnum)Np;
    for ( size_t j = 0; j < RPoints; j++ ){// For fermion positions
      size_t idx = this->DetermineTotalIndex( vec<size_t>(j, state_id) );
      MatElemts.push_back( std::make_tuple(idx, idx, val) );
    }
    state_id++;
  }
  H_Phonon = BuildSparseHamiltonian( this->GetTotalHilbertSpace(), MatElemts );
}

template<typename Tnum>
void Holstein<Tnum>::FermionR( const int& RPoints, const RealType& J, const Basis& bs ){
  H_Kinetic.zeros();
  std::vector<std::tuple<int, int, Tnum> > MatElemts;
  for ( size_t j = 0; j < RPoints; j++ ){
    size_t jp = (j == RPoints - 1) ? 0 : j + 1;
    for ( size_t cnt = 0; cnt < bs.GetHilbertSpace(); cnt++ ){
      size_t idx1 = this->DetermineTotalIndex( vec<size_t>(j, cnt) );
      size_t idx2 = this->DetermineTotalIndex( vec<size_t>(jp, cnt) );
      Tnum val = -1.0 * (Tnum)J;
      MatElemts.push_back( std::make_tuple(idx1, idx2, val) );
      MatElemts.push_back( std::make_tuple(idx2, idx1, val) );
    }
  }
  H_Kinetic = BuildSparseHamiltonian( this->GetTotalHilbertSpace(), MatElemts );
}

template<typename Tnum>
void Holstein<Tnum>::FermionPhononR( const int& RPoints, const RealType& G, const Basis& bs ){
  H_Couple.zeros();
  std::vector<std::tuple<int, int, Tnum> > MatElemts;
  std::vector< std::vector<int> > States = bs.BStates;
  size_t cid = 0;
  typename std::vector< std::vector<int> >::const_iterator it = States.begin();
  for (; it != States.end(); ++it ){
    for ( int Floc = 0; Floc< RPoints; Floc++ ){
      std::vector<int> state = *it;
      RealType tg = bs.CreatePhonon(state, Floc);
      RealType Npf = state.at(Floc);
      if ( tg > -1.0e-5 ){
        size_t rid = bs.GetIndexFromTag(tg);
        if ( state == bs.BStates.at(cid) ){
          size_t ridx = this->DetermineTotalIndex( vec<size_t>(Floc, rid) );
          size_t cidx = this->DetermineTotalIndex( vec<size_t>(Floc, cid) );
          MatElemts.push_back( std::make_tuple(ridx, cidx, -1.0e0 * G * sqrt(Npf) ) );
        }
      }
      state = *it;
      if ( state.at(Floc) ){
        RealType Npi = state.at(Floc);
        tg = bs.DestroyPhonon(state, Floc);
        if ( tg > -1.0e-5 ){
          size_t rid = bs.GetIndexFromTag(tg);
          if ( state == bs.BStates.at(cid) ) {
            size_t ridx = this->DetermineTotalIndex( vec<size_t>(Floc, rid) );
            size_t cidx = this->DetermineTotalIndex( vec<size_t>(Floc, cid) );
            MatElemts.push_back( std::make_tuple(ridx, cidx, -1.0e0 * G * sqrt(Npi) ) );
          }
        }
      }
    }
    cid++;
  }
  H_Couple = BuildSparseHamiltonian( this->GetTotalHilbertSpace(), MatElemts );
}

template<typename Tnum>
void Holstein<Tnum>::HolsteinModelR( const int& RPoints, const std::vector<Basis>& bs, const RealType& W, const RealType& G, const RealType& J ){
  PhononR( RPoints, W, bs.at(0) );
  FermionR( RPoints, J, bs.at(0) );
  FermionPhononR( RPoints, G, bs.at(0) );
  //* Build H_total */
  this->H_total = H_Phonon + H_Kinetic + H_Couple;
}


//* The following function designed to work with Holstein model in momentum space.
template<typename Tnum>
void Holstein<Tnum>::PhononK( const Basis& bs, const RealType& W ){
  H_Phonon.zeros();
  std::vector<std::tuple<int, int, Tnum> > MatElemts;
  MatElemts.clear();
  int state_id = 0;
  for ( std::vector<int> b : bs.BStates ){
    int Np = std::accumulate(b.begin()+1, b.end(), 0);//! First element is the fermion momentum
    Tnum val = (Tnum)W * (Tnum)Np;
    MatElemts.push_back( std::make_tuple(state_id, state_id, val) );
    state_id++;
  }
  H_Phonon = BuildSparseHamiltonian( this->GetTotalHilbertSpace(), MatElemts );
}

template<typename Tnum>
void Holstein<Tnum>::FermionK( const std::vector<int>& KPoints, const Basis& bs, const RealType& J, const RealType& Phi){
  H_Kinetic.zeros();
  std::vector<std::tuple<int, int, Tnum> > MatElemts;
  MatElemts.clear();
  int state_id = 0;
  for ( std::vector<int> b : bs.BStates ){
    RealType Kf = KPoints.at(b.at(0)) * PI / RealType(KPoints.size()/2);
    Tnum val = -2.0e0 * J * cos( Kf + Phi * PI );
    MatElemts.push_back( std::make_tuple(state_id, state_id, val) );
    state_id++;
  }
  H_Kinetic = BuildSparseHamiltonian( this->GetTotalHilbertSpace(), MatElemts );
}

template<typename Tnum>
void Holstein<Tnum>::FermionPhononK( const std::vector<int>& KPoints, const Basis& bs, const RealType& G, const int& Nmax, const int WithoutK0Phonon){
  H_Couple.zeros();
  std::vector<std::tuple<int, int, Tnum> > MatElemts;
  MatElemts.clear();
  int L = bs.GetL();
  size_t cid = 0;
  for ( std::vector<int> b : bs.BStates ){
    // //std::cout << "Original" << std::flush;
    // //PrintVector(b, 3, " ");
    int Kfi = b.at(0);
    for ( int Ki = 0; Ki < KPoints.size(); Ki++ ){
      if ( Kfi == Ki ) continue;//! No momentum changes
      // Create -q phonon
      int DeltaK = KPoints.at(Kfi) - KPoints.at(Ki);//* This is -q
      int DeltaKi = DeltaKIndex( DeltaK, KPoints );
      int Idx = DeltaKi+1-WithoutK0Phonon;//! The index in the basis vector
      std::vector<int> state = b;
      state.at(0) = Ki;//* Fermion from k to k+q
      state.at(Idx) += 1;
      if ( std::accumulate(state.begin()+1, state.end(), 0) <= Nmax ){//* state exist in basis
        // //std::cout << "\tCreate -q  " << std::flush;
        // //PrintVector(state, 3, " ");
        RealType tg = BosonBasisTag(state);
        if ( !(bs.DummyCheckState(tg, state)) ) throw std::runtime_error("State is not match!");//! Dummy check
        size_t rid = bs.GetIndexFromTag(tg);
        // //std::cout << "\tFrom Tag -q" << std::flush;
        // //PrintVector(bs.BStates.at(rid), 3, " ");
        RealType Np = RealType(state.at(Idx)) / RealType(L);
        MatElemts.push_back( std::make_tuple(rid, cid, -1.0e0 * (Tnum)G * sqrt(Np) ) );
      }
      // Destory q phonon
      DeltaK = KPoints.at(Ki) - KPoints.at(Kfi);//* This is +q
      DeltaKi = DeltaKIndex( DeltaK, KPoints );
      Idx = DeltaKi+1-WithoutK0Phonon;//! The index in the basis vector
      state = b;
      if ( state.at(Idx) ){
        state.at(0) = Ki;//* Fermion from k to k+q
        RealType Np = RealType(state.at(Idx)) / RealType(L);
        state.at(Idx) -= 1;
        // //std::cout << "\tDestory +q " << std::flush;
        // //PrintVector(state, 3, " ");
        RealType tg = BosonBasisTag(state);
        if ( !(bs.DummyCheckState(tg, state)) ) throw std::runtime_error("State is not match!");//! Dummy check
        size_t rid = bs.GetIndexFromTag(tg);
        // //std::cout << "\tFrom Tag +q" << std::flush;
        // //PrintVector(bs.BStates.at(rid), 3, " ");
        MatElemts.push_back( std::make_tuple(rid, cid, -1.0e0 * (Tnum)G * sqrt(Np) ) );
      }
    }
    cid++;
  }
  H_Couple = BuildSparseHamiltonian( this->GetTotalHilbertSpace(), MatElemts );
}

template<typename Tnum>
void Holstein<Tnum>::HolsteinModelK( const std::vector<int>& KPoints, const std::vector<Basis>& bs, const RealType& W, const RealType& G, const RealType& J, const int& Nmax, const int WithoutK0Phonon ){
  PhononK( bs.at(0), W );
  FermionK( KPoints, bs.at(0), J );
  FermionPhononK( KPoints, bs.at(0), G, Nmax, WithoutK0Phonon );
  //* Build H_total
  this->H_total = H_Phonon + H_Kinetic + H_Couple;
}

template class Holstein<RealType>;
template class Holstein<ComplexType>;
