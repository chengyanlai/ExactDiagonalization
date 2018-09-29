#include <stddef.h>
#include <cmath>
#include <numeric>//* std::accumulate
#include "src/EDType.hpp"
#include "src/Basis/Basis.hpp"

//* DEBUG use */
static void PrintBosonBasis( const std::vector<int> state ){
  typename std::vector<int>::const_iterator it = state.begin();
  for (;it != state.end(); ++it ){
    std::cout << *it << ", " << std::flush;
  }
  std::cout << std::endl;
};

RealType Basis::CreatePhonon( std::vector<int>& state, const int site )const{
  if ( std::accumulate(state.begin(), state.end(), 0) < N ){
    state.at(site) += 1;
    return BosonBasisTag(state);
  }else{
    return -1;
  }
}

RealType Basis::DestroyPhonon( std::vector<int>& state, const int site )const{
  if ( state.at(site) ){
    state.at(site) -= 1;
    return BosonBasisTag(state);
  }else{
    return -1;
  }
}

//! This is specifically for Limited functional space which applies translational invariant
RealType Basis::FermionJumpRight( std::vector<int>& state, const int NumJumps )const{
  int L = state.size();
  std::vector<int> ns;
  ns.clear();
  ns.insert(ns.begin(), state.begin()+NumJumps, state.end() );
  ns.insert(ns.end(), state.begin(), state.begin()+NumJumps );
  state = ns;
  assert( state.size() == L );
  return BosonBasisTag(state);
}

//! This is specifically for Limited functional space which applies translational invariant
RealType Basis::FermionJumpLeft( std::vector<int>& state, const int NumJumps )const{
  int L = state.size();
  std::vector<int> ns;
  ns.clear();
  ns.insert(ns.begin(), state.end()-NumJumps, state.end() );
  ns.insert(ns.end(), state.begin(), state.end()-NumJumps);
  state = ns;
  assert( state.size() == L );
  return BosonBasisTag(state);
}

void Basis::DummyCheckState(){
  for (size_t i1 = 0; i1 < BStates.size(); i1++){
    std::vector<int> st = BStates.at(i1);
    for (size_t cp = 0; cp < BStates.size(); cp ++){
      if ( st == BStates.at(cp) && i1 != cp ){
        std::cout << i1 <<  " " << BTags.at(i1) << " " << cp << " " << BTags.at(cp) << std::endl;
        PrintBosonBasis(st);
        PrintBosonBasis(BStates.at(cp));
        std::cin.get();
      }
    }
  }
}

//! This is specifically for Limited functional space which applies translational invariant
std::vector<std::vector<int> > Basis::ApplyOffdiagonal( const std::vector<std::vector<int> >& InputStates ){
  //* States has phonon configurations with fermion at site-0
  //* Btags need to be sorted.
  std::vector<std::vector<int> > NewStates;
  NewStates.clear();
  typename std::vector<std::vector<int> >::const_iterator it = InputStates.begin();
  for (; it != InputStates.end(); ++it ){
    //* Create Phonon
    std::vector<int> State = *it;
    RealType tg = CreatePhonon(State);
    assert( State.size() == L );
    size_t Idx = this->GetIndexFromTag(tg);
    if ( !(State == BStates.at(Idx)) ){
      NewStates.push_back(State);
      if ( tg > BTags.at(Idx) ){
        BTags.insert( BTags.begin() + Idx + 1, tg );
        BStates.insert( BStates.begin() + Idx + 1, State );
      }else{
        BTags.insert( BTags.begin() + Idx, tg );
        BStates.insert( BStates.begin() + Idx, State );
      }
    }
    //* Fermion jump right
    State = *it;
    tg = FermionJumpRight( State );
    assert( State.size() == L );
    Idx = this->GetIndexFromTag(tg);
    if ( !(State == BStates.at(Idx)) ){
      NewStates.push_back(State);
      if ( tg > BTags.at(Idx) ){
        BTags.insert( BTags.begin() + Idx + 1, tg );
        BStates.insert( BStates.begin() + Idx + 1, State );
      }else{
        BTags.insert( BTags.begin() + Idx, tg );
        BStates.insert( BStates.begin() + Idx, State );
      }
    }
    //* Fermiom jump left
    State = *it;
    tg = FermionJumpLeft( State );
    assert( State.size() == L );
    Idx = this->GetIndexFromTag(tg);
    if ( !(State == BStates.at(Idx)) ){
      NewStates.push_back(State);
      if ( tg > BTags.at(Idx) ){
        BTags.insert( BTags.begin() + Idx + 1, tg );
        BStates.insert( BStates.begin() + Idx + 1, State );
      }else{
        BTags.insert( BTags.begin() + Idx, tg );
        BStates.insert( BStates.begin() + Idx, State );
      }
    }
  }
  return NewStates;
}

//* Limited function space with one fermion at site-0 */
void Basis::PhononLFS(){
  HaveU1 = false;
  assert( !(isFermion) );
  BStates.clear();
  BTags.clear();
  std::vector<int> State(L, 0);
  std::vector<std::vector<int> > NewStates;
  NewStates.push_back( State );
  BStates.push_back( State );// add this vacuum state
  BTags.push_back( BosonBasisTag(State) );// 0
  for ( int cnt = 0; cnt < N; cnt++ ){
    NewStates = ApplyOffdiagonal( NewStates );
    std::cout << "N_h = " << cnt << " finished, and number of basis = " << BTags.size() << std::endl;
  }
  // DummyCheckState();
  assert( BStates.size() == BTags.size() );
}

//* Holstein phonon in real space
void Basis::PhononR(){
  HaveU1 = false;
  assert( !(isFermion) );
  std::vector< std::vector<int> > work;//* This is final basis
  int IdxK = 0, NWork = 0;
  for ( int Np = 0; Np <= N; Np++ ){//* Loop over all N for phonon
    std::vector<int> Ivec(L, 0);
    Ivec.at(0) = Np;
    work.push_back( Ivec );
    BTags.push_back( BosonBasisTag(Ivec) );
    while ( Ivec[L - 1] < Np ) {//* Go over all possible phonon configurations in Np phonon sector
      for (ptrdiff_t cnt = L - 2; cnt > -1; cnt--) {
        if ( Ivec[cnt] != 0 ){
          IdxK = cnt;
          break;
        }
      }
      Ivec[IdxK] = Ivec[IdxK] - 1;
      NWork = std::accumulate(Ivec.begin(), Ivec.begin() + IdxK + 1, 0);
      Ivec.at(IdxK+1) = (int)(Np - NWork);
      if ( IdxK < L - 2 ){
        for (ptrdiff_t cnt = IdxK + 2; cnt < L; cnt++) {
          Ivec.at(cnt) = 0;
        }
      }
      NWork = std::accumulate(Ivec.begin(), Ivec.end(), 0);
      assert( Np == NWork );
      work.push_back( Ivec );
      BTags.push_back( BosonBasisTag(Ivec) );
    }
  }
  //* sort BTags
  BStates = SortBTags( work, BTags );
  assert( BStates.size() == BTags.size() );
}

//* Holstein phonon in momentum space (total momentum conserved)
void Basis::PhononK(const std::vector<int> Kn, const int TargetK, const int WithoutK0Phonon ){
  HaveU1 = false;
  assert( !(isFermion) );
  assert( WithoutK0Phonon == 0 || WithoutK0Phonon == 1);
  std::vector< std::vector<int> > work;//* This is final basis
  int IdxK = 0, NWork = 0;
  for ( int Kf = 0; Kf < Kn.size(); Kf++ ){//* Fermion momentum index
    for ( int Np = 0; Np <= N; Np++ ){//* Loop over all N for phonon
      std::vector<int> Ivec(L - WithoutK0Phonon, 0);
      Ivec.at(0) = Np;
      int KTotal = Kn.at(Kf);
      for ( int kpi = 0; kpi < Ivec.size(); kpi++ ) KTotal += Kn.at(kpi+WithoutK0Phonon) * Ivec.at(kpi);//* If WithoutK0Phono, skip k=0
      // std::cout << Kf << " " << std::flush;
      // PrintVector(Ivec, 2, " ");
      // std::cout << "Kt = " << KTotal << std::flush;
      int KIndex = DeltaKIndex( KTotal, Kn );//* FBZ
      // std::cout << " to " << KTotal << std::endl;
      if ( KTotal == TargetK ){
        std::vector<int> tmp = Ivec;
        std::vector<int>::iterator it = tmp.begin();
        tmp.insert ( it , Kf );
        assert( tmp.size() == L + 1 - WithoutK0Phonon );
        work.push_back( tmp );
        BTags.push_back( BosonBasisTag(tmp) );
        // PrintVector(tmp, 3, " ");
      }
      while ( Ivec[L - WithoutK0Phonon - 1] < Np ) {//* Go over all possible phonon configurations in Np phonon sector
        //? For instance, L = 4, N = 2
        //? (2,0,0,0)
        //? (1,1,0,0)
        //? (1,0,1,0)
        //? (1,0,0,1)
        //? (0,2,0,0)
        //? (0,1,1,0)
        //? (0,1,0,1)
        //? (0,0,2,0)
        //? (0,0,1,1)
        //? (0,0,0,2)
        for (ptrdiff_t cnt = L - WithoutK0Phonon - 2; cnt > -1; cnt--) {
          if ( Ivec[cnt] != 0 ){
            IdxK = cnt;
            break;
          }
        }
        ////std::cout << "IdxK=" << IdxK << std::endl;
        Ivec[IdxK] = Ivec[IdxK] - 1;
        NWork = std::accumulate(Ivec.begin(), Ivec.begin() + IdxK + 1, 0);
        Ivec.at(IdxK+1) = (int)(Np - NWork);
        if ( IdxK < L - WithoutK0Phonon - 2 ){
          for (ptrdiff_t cnt = IdxK + 2; cnt < L - WithoutK0Phonon; cnt++) {
            Ivec.at(cnt) = 0;
          }
        }
        ////std::cout << Kn.size() << " " << Ivec.size() << std::endl;
        NWork = std::accumulate(Ivec.begin(), Ivec.end(), 0);
        assert( Np == NWork );
        KTotal = Kn.at(Kf);
        for ( int kpi = 0; kpi < Ivec.size(); kpi++ ) KTotal += Kn.at(kpi+WithoutK0Phonon) * Ivec.at(kpi);//* If WithoutK0Phono, skip k=0
        int KIndex = DeltaKIndex( KTotal, Kn );//* FBZ
        if ( KTotal == TargetK ){
          std::vector<int> tmp = Ivec;
          std::vector<int>::iterator it = tmp.begin();
          tmp.insert ( it , Kf );
          assert( tmp.size() == L - WithoutK0Phonon + 1 );
          work.push_back( tmp );
          BTags.push_back( BosonBasisTag(tmp) );
          // PrintVector(tmp, 3, " ");
        }
      }
    }
  }
  //* sort BTags
  BStates = SortBTags( work, BTags );
  assert( BStates.size() == BTags.size() );
}
