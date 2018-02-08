#include <stddef.h>
#include <cmath>
#include "src/EDType.hpp"
#include "src/Basis/Basis.hpp"

/* DEBUG use */
static void PrintBosonBasis( const std::vector<int> state ){
  typename std::vector<int>::const_iterator it = state.begin();
  for (;it != state.end(); ++it ){
    std::cout << *it << ", " << std::flush;
  }
  std::cout << std::endl;
};

RealType Basis::CreatePhonon( std::vector<int>& state, const int site )const{
  state.at(site) += 1;
  return BosonBasisTag(state);
}

RealType Basis::DestroyPhonon( std::vector<int>& state, const int site )const{
  if ( state.at(site) ){
    state.at(site) = state.at(site) - 1;
    return BosonBasisTag(state);
  }else{
    return -1;
  }
}

RealType Basis::FermionJumpRight( std::vector<int>& state, const int NumJumps )const{
  int L = state.size();
  std::vector<int> ns;
  ns.clear();
  ns.insert(ns.begin(), state.begin()+NumJumps, state.end() );
  ns.insert(ns.end(), state.begin(), state.begin()+NumJumps );
  // ns.push_back( state.at(0) );// if NumJumps = 1
// std::cout << "R - " << NumJumps << "\n";
// PrintBosonBasis(state);
// PrintBosonBasis(ns);
  state = ns;
  assert( state.size() == L );
  return BosonBasisTag(state);
}

RealType Basis::FermionJumpLeft( std::vector<int>& state, const int NumJumps )const{
  int L = state.size();
  std::vector<int> ns;
  ns.clear();
  // ns.push_back( state.back() );// if NumJumps = 1
  ns.insert(ns.begin(), state.end()-NumJumps, state.end() );
  ns.insert(ns.end(), state.begin(), state.end()-NumJumps);
// std::cout << "L - " << NumJumps << "\n";
// PrintBosonBasis(state);
// PrintBosonBasis(ns);
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
        // std::cout << std::endl;
        std::cin.get();
      }
    }
  }
}

std::vector<std::vector<int> > Basis::ApplyOffdiagonal( const std::vector<std::vector<int> >& InputStates ){
  // States has phonon configurations with fermion at site-0
  // Btags need to be sorted.
  std::vector<std::vector<int> > NewStates;
  NewStates.clear();
  typename std::vector<std::vector<int> >::const_iterator it = InputStates.begin();
  for (; it != InputStates.end(); ++it ){
    // Create Phonon
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
    // Fermion jump right
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
    // Fermiom jump left
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

void Basis::Phonon(){
  /* Limited function space with one fermion at site-0 */
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
std::cout << cnt << " " << BTags.size() << std::endl;;
  }
  // DummyCheckState();
  assert( BStates.size() == BTags.size() );
}
