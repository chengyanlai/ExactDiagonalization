#include <stddef.h>
#include <cmath>
#include <numeric>//accumulate, iota
#include <tuple>
#include "src/EDType.hpp"
#include "src/Basis/Basis.hpp"

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

bool Basis::CheckExist2( const std::vector<int>& State, const std::vector<std::vector<int> > States ){
  typename std::vector<std::vector<int> >::const_reverse_iterator it = States.rbegin();
  for(; it != States.rend(); ++it){
    if ( State == *it ) return true;
  }
  return false;
}

bool Basis::CheckExist1( const std::vector<int>& State ){
  RealType tg = BosonBasisTag(State);
  size_t Idx = this->GetIndexFromTag( tg );
  return State == BStates.at(Idx);
//   if ( State == BStates.at(Idx) ){
// // std::cout << "new\n";
// // PrintBosonBasis(State);
// // std::cout << "cp to\n";
// // PrintBosonBasis(BStates.at(Idx));
// // std::cout << (State == BStates.at(Idx)) << std::endl;
// // std::cin.get();
//     return true;
//   }else{
// std::cout << "new\n";
// PrintBosonBasis(State);
// std::cout << "cp to\n";
// PrintBosonBasis(BStates.at(Idx));
// std::cout << (State == BStates.at(Idx)) << std::endl;
// // std::cin.get();
//     return false;
//   }
}

std::vector<std::vector<int> > Basis::ApplyOffdiagonal( const std::vector<std::vector<int> > InputStates ){
  // States has phonon configurations with fermion at site-0
  // Btags need to be sorted.
  std::vector<std::vector<int> > NewStates;
  NewStates.clear();
  typename std::vector<std::vector<int> >::const_iterator it = InputStates.begin();
  for (; it != InputStates.end(); ++it ){
    // Create Phonon
    std::vector<int> State1 = *it;
    State1.at(0) += 1;
    assert( State1.size() == L );
    if ( !CheckExist1( State1 ) && !CheckExist2( State1, NewStates ) ) NewStates.push_back(State1);
    // Fermion jump right
    std::vector<int> State2;
    State2.insert(State2.begin(), it->begin()+1, it->end());
    State2.push_back( it->at(0) );
    assert( State2.size() == L );
    if ( !CheckExist1( State2 ) && !CheckExist2( State2, NewStates ) ) NewStates.push_back(State2);
    // Fermiom jump left
    std::vector<int> State3;
    State3.push_back( it->back() );
    State3.insert(State3.end(), it->begin(), it->end()-1);
    assert( State3.size() == L );
    if ( !CheckExist1( State3 ) && !CheckExist2( State3, NewStates ) ) NewStates.push_back(State3);
  }
  return NewStates;
}

void Basis::Phonon(){
  /* Limited function space with one fermion at site-0 */
  HaveU1 = false;
  assert( !(isFermion) );
  BTags.clear();

  // std::vector< std::tuple<int, std::vector<int> > > lfs;// limited functional space
  // lfs.push_back( std::make_tuple(0, std::vector<int>(L,0)) );
  std::vector<int> State(L, 0);
  std::vector<std::vector<int> > NewStates;
  NewStates.push_back( State );
  // add this vacuum state
  BStates.push_back( State );
  BTags.push_back( BosonBasisTag(State) );
  for ( int cnt = 0; cnt < N; cnt++ ){
    NewStates = ApplyOffdiagonal( NewStates );
    if ( NewStates.size() ){
      typename std::vector<std::vector<int> >::const_iterator it = NewStates.begin();
      for (; it != NewStates.end(); ++it ){
        BStates.push_back( *it );
        BTags.push_back( BosonBasisTag(*it) );
      }
      BStates = SortBTags( BStates, BTags );
    }
  }
  DummyCheckState();
  assert( BStates.size() == BTags.size() );
}
