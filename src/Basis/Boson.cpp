#include <stddef.h>
#include <cmath>
#include <numeric>//accumulate, iota
#include "src/EDType.hpp"
#include "src/Basis/Basis.hpp"

void Basis::Boson(){
  HaveU1 = true;
  assert( !(isFermion) );
  std::vector< std::vector<int> > work;
  int k, Nres;
  std::vector<int> Ivec(L, 0);
  Ivec.at(0) = N;
  work.push_back( Ivec );
  BTags.push_back( BosonBasisTag(Ivec) );
  while ( Ivec[L - 1] < N ) {
    for (ptrdiff_t cnt = L - 2; cnt > -1; cnt--) {
      if ( Ivec[cnt] != 0 ){
        k = cnt;
        break;//break for
      }
    }
    Ivec[k] = Ivec[k] - 1;
    Nres = std::accumulate(Ivec.begin(), Ivec.begin() + k + 1, 0);
    Ivec.at(k+1) = (int)(N - Nres);
    if ( k < L - 2 ){
      for (ptrdiff_t cnt = k + 2; cnt < L; cnt++) {
        Ivec.at(cnt) = 0;
      }
    }
    work.push_back( Ivec );
    BTags.push_back( BosonBasisTag(Ivec) );
    Nres = std::accumulate(Ivec.begin(), Ivec.end(), 0);
    assert( N == Nres );
  }
  // sort BTags
  std::vector<size_t> NewIdx = SortBTags( BTags );
  // put work in order.
  for ( auto i : NewIdx){
    BStates.push_back( work.at(i) );
  }
  assert( BStates.size() == BTags.size() );
}

void Basis::Boson( const int MaxLocalB ){
  HaveU1 = false;
  /* Build bosonic basis from zero to maximum boson number allowed locally.
    Total boson number is not larger than N */
  assert( !(isFermion) );
  assert( L * MaxLocalB >= N );
  std::vector< std::vector<int> > work;
  int k, Nres;
  /* NOTE: Include all U(1) sector which has particle number smaller than N */
  std::vector<int> Ivec(L,0);
  work.push_back( Ivec );
  BTags.push_back( BosonBasisTag(Ivec) );
  for (ptrdiff_t cntN = N; cntN > 0; cntN--) {
    Ivec.assign(L, 0);
    Ivec.at(0) = cntN;
    if( cntN <= MaxLocalB ){
      work.push_back( Ivec );
      BTags.push_back( BosonBasisTag(Ivec) );
    }
    while ( Ivec[L - 1] < cntN ) {
      for (ptrdiff_t cnt = L - 2; cnt > -1; cnt--) {
        if ( Ivec[cnt] != 0 ){
          k = cnt;
          break;//break for-loop
        }
      }
      Ivec[k] = Ivec[k] - 1;
      Nres = std::accumulate(Ivec.begin(), Ivec.begin() + k + 1, 0);
      Ivec.at(k+1) = (int)(cntN - Nres);
      if ( k < L - 2 ){
        for (ptrdiff_t cnt = k + 2; cnt < L; cnt++) {
          Ivec.at(cnt) = 0;
        }
      }
      bool isBasis = true;
      for ( auto n : Ivec){
        if ( n > MaxLocalB ){
          isBasis = false;
          break;
        }
      }
      if ( isBasis ) {
        work.push_back( Ivec );
        BTags.push_back( BosonBasisTag(Ivec) );
        Nres = std::accumulate(Ivec.begin(), Ivec.end(), 0);
        assert( cntN == Nres );
      }
    }
  }
  // sort BTags
  std::vector<size_t> NewIdx = SortBTags( BTags );
  // put work in order.
  for ( auto i : NewIdx ){
    BStates.push_back( work.at(i) );
  }
  assert( BStates.size() == BTags.size() );
}

void Basis::BosonTB( const size_t TBloc, const bool HARD_CUT ){
  assert( !(isFermion) );
  /* NOTE: Terminator Beam can not be located in first site or larger than system size */
  assert( TBloc > 0 );
  assert( TBloc < L );
  std::vector< std::vector<int> > work;
  int k, Nres;
  /* NOTE: Include all U(1) sector which has particle number smaller than N */
  std::vector<int> Ivec(L,0);
  work.push_back( Ivec );
  BTags.push_back( BosonBasisTag(Ivec) );
  for (ptrdiff_t cntN = N; cntN > 0; cntN--) {
    Ivec.assign(L, 0);
    Ivec.at(0) = cntN;
    work.push_back( Ivec );
    BTags.push_back( BosonBasisTag(Ivec) );
    while ( Ivec[L - 1] < cntN ) {
      for (ptrdiff_t cnt = L - 2; cnt > -1; cnt--) {
        if ( Ivec[cnt] != 0 ){
          k = cnt;
          break;//break for-loop
        }
      }
      Ivec[k] = Ivec[k] - 1;
      Nres = std::accumulate(Ivec.begin(), Ivec.begin() + k + 1, 0);
      Ivec.at(k+1) = (int)(cntN - Nres);
      if ( k < L - 2 ){
        for (ptrdiff_t cnt = k + 2; cnt < L; cnt++) {
          Ivec.at(cnt) = 0;
        }
      }
      /* NOTE: Rule out the state which is not existed in terminator beam setup */
      if ( !(HARD_CUT) || (Ivec.at(TBloc) < 2 && HARD_CUT) ) {
        work.push_back( Ivec );
        BTags.push_back( BosonBasisTag(Ivec) );
        Nres = std::accumulate(Ivec.begin(), Ivec.end(), 0);
        assert( cntN == Nres );
      }
    }
  }
  // sort BTags
  std::vector<size_t> NewIdx = SortBTags( BTags );
  // put work in order.
  for ( auto i : NewIdx ){
    BStates.push_back( work.at(i) );
  }
  assert( BStates.size() == BTags.size() );
}

/* This has to be unique number to tag boson basis state. */
RealType BosonBasisTag( const std::vector<int> vec ){
  RealType tag = 0.0e0;
  int cnt = 0;
  for (auto val : vec) {
    tag += sqrt(1.0e+2 * (RealType)cnt + 3.0e0 ) * (RealType)val;
    cnt++;
  }
  return tag;
}

template <typename T>
std::vector<size_t> SortBTags( std::vector<T> &v ) {
  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota( idx.begin(), idx.end(), 0 );

  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  // put v in order.
  std::vector<T> new_v;
  for ( auto i : idx ){
    new_v.push_back( v.at(i) );
  }
  v = new_v;
  assert( v.size() == idx.size() );
  return idx;
}
template std::vector<size_t> SortBTags( std::vector<RealType> &v );

void Basis::printBosonBasis( const std::vector<int> state )const
{
  for( auto i : state ){
    std::cout << i << ", " << std::flush;
  }
  std::cout << std::endl;
}

