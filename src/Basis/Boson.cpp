#include <cmath>
#include <numeric>//accumulate, iota
#include "src/EDType.hpp"
#include "src/Basis/Basis.hpp"

void Basis::Boson()
{
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
