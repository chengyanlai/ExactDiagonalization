#include <stddef.h>
#include <cmath>
#include <numeric>//accumulate, iota
#include "src/EDType.hpp"
#include "src/Basis/Basis.hpp"

void Basis::Boson(){
  /* This is boson basis with single U(1) sector */
  HaveU1 = true;
  assert( !(isFermion) );
  std::vector< std::vector<int> > work;
  int k, NWork;
  std::vector<int> Ivec(L, 0);
  /* NOTE: Start with (N, 0, 0, 0, ...) */
  Ivec.at(0) = N;
  work.push_back( Ivec );
  BTags.push_back( BosonBasisTag(Ivec) );
  while ( Ivec[L - 1] < N ) {
    /* NOTE: Stop at (0, 0, 0, ..., N) */
    for (ptrdiff_t cnt = L - 2; cnt > -1; cnt--) {
      if ( Ivec[cnt] != 0 ){
        k = cnt;
        break;//break for
        /* ( ......, ?, 0, 0, ?)
                     ^ this is `k`
          1st time, k = 0.
          2nd time, k = 1.
              .
              .
              .
                  , k = 0.
        */
      }
    }
    Ivec[k] = Ivec[k] - 1;
    /* 1st time, (N-1, 0, 0, 0, .... )
       2nd time, (N-1, 0, 0, 0, .... )
        .
        .
        .
               , (N-2, 0, 0, 0, .... )
    */
    NWork = std::accumulate(Ivec.begin(), Ivec.begin() + k + 1, 0);
    /* NWork is the rest of the particle not in the last site - 1
       1st time, NWork = N-1
       2nd time, NWork = N-1
        .
        .
        .
               , NWork = N-2
    */
    Ivec.at(k+1) = (int)(N - NWork);
    /* 1st time, (N-1, 1, 0, 0, .... )
       2nd time, (N-1, 0, 1, 0, .... )
        .
        .
        .      , (N-1, 0, 0, .... , 1)
               , (N-2, 2, 0, .... , 1)
    */
    if ( k < L - 2 ){
      for (ptrdiff_t cnt = k + 2; cnt < L; cnt++) {
        Ivec.at(cnt) = 0;
      }
      /* 1st time, (N-1, 1, 0, 0, .... )
         2nd time, (N-1, 0, 1, 0, .... )
          .
          .
          .      , (N-1, 0, 0, .... , 1)
                 , (N-2, 2, 0, .... , 0)
      */
    }
    work.push_back( Ivec );
    // PrintBosonBasis(Ivec);
    BTags.push_back( BosonBasisTag(Ivec) );
    NWork = std::accumulate(Ivec.begin(), Ivec.end(), 0);
    assert( N == NWork );
  }
  // sort BTags
  BStates = SortBTags( work, BTags );
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

template <typename T>// Small to big
std::vector<std::vector<int> > SortBTags( const std::vector<std::vector<int> >& st, std::vector<T> &v ) {
  assert( st.size() == v.size() );
  std::vector<std::vector<int> > NewSt;
  NewSt.clear();
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
  for ( auto i : idx){
    NewSt.push_back( st.at(i) );
  }
  return NewSt;
}
template std::vector<std::vector<int> > SortBTags( const std::vector<std::vector<int> >& st, std::vector<RealType> &v );
