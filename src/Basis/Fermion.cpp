#include <cmath>//pow
#include "src/EDType.hpp"
#include "src/bitwise.h"
#include "src/Basis/Basis.hpp"

void Basis::Fermion()
{
  assert( isFermion );
  size_t minrange = 0;
  size_t maxrange = 0;
  for (size_t cnt = 0; cnt < N; cnt++) {
    minrange += pow(2, cnt);
    maxrange += pow(2, L - cnt - 1);
  }
  size_t sub_cnt = 0;
  std::vector<size_t> windex (maxrange+1, 0);
  for (size_t cnt1 = minrange; cnt1 <= maxrange; cnt1++) {
    int nbit = 0;
    for (size_t cnt2 = 0; cnt2 < L; cnt2++) {
      if ( btest(cnt1, cnt2) ){
        nbit += 1;
      }
    }
    if (nbit == N){
      sub_cnt += 1;
      FStates.push_back(cnt1);
      windex.at(cnt1) = sub_cnt - 1;//NOTE: I start from 0!
      // INFO(cnt1 << " " << sub_cnt << " " << windex.at(cnt1));
    }
  }
  FTags = windex;
}

void Basis::Fermion( const int MinF )
{
  // This will include all possible fermion numbers, the input parameter has no meaning for now.
  assert( isFermion );
  // assert( L % 2 == 0 );
  size_t minrange = 0;
  size_t maxrange = pow(2, L);
  NTotal.clear();
  for (size_t cnt1 = minrange; cnt1 < maxrange; cnt1++) {
    FStates.push_back(cnt1);
    FTags.push_back(cnt1);
    int nbit = 0;
    for (size_t cnt2 = 0; cnt2 < L; cnt2++) {
      if ( btest(cnt1, cnt2) ){
        nbit += 1;
      }
    }
    NTotal.push_back(nbit);
  }
}

void Basis::printFermionBasis( const int state )const{
  for (size_t cnt = 0; cnt < L; cnt++) {
    INFO_NONEWLINE( btest(state, cnt) << ", " );
  }
}
