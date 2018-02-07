#include <stddef.h>
#include <cmath>
#include <numeric>//accumulate, iota
#include "src/EDType.hpp"
#include "src/Basis/Basis.hpp"

static bool CheckState( const int N, const std::vector<int> Input ){
  bool Allow = true;
  int Nt = std::accumulate(Input.begin(), Input.end(), 0);
  int Nused = 0;
  int Nc = 0;
  typename std::vector<int>::const_iterator it = Input.begin();
  for (;it != Input.end(); ++it){
    Nused += *it;
    if ( Nused > N ){
      Allow = false;
      break;
    }
    Nused += 1;
    Nc += *it;
    if ( Nc == Nt ) break;
  }
  return Allow;
}

void Basis::Phonon(){
  /* Limited function space with one fermion at site-0 */
  HaveU1 = false;
  assert( !(isFermion) );
  std::vector< std::vector<int> > work;
  int k, NWork;
  for ( size_t Np = 0; Np <= N; Np++ ){
    std::vector<int> Ivec(L, 0);
    Ivec.at(0) = Np;
    work.push_back( Ivec );
PrintBosonBasis(Ivec);
    BTags.push_back( BosonBasisTag(Ivec) );
    while ( Ivec[L - 1] < Np ) {
      for (ptrdiff_t cnt = L - 2; cnt > -1; cnt--) {
        if ( Ivec[cnt] != 0 ){
          k = cnt;
          break;
        }
      }
      Ivec[k] = Ivec[k] - 1;
      NWork = std::accumulate(Ivec.begin(), Ivec.begin() + k + 1, 0);
      Ivec.at(k+1) = (int)(Np - NWork);
      if ( k < L - 2 ){
        for (ptrdiff_t cnt = k + 2; cnt < L; cnt++) {
          Ivec.at(cnt) = 0;
        }
      }
      if ( CheckState(N, Ivec) ){
        work.push_back( Ivec );
PrintBosonBasis(Ivec);
        BTags.push_back( BosonBasisTag(Ivec) );
        NWork = std::accumulate(Ivec.begin(), Ivec.end(), 0);
        assert( Np == NWork );
      }
    }
  }
  // sort BTags
  std::vector<size_t> NewIdx = SortBTags( BTags );
  // put work in order.
  for ( auto i : NewIdx){
    BStates.push_back( work.at(i) );
  }
  assert( BStates.size() == BTags.size() );
}
