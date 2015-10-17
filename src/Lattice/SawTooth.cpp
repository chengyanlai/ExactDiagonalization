#include <cassert>
#include "src/EDType.hpp"
#include "src/Lattice/Preset.hpp"

std::vector< Node<RealType, int>* > SawTooth(const int &L,
  const std::vector<RealType> JAB, const std::vector<RealType> JAA,
  const bool OBC){
  if ( OBC ){
    assert( L % 2 == 1 );
    assert( JAB.size() == L - 1);
    assert( JAA.size() == (L - 1) / 2);
  }else{
    assert( L % 2 == 0 );
    assert( JAB.size() == L );
    assert( JAA.size() == L/2 );
  }
  std::vector< Node<RealType, int>* > lattice;
  int cnt = 0;
  int cntAA = 0;
  while ( cnt < L ) {
    if ( lattice.size() == 0 ) {
      Node<RealType, int> *A = new Node<RealType, int>(cnt);
      lattice.push_back(A);
    } else if ( cnt % 2 == 1 ) {
      Node<RealType, int> *B = new Node<RealType, int>(cnt, lattice[cnt-1], JAB[cnt-1]);
      lattice.push_back(B);
    } else {
      Node<RealType, int> *A = new Node<RealType, int>(cnt, lattice[cnt-1], JAB[cnt-1]);
      A->LinkTo(lattice[cnt-2], JAA[cntAA]);
      lattice.push_back(A);
      cntAA++;
    }
    cnt++;
  }
  if ( !(OBC) ){
    lattice[0]->LinkTo( lattice[L-1], JAB[L-1] );
    lattice[0]->LinkTo( lattice[L-2], JAA[L/2-1] );
  }
  assert( lattice.size() == L );
  return lattice;
}
