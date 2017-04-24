#include <cassert>
#include "src/EDType.hpp"
#include "src/Lattice/preset.hpp"

std::vector< Node<RealType>* > CreutzLadder(const int &LUnit,
  const std::vector<RealType> JAA, const std::vector<RealType> JBB,
  const std::vector<RealType> JdAB, const std::vector<RealType> JdBA,
  const std::vector<RealType> JvAB,
  const bool OBC){
  if ( OBC ){
    assert( JAA.size() == LUnit - 1);
    assert( JBB.size() == LUnit - 1);
    assert( JdAB.size() == LUnit - 1);
    assert( JdBA.size() == LUnit - 1);
  }else{
    assert( JAA.size() == LUnit);
    assert( JBB.size() == LUnit);
    assert( JdAB.size() == LUnit);
    assert( JdBA.size() == LUnit);
  }
  assert( JvAB.size() == LUnit);
  std::vector< Node<RealType>* > lattice;
  int cnt = 0;
  while ( cnt < LUnit ) {
    if ( cnt == 0 ){
      Node<RealType> *A = new Node<RealType>(2 * cnt, NULL, 0.0e0, "A");
      lattice.push_back(A);
      Node<RealType> *B = new Node<RealType>(2 * cnt + 1, NULL, 0.0e0, "B");
      B->LinkTo( lattice[2*cnt], JvAB[cnt]);//link vertical AB bond
      lattice.push_back(B);
    }else if ( cnt > 0 ) {
      Node<RealType> *A = new Node<RealType>(2 * cnt, lattice[2*(cnt-1)], JAA[cnt-1], "A");//link NN AA bond
      A->LinkTo( lattice[2*(cnt-1)+1], JdAB[cnt-1]);//link diagonal AB(look from right to left) bond
      lattice.push_back(A);
      Node<RealType> *B = new Node<RealType>(2 * cnt + 1, lattice[2*(cnt-1)+1], JBB[cnt-1], "B");//link NN BB bond
      B->LinkTo( lattice[2*(cnt-1)], JdBA[cnt-1]);//link diagonal BA(look from right to left) bond
      B->LinkTo( lattice[2*cnt], JvAB[cnt]);//link vertical AB bond
      lattice.push_back(B);
    }
    cnt++;
  }
  if ( !(OBC) ){
    lattice[0]->LinkTo( lattice[2*(LUnit-1)], JAA[LUnit-1] );//PBC AA bond
    lattice[1]->LinkTo( lattice[2*(LUnit-1)+1], JBB[LUnit-1] );//PBC BB bond
    lattice[0]->LinkTo( lattice[2*(LUnit-1)+1], JdAB[LUnit-1] );//PBC AB bond
    lattice[1]->LinkTo( lattice[2*(LUnit-1)], JdBA[LUnit-1] );//PBC BA bond
  }
  assert( lattice.size() == 2 * LUnit );
  return lattice;
}

std::vector< Node<ComplexType>* > CreutzLadder(const int &LUnit,
  const std::vector<ComplexType> JAA, const std::vector<ComplexType> JBB,
  const std::vector<ComplexType> JdAB, const std::vector<ComplexType> JdBA,
  const std::vector<ComplexType> JvAB,
  const bool OBC){
  if ( OBC ){
    assert( JAA.size() == LUnit - 1);
    assert( JBB.size() == LUnit - 1);
    assert( JdAB.size() == LUnit - 1);
    assert( JdBA.size() == LUnit - 1);
  }else{
    assert( JAA.size() == LUnit);
    assert( JBB.size() == LUnit);
    assert( JdAB.size() == LUnit);
    assert( JdBA.size() == LUnit);
  }
  assert( JvAB.size() == LUnit);
  std::vector< Node<ComplexType>* > lattice;
  int cnt = 0;
  int cntAA = 0;
  while ( cnt < LUnit ) {
    if ( cnt == 0 ){
      Node<ComplexType> *A = new Node<ComplexType>(2 * cnt, NULL, 0.0e0, "A");
      lattice.push_back(A);
      Node<ComplexType> *B = new Node<ComplexType>(2 * cnt + 1, NULL, 0.0e0, "B");
      B->LinkTo( lattice[2*cnt], JvAB[cnt]);//link vertical AB bond
      lattice.push_back(B);
    }else if ( cnt > 0 ) {
      Node<ComplexType> *A = new Node<ComplexType>(2 * cnt, lattice[2*(cnt-1)], JAA[cnt-1], "A");//link NN AA bond
      A->LinkTo( lattice[2*(cnt-1)+1], JdAB[cnt-1]);//link diagonal AB(look from right to left) bond
      lattice.push_back(A);
      Node<ComplexType> *B = new Node<ComplexType>(2 * cnt + 1, lattice[2*(cnt-1)+1], JBB[cnt-1], "B");//link NN BB bond
      B->LinkTo( lattice[2*(cnt-1)], JdBA[cnt-1]);//link diagonal BA(look from right to left) bond
      B->LinkTo( lattice[2*cnt], JvAB[cnt]);//link vertical AB bond
      lattice.push_back(B);
    }
    cnt++;
  }
  if ( !(OBC) ){
    // NOTE: This part need to be tested with non-zero imaginary number.
    lattice[0]->LinkTo( lattice[2*(LUnit-1)], JAA[LUnit-1] );//PBC AA bond
    lattice[1]->LinkTo( lattice[2*(LUnit-1)+1], JBB[LUnit-1] );//PBC BB bond
    lattice[0]->LinkTo( lattice[2*(LUnit-1)+1], JdAB[LUnit-1] );//PBC AB bond
    lattice[1]->LinkTo( lattice[2*(LUnit-1)], JdBA[LUnit-1] );//PBC BA bond
  }
  assert( lattice.size() == 2 * LUnit );
  return lattice;
}
