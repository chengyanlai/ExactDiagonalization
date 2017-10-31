#include <cassert>
#include "src/EDType.hpp"
#include "src/Lattice/preset.hpp"

template<typename T>
std::vector< Node<T>* > CreutzLadder(const int &LUnit, const std::vector<T> JAA, const std::vector<T> JBB, const std::vector<T> JdAB, const std::vector<T> JdBA, const std::vector<T> JvAB, const bool OBC){
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
  std::vector< Node<T>* > lattice;
  int cnt = 0;
  while ( cnt < LUnit ) {
    if ( cnt == 0 ){
      Node<T> *A = new Node<T>(2 * cnt, NULL, 0.0e0, "A");
      lattice.push_back(A);
      Node<T> *B = new Node<T>(2 * cnt + 1, NULL, 0.0e0, "B");
      B->LinkTo( lattice[2*cnt], JvAB[cnt]);//link vertical AB bond
      lattice.push_back(B);
    }else if ( cnt > 0 ) {
      Node<T> *A = new Node<T>(2 * cnt, lattice[2*(cnt-1)], JAA[cnt-1], "A");//link NN AA bond
      A->LinkTo( lattice[2*(cnt-1)+1], JdAB[cnt-1]);//link diagonal AB(look from right to left) bond
      lattice.push_back(A);
      Node<T> *B = new Node<T>(2 * cnt + 1, lattice[2*(cnt-1)+1], JBB[cnt-1], "B");//link NN BB bond
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
template std::vector< Node<RealType>* > CreutzLadder(const int &LUnit, const std::vector<RealType> JAA, const std::vector<RealType> JBB, const std::vector<RealType> JdAB, const std::vector<RealType> JdBA, const std::vector<RealType> JvAB, const bool OBC);
template std::vector< Node<ComplexType>* > CreutzLadder(const int &LUnit, const std::vector<ComplexType> JAA, const std::vector<ComplexType> JBB, const std::vector<ComplexType> JdAB, const std::vector<ComplexType> JdBA, const std::vector<ComplexType> JvAB, const bool OBC);
