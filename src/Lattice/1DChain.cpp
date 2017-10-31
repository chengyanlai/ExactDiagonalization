#include <cassert>
#include "src/EDType.hpp"
#include "src/Lattice/preset.hpp"

template<typename T>
std::vector< Node<T>* > NN_1D_Chain(const int &L, const std::vector<T> J, const bool OBC){
  if ( OBC ){
    assert( J.size() == L - 1);
  }else{
    assert( J.size() == L );
  }
  std::vector< Node<T>* > lattice;
  int cnt = 0;
  while ( cnt < L ) {
    if ( lattice.size() == 0 ) {
      Node<T> *A = new Node<T>(cnt);
      lattice.push_back(A);
    } else {
      Node<T> *A = new Node<T>(cnt, lattice[cnt-1], J[cnt-1]);
      lattice.push_back(A);
    }
    cnt++;
  }
  if ( not OBC ) lattice[0]->LinkTo( lattice[L-1], J[L-1]);
  assert( lattice.size() == L );
  return lattice;
}
template std::vector< Node<RealType>* > NN_1D_Chain(const int &L, const std::vector<RealType> J, const bool OBC);
template std::vector< Node<ComplexType>* > NN_1D_Chain(const int &L, const std::vector<ComplexType> J, const bool OBC);
