#ifndef __Node_HPP__
#define __Node_HPP__
#include <vector>
#include "src/EDType.hpp"

template<typename Tnum = RealType, class T = int>
class Node
{
public:
  Node();
  Node(const T& item, Node<Tnum, T>* ptrnext = NULL, Tnum J = 1.0e0);
  virtual ~Node();
  T data;
  void LinkTo(Node<Tnum, T>* p, Tnum J);
  // inline int TotalNodes()const{return this->NumNodes;};
  inline int NumNeighbors()const{return this->NumLinks;};
  inline std::vector< Tnum > getJval()const{return this->Jval;};
  inline std::vector< Node<Tnum, T>* > getNeighbors()const{return this->Neighbor;};
  inline bool VerifySite()const{return this->Neighbor.size() == this->Jval.size();};
  inline void Label(const T &la){this->data = la;};
private:
  // static int NumNodes;
  int NumLinks;
  std::vector< Node<Tnum, T>* > Neighbor;
  std::vector< Tnum > Jval;
};

#endif//__Node_HPP__
