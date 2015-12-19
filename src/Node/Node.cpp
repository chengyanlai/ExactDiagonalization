#include "src/Node/Node.hpp"

template<typename Tnum>
Node<Tnum>::Node()
{
  this->NumLinks = 0;
}

template<typename Tnum>
Node<Tnum>::Node(const size_t &item, Node<Tnum>* p, Tnum J)
{
  this->NumLinks = 0;
  // this->NumNodes++;
  this->data = item;
  if ( p != NULL ){
    LinkTo(p, J);
  }
}

template<typename Tnum>
Node<Tnum>::~Node(){}

template<>
void Node<RealType>::LinkTo(Node<RealType>* p, RealType J)
{
  this->NumLinks++;
  this->Neighbor.push_back(p);
  this->Jval.push_back(J);
  p->NumLinks++;
  p->Neighbor.push_back(this);
  p->Jval.push_back(J);
}

template<>
void Node<ComplexType>::LinkTo(Node<ComplexType>* p, ComplexType J)
{
  this->NumLinks++;
  this->Neighbor.push_back(p);
  this->Jval.push_back(J);
  p->NumLinks++;
  p->Neighbor.push_back(this);
  p->Jval.push_back(std::conj(J));
}

template class Node<RealType>;
template class Node<ComplexType>;
