#include "src/Node/Node.hpp"

// template<typename Tnum, class T>
// int Node<Tnum, T>::NumNodes = 0;

template<typename Tnum, class T>
Node<Tnum, T>::Node()
{
  // this->data = NumNodes;
  this->NumLinks = 0;
  // this->NumNodes++;
}

template<typename Tnum, class T>
Node<Tnum, T>::Node(const T& item, Node<Tnum, T>* p, Tnum J)
{
  this->NumLinks = 0;
  // this->NumNodes++;
  this->data = item;
  if ( p != NULL ){
    this->NumLinks++;
    this->Neighbor.push_back(p);
    this->Jval.push_back(J);
    p->NumLinks++;
    p->Neighbor.push_back(this);
    p->Jval.push_back(J);
  }
}

template<typename Tnum, class T>
Node<Tnum, T>::~Node(){}

template<typename Tnum, class T>
void Node<Tnum, T>::LinkTo(Node<Tnum, T>* p, Tnum J)
{
  this->NumLinks++;
  this->Neighbor.push_back(p);
  this->Jval.push_back(J);
  p->NumLinks++;
  p->Neighbor.push_back(this);
  p->Jval.push_back(J);
}

template class Node<RealType, int>;
