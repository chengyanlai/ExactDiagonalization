#ifndef __Node_HPP__
#define __Node_HPP__
#include <vector>
#include "src/EDType.hpp"

template<typename Tnum = RealType>
class Node
{
public:
  Node();
  Node(const size_t &item, Node<Tnum>* ptrnext = NULL, Tnum J = 1.0e0);
  virtual ~Node();
  size_t data;
  void LinkTo(Node<Tnum>* p, Tnum J);
  // inline int TotalNodes()const{return this->NumNodes;};
  inline int NumNeighbors()const{return this->NumLinks;};
  inline std::vector< Tnum > getJval()const{return this->Jval;};
  inline std::vector< Node<Tnum>* > getNeighbors()const{return this->Neighbor;};
  inline bool VerifySite()const{
    INFO("Site - " << this->data << " has " << this->NumLinks << " neighbors.");
    for (size_t cnt = 0; cnt < Neighbor.size(); cnt++) {
      INFO("\t Linked to " << Neighbor.at(cnt)->data <<
        " with hopping coefficient " << Jval.at(cnt) );
    }
    return this->Neighbor.size() == this->Jval.size();
  };
  inline void Label(const size_t &la){this->data = la;};
private:
  // static int NumNodes;
  int NumLinks;
  std::vector< Node<Tnum>* > Neighbor;
  std::vector< Tnum > Jval;
};

#endif//__Node_HPP__
