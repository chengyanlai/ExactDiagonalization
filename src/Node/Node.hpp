#ifndef __Node_HPP__
#define __Node_HPP__
#include <fstream>
#include <vector>
#include "src/EDType.hpp"

template<typename Tnum = RealType>
class Node
{
public:
  Node();
  Node(const size_t &item, Node<Tnum>* ptrnext = NULL, Tnum J = 1.0e0, std::string sublattice = "A");
  virtual ~Node();
  size_t data;
  void LinkTo(Node<Tnum>* p, Tnum J);
  // inline int TotalNodes()const{return this->NumNodes;};
  inline int NumNeighbors()const{return this->NumLinks;};
  inline std::vector< Tnum > GetJval()const{return this->Jval;};
  inline std::vector< Node<Tnum>* > GetNeighbors()const{return this->Neighbor;};
  inline bool VerifySite(std::ofstream& file)const{
    file << "Site - " << this->data << "(" << this->label << ") has " << this->NumLinks << " neighbors.\n";
    for (size_t cnt = 0; cnt < Neighbor.size(); cnt++) {
      file << "\t Linked to " << Neighbor.at(cnt)->data << " with hopping coefficient " << Jval.at(cnt) << std::endl;;
    }
    return this->Neighbor.size() == this->Jval.size();
  };
  inline void Label(const size_t &la){this->data = la;};
private:
  // static int NumNodes;
  int NumLinks;
  std::string label;
  std::vector< Node<Tnum>* > Neighbor;
  std::vector< Tnum > Jval;
};

#endif//__Node_HPP__
