#include <iostream>
#include <vector>
#include "src/Node/Node.hpp"

int main(int argc, char const *argv[]) {
  Node<double, int> *A = new Node<double, int>(1);
  Node<double, int> *B = new Node<double, int>(2, A);

  std::cout << A->data << " " << A->NumNeighbors() << " " <<
    A->VerifySite() << " " << std::flush;
  for (auto &j : A->getNeighbors() ){
    std::cout << j->data << std::endl;
  }
  std::cout << B->data << " " << B->NumNeighbors() << " " <<
    B->VerifySite() << " " << std::flush;
  for (auto &j : B->getNeighbors() ){
    std::cout << j->data << std::endl;
  }
  return 0;
}
