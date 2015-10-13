#ifndef __PRESET_HPP__
#define __PRESET_HPP__
#include <vector>
#include "src/EDType.hpp"
#include "src/Node/Node.hpp"

std::vector< Node<RealType, int>* > NN_1D_Chain(const int &L,
  const std::vector<RealType> J, const bool OBC = true);
#endif//__PRESET_HPP__
