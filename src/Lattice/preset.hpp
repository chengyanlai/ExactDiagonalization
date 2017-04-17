#ifndef __PRESET_HPP__
#define __PRESET_HPP__
#include <vector>
#include "src/EDType.hpp"
#include "src/Node/Node.hpp"

std::vector< Node<RealType>* > NN_1D_Chain(const int &L,
  const std::vector<RealType> J, const bool OBC = true);
std::vector< Node<ComplexType>* > NN_1D_Chain(const int &L,
  const std::vector<ComplexType> J, const bool OBC = true);
std::vector< Node<RealType>* > SawTooth(const int &L,
  const std::vector<RealType> J1, const std::vector<RealType> J2,
  const bool OBC = true);
std::vector< Node<ComplexType>* > SawTooth(const int &L,
  const std::vector<ComplexType> J1, const std::vector<ComplexType> J2,
  const bool OBC = true);
std::vector< Node<RealType>* > CreutzLadder(const int &LUnit,
  const std::vector<RealType> JAA, const std::vector<RealType> JBB,
  const std::vector<RealType> JdAB, const std::vector<RealType> JdBA,
  const std::vector<RealType> JvAB, const bool OBC = true);
std::vector< Node<ComplexType>* > CreutzLadder(const int &LUnit,
  const std::vector<ComplexType> JAA, const std::vector<ComplexType> JBB,
  const std::vector<ComplexType> JdAB, const std::vector<ComplexType> JdBA,
  const std::vector<ComplexType> JvAB, const bool OBC = true);
#endif//__PRESET_HPP__
