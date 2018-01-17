#ifndef __PRESET_HPP__
#define __PRESET_HPP__
#include <vector>
#include "src/EDType.hpp"
#include "src/Node/Node.hpp"

template<typename T>
std::vector< Node<T>* > NN_1D_Chain(const int &L, const std::vector<T> J, const bool OBC = true);
template<typename T>
std::vector< Node<T>* > NNN_1D_Chain(const int &L, const std::vector<T> J, const bool OBC = true);

template<typename T>
std::vector< Node<T>* > SawTooth(const int &L, const std::vector<T> J1, const std::vector<T> J2, const bool OBC = true);

template<typename T>
std::vector< Node<T>* > CreutzLadder(const int &LUnit, const std::vector<T> JAA, const std::vector<T> JBB, const std::vector<T> JdAB, const std::vector<T> JdBA, const std::vector<T> JvAB, const bool OBC = true);
#endif//__PRESET_HPP__
