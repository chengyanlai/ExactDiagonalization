#ifndef __RUNGE_KUTTA_HPP__
#define __RUNGE_KUTTA_HPP__
#include "src/EDType.hpp"
#include "src/EigenMatrix.hpp"
#include "src/Hamiltonian/Hamiltonian.hpp"

void Lindblad_RK4( const RealType &dt, const RealType &gamma,
  const size_t TBloc, const std::vector<Basis> &bas,
  const std::vector<Hamiltonian<ComplexType> > &ham,
  const std::vector<std::vector<size_t> > &CIdx,
  std::vector<ComplexMatrixType> &Rhos);

void Lindblad_Newton( const RealType &dt, const RealType &gamma,
  const size_t TBloc, const std::vector<Basis> &bas,
  const std::vector<Hamiltonian<ComplexType> > &ham,
  const std::vector<std::vector<size_t> > &CIdx,
  std::vector<ComplexMatrixType> &Rhos);

ComplexMatrixType Lindblad1(const size_t TBloc, const ComplexMatrixType &MapMat,
  const Basis &bs, const std::vector<size_t> &CIdx);

ComplexMatrixType Lindblad2( const size_t TBloc, const Basis &bs,
  const ComplexMatrixType &rho);

#endif /* end of include guard: __RUNGE_KUTTA_HPP__ */
