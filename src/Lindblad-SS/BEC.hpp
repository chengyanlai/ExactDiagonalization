#ifndef __RUNGE_KUTTA_HPP__
#define __RUNGE_KUTTA_HPP__
#include <tuple>
#include "src/EDType.hpp"
#include "src/EigenMatrix.hpp"
#include "src/Hamiltonian/Hamiltonian.hpp"

void Lindblad_RK4( const RealType &dt,
  const std::vector< std::vector< std::tuple<ptrdiff_t,RealType> > > Idx1,
  const std::vector< std::tuple<int,int> > Idx2,
  const std::vector<RealType> Gammas, const std::vector<Basis> &bas,
  const Hamiltonian<ComplexType> &ham, ComplexMatrixType &Rhos);

void Lindblad_Newton( const RealType &dt,
  const std::vector< std::vector< std::tuple<ptrdiff_t,RealType> > > Idx1,
  const std::vector< std::tuple<int,int> > Idx2,
  const std::vector<RealType> Gammas,
  const std::vector<Basis> &bas,
  const Hamiltonian<ComplexType> &ham,
  ComplexMatrixType &Rhos);

ComplexMatrixType Lindblad1( const std::vector< std::vector< std::tuple<ptrdiff_t,RealType> > > Idx,
  const std::vector<RealType> Gammas,
  const Basis &bs, const ComplexMatrixType &rho);

ComplexMatrixType Lindblad2( const std::vector< std::tuple<int,int> > Idx,
  const std::vector<RealType> Gammas,
  const Basis &bs, const ComplexMatrixType &rho);

#endif /* end of include guard: __RUNGE_KUTTA_HPP__ */
