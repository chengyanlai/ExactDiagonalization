#ifndef __RUNGE_KUTTA_HPP__
#define __RUNGE_KUTTA_HPP__
#include <tuple>
#include "src/bitwise.h"
#include "src/EDType.hpp"
#include "src/EigenMatrix.hpp"
#include "src/Hamiltonian/Hamiltonian.hpp"

void FRK4( const RealType &dt, const std::vector<RealType> &gammas,
  const std::vector<std::tuple<int,int,int> > &SiteTypesSpin,
  const std::vector<std::vector<std::pair<int,int> > > &BasisIds,
  const std::vector<std::vector<std::vector<std::pair<size_t, size_t> > > > &CollapseIds,
  const std::vector<std::vector<Basis> > &bas,
  const std::vector<Hamiltonian<ComplexType> > &hams,
  std::vector<ComplexMatrixType> &Rhos );
void FNewton( const RealType &dt, const std::vector<RealType> &gammas,
  const std::vector<std::tuple<int,int,int> > &SiteTypesSpin,
  const std::vector<std::vector<std::pair<int,int> > > &BasisIds,
  const std::vector<std::vector<std::vector<std::pair<size_t, size_t> > > > &CollapseIds,
  const std::vector<std::vector<Basis> > &bas,
  const std::vector<Hamiltonian<ComplexType> > &hams,
  std::vector<ComplexMatrixType> &Rhos );
ComplexMatrixType LindbladTerm( const size_t dim, const std::vector<std::pair<size_t, size_t> > &CIdx, const ComplexMatrixType &MapMat);
ComplexMatrixType AntiCommutatorF( const int type, const size_t Site, const size_t Spin,
  const std::vector<Basis> Bs, Hamiltonian<ComplexType> ham, const ComplexMatrixType &rho);
void Cfdagger(const size_t Site, const int Spin,
  const std::vector<std::vector<Basis> > &Bs, std::vector<Hamiltonian<ComplexType> > &hams,
  const std::map<std::pair<int,int>, int > PairIndex1,
  const std::map<int, std::pair<int,int> > PairIndex2,
  std::vector<std::pair<int,int> > &BasisIdx, std::vector<std::vector<std::pair<size_t, size_t> > > &CollapseIdx);
void Cf(const size_t Site, const int Spin,
  const std::vector<std::vector<Basis> > &Bs, std::vector<Hamiltonian<ComplexType> > &hams,
  const std::map<std::pair<int,int>, int > PairIndex1,
  const std::map<int, std::pair<int,int> > PairIndex2,
  std::vector<std::pair<int,int> > &BasisIdx, std::vector<std::vector<std::pair<size_t, size_t> > > &CollapseIdx);
#endif /* end of include guard: __RUNGE_KUTTA_HPP__ */
