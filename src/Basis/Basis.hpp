#ifndef __BASIS_HPP__
#define __BASIS_HPP__
#include <iostream>
#include <vector>
#include <cassert>
#include <algorithm>//sort, distance
#include "src/EDType.hpp"
#include "src/bitwise.h"


template<class RandomIt, class T>
inline RandomIt binary_locate(RandomIt first, RandomIt last, const T& val) {
  if(val == *first) return first;
  auto d = std::distance(first, last);
  if(d==1) return first;
  auto center = (first + (d/2));
  if(val < *center) return binary_locate(first, center, val);
  return binary_locate(center, last, val);
}

class Basis
{
template<typename Tnum>
  friend class Hamiltonian;
public:
  Basis(const bool _isFermion = false);
  Basis(const size_t _L, const size_t _N, const bool _isFermion = false);
  virtual ~Basis();
  void Save( const std::string filename, const std::string gname );
  void Load( const std::string filename, const std::string gname );
  inline bool GetType()const{return isFermion;};
  inline size_t GetL()const{return L;};
  inline size_t GetN()const{return N;};
  inline size_t GetHilbertSpace()const{
    if (isFermion) {
      return FStates.size();
    } else {
      return BStates.size();
    }
  };

  /* Fermion functions */
  void Fermion();
  void Fermion( const int IncludeAllU1 );// This has no U(1) symmetry.
  /* Spin - 1/2, share the Fermion functions */
  void SpinOneHalf();
  void TIsing();

  inline std::vector<int> GetFStates()const{
    return FStates;
  };//for Fermion and SpinOneHalf

  inline std::vector<size_t> GetFTags()const{
    return FTags;
  };//for Fermion and SpinOneHalf

  inline int GetNfTotal(const size_t idx)const{
    return NTotal.at(idx);
  };// for fermion without U(1)

  inline int GetSzTotal(const size_t idx)const{
    return NTotal.at(idx);
  };// for Transverse Ising spin

  inline void PrintFermionBasis( const int state )const{
    for (size_t cnt = 0; cnt < L; cnt++) {
      std::cout << btest(state, cnt) << ", " << std::flush;
    }
    std::cout << std::endl;
  };

  inline void PrintSpinOneHalfBasis( const int state )const{
    PrintFermionBasis( state );
  };

  /* Boson functions */
  void Boson();
  void Phonon();
  inline std::vector< std::vector<int> > GetBStates()const{
    return BStates;
  };

  inline std::vector<RealType> GetBTags()const{
    return BTags;
  };

  inline void PrintBosonBasis( const std::vector<int> state )const{
    typename std::vector<int>::const_iterator it = state.begin();
    for (;it != state.end(); ++it ){
      std::cout << *it << ", " << std::flush;
    }
    std::cout << std::endl;
  };

  inline void PrintAllBosonBasis()const{
    typename std::vector< std::vector<int> >::const_iterator it = BStates.begin();
    for ( ;it != BStates.end(); ++it ){
      PrintBosonBasis(*it);
    };
  };

  inline size_t GetIndexFromTag( const RealType tg )const{
    assert( !(isFermion) );
    auto it = binary_locate(BTags.begin(), BTags.end(), tg);
    size_t idx = std::distance(BTags.begin(), it);
    if ( std::abs(BTags.at(idx) - tg) < 1.0e-12 ){
      return idx;
    } else{
      return GetHilbertSpace() - 1;// To be a valid index
    }
  };

private:
  size_t L;
  size_t N;
  bool isFermion;
  bool HaveU1;
  std::vector< std::vector<int> > BStates;//for Boson
  std::vector<RealType> BTags;//for Boson
  std::vector<int> FStates;//for Fermion, and SpinOneHalf.
  std::vector<size_t> FTags;//for Fermion, and SpinOneHalf.
  std::vector<int> NTotal;//for Fermion without U(1).

  /* Holstein Phonon */
  std::vector<std::vector<int> > ApplyOffdiagonal( const std::vector<std::vector<int> > InputStates );
  bool CheckExist1( const std::vector<int>& State );
  bool CheckExist2( const std::vector<int>& State, const std::vector<std::vector<int> > States );
  void DummyCheckState();
};

RealType BosonBasisTag( const std::vector<int> vec );
template<typename T>
  std::vector<std::vector<int> > SortBTags( const std::vector<std::vector<int> >& st, std::vector<T> &v );

#endif//__BASIS_HPP__
