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
  template<typename Tnum> friend class Hamiltonian;
  template<typename Tnum> friend class FHM;
  template<typename Tnum> friend class BHM;
  template<typename Tnum> friend class Holstein;
public:
  Basis(const bool _isFermion = false);
  Basis(const size_t _L, const size_t _N, const bool _isFermion = false);
  virtual ~Basis();
  void Save( const std::string filename, const std::string gname );
  void Load( const std::string filename, const std::string gname );

  inline bool GetType()const{
    return isFermion;
  };

  inline size_t GetL()const{
    return L;
  };

  inline size_t GetN()const{
    return N;
  };

  inline size_t GetHilbertSpace()const{
    if (isFermion) {
      return FStates.size();
    } else {
      return BStates.size();
    }
  };

  inline size_t size()const{
    return GetHilbertSpace();
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

  /* Boson functions */
  void Boson();
  void Phonon();
  inline std::vector< std::vector<int> > GetBStates()const{
    return BStates;
  };

  inline std::vector<RealType> GetBTags()const{
    return BTags;
  };

  inline size_t GetIndexFromTag( const RealType tg )const{
    assert( !(isFermion) );
    auto it = binary_locate(BTags.begin(), BTags.end(), tg);
    size_t idx = std::distance(BTags.begin(), it);
    // if ( std::abs(BTags.at(idx) - tg) < 1.0e-12 ){
      return idx;
    // } else{
    //   return GetHilbertSpace() - 1;// To be a valid index
    // }
  };

  /* Holstein Phonon */
  RealType CreatePhonon( std::vector<int>& state, const int site=0 )const;
  RealType DestroyPhonon( std::vector<int>& state, const int site=0 )const;
  RealType FermionJumpRight( std::vector<int>& state, const int NumJumps = 1 )const;
  RealType FermionJumpLeft( std::vector<int>& state, const int NumJumps = 1 )const;

  friend std::ostream& operator<<(std::ostream& os, const Basis& st);
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
  std::vector<std::vector<int> > ApplyOffdiagonal( const std::vector<std::vector<int> >& InputStates );
  void DummyCheckState();

};

RealType BosonBasisTag( const std::vector<int> vec );
template<typename T>
  std::vector<std::vector<int> > SortBTags( const std::vector<std::vector<int> >& st, std::vector<T> &v );

#endif//__BASIS_HPP__
