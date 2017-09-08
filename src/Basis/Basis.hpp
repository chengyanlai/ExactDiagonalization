#ifndef __BASIS_HPP__
#define __BASIS_HPP__
#include <vector>
#include <cassert>
#include <algorithm>//sort, distance
#include "src/EDType.hpp"


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
  void Boson();
  void Boson( const int MinB, const int MaxB );
  void BosonTB( const size_t TBloc, const bool HARD_CUT = false );
  void Fermion();
  void Fermion( const int MinF );
  void SpinOneHalf();
  void TIsing();
  void Save( const std::string filename, const std::string gname );
  void Load( const std::string filename, const std::string gname );
  inline bool getType()const{return isFermion;};
  inline size_t getL()const{return L;};
  inline size_t getN()const{return N;};
  inline size_t getHilbertSpace()const{
    if (isFermion) {
      return FStates.size();
    } else {
      return BStates.size();
    }
  };
  inline std::vector<int> getFStates()const{return FStates;};//for Fermion and SpinOneHalf
  inline std::vector<size_t> getFTags()const{return FTags;};//for Fermion and SpinOneHalf
  inline int getSzTotal(const size_t idx)const{return SzTotal.at(idx);};//for Ising spin
  void printFermionBasis( const int state )const;
  void printSpinOneHalfBasis( const int state )const;
  inline std::vector< std::vector<int> > getBStates()const{return BStates;};//for Boson
  inline std::vector<RealType> getBTags()const{return BTags;};//for Boson
  void printBosonBasis( const std::vector<int> state )const;//for Boson
  inline size_t getIndexFromTag( const RealType tg )const{
    assert( !(isFermion) );
    auto it = binary_locate(BTags.begin(), BTags.end(), tg);
    size_t idx = std::distance(BTags.begin(), it);
    if ( std::abs(BTags.at(idx) - tg) < 1.0e-12 ){
      return idx;
    } else{
      return getHilbertSpace();
    }
  };//for Boson
private:
  size_t L;
  size_t N;
  bool isFermion;
  std::vector< std::vector<int> > BStates;//for Boson
  std::vector<RealType> BTags;//for Boson
  std::vector<int> FStates;//for Fermion, and SpinOneHalf.
  std::vector<size_t> FTags;//for Fermion, and SpinOneHalf.
  std::vector<int> NfTotal;//for Fermion without U(1).
  std::vector<int> SzTotal;//for Transverse Ising
};

RealType BosonBasisTag( const std::vector<int> vec );
template <typename T> std::vector<size_t> SortBTags( std::vector<T> &v );
#endif//__BASIS_HPP__
