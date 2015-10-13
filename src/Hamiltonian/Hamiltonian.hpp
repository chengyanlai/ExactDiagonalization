#ifndef __HAMILTONIAN_HPP__
#define __HAMILTONIAN_HPP__
#include <vector>
#include "src/EDType.hpp"
#include "src/EigenMatrix.hpp"
#include "src/Node/Node.hpp"
#include "src/Basis/Basis.hpp"

template<typename Tnum = RealType, typename Tlabel = int>
class Hamiltonian{
public:
  typedef Eigen::Matrix<Tnum, Eigen::Dynamic, Eigen::Dynamic,
    Eigen::AutoAlign|Eigen::RowMajor> MatrixType;
  typedef Eigen::SparseMatrix<Tnum, Eigen::AutoAlign|Eigen::RowMajor>
    SparseMatrixType;
  typedef Eigen::Matrix<Tnum, Eigen::Dynamic, 1, Eigen::AutoAlign> VectorType;
  typedef Eigen::Triplet<Tnum> Triplet;
  Hamiltonian( const std::vector<Basis> &bs );
  virtual ~Hamiltonian();
  inline size_t getTotalHilbertSpace()const{
    size_t tmp = 1;
    for (auto &j : HilbertSpaces){
      tmp *= j;
    }
    return tmp;
  };
  inline void BuildTotalHamiltonian(){H_total = H_hop + H_local;};
  void BuildIntraLocalHamiltonian(
    const std::vector< std::vector<Tnum> > &Vloc,
    const std::vector< std::vector<Tnum> > &Uloc,
    const std::vector<Basis> &bs );
  void BuildIntraHoppingHamiltonian(
    const std::vector<Basis> &bs, const std::vector< Node<Tnum, Tlabel>* > &lt );
  void BosonIntraLocalPart( const size_t species_id,
    const std::vector<Tnum> &Vloc, const std::vector<Tnum> &Uloc,
    const Basis &bs, std::vector<Triplet> &hloc );
  void BosonIntraHoppingPart( const size_t species_id,
    const std::vector< Node<Tnum, Tlabel>* > &lt,
    const Basis &bs, std::vector<Triplet> &hhop );
  void eigh( const bool FullDiagonalization = false );
  // void exp( const bool FullDiagonalization = false );
  // inline SparseMatrixType getTotalHamiltonian()const{return H_total;};
private:
  std::vector<size_t> HilbertSpaces;
  SparseMatrixType H_hop;
  SparseMatrixType H_local;
  SparseMatrixType H_total;
  inline size_t DetermineTotalIndex( size_t species_idx, size_t state_idx){
    size_t tidx = 0;
    for (size_t cnt = 0; cnt < HilbertSpaces.size(); cnt++) {
      if ( species_idx > cnt ){
        tidx += (species_idx + 1) * HilbertSpaces.at(cnt);
      }else{
        assert( state_idx < HilbertSpaces.at(cnt) );
        tidx += state_idx;
        break;
      }
    }
    return tidx;
  };
};
#endif//__HAMILTONIAN_HPP__
