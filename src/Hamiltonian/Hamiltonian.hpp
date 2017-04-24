#ifndef __HAMILTONIAN_HPP__
#define __HAMILTONIAN_HPP__
#include <vector>
#include "src/EDType.hpp"
#include "src/EigenMatrix.hpp"
#include "src/Node/Node.hpp"
#include "src/Basis/Basis.hpp"

template<typename Tnum = RealType>
class Hamiltonian{
public:
  typedef Eigen::Matrix<Tnum, Eigen::Dynamic, Eigen::Dynamic,
    Eigen::AutoAlign|Eigen::RowMajor> MatrixType;
  typedef Eigen::SparseMatrix<Tnum, Eigen::AutoAlign|Eigen::RowMajor>
    SparseMatrixType;
  typedef Eigen::Matrix<Tnum, Eigen::Dynamic, 1, Eigen::AutoAlign> VectorType;
  typedef Eigen::Triplet<Tnum> Triplet;
  typedef Eigen::Map<MatrixType> MapMatrix;
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
  // inline void CheckTotalHamiltonian(){H_total.isApprox(H_total)};
  void BuildLocalHamiltonian(
    const std::vector< std::vector<Tnum> > &Vloc,
    const std::vector< std::vector<Tnum> > &Uloc,
    const std::vector<Basis> &bs );
  void BuildHoppingHamiltonian(
    const std::vector<Basis> &bs, const std::vector< Node<Tnum>* > &lt );
  void BuildHoppingHamiltonian(
    const std::vector<Basis> &bs, const std::vector< std::vector< Node<Tnum>* > > &lt );
  void BuildXXZHamiltonian(const Tnum Delta, const std::vector<Basis> &bs,
    const std::vector< Node<Tnum>* > &lt );
  void BuildTIsingHamiltonian(const Tnum hz, const std::vector<Basis> &bs,
    const std::vector< Node<Tnum>* > &lt );
  /* vvvvvvv Boson Functions vvvvvvv */
  void BosonIntraLocalPart( const size_t species_id,
    const std::vector<Tnum> &Vloc, const std::vector<Tnum> &Uloc,
    const Basis &bs, std::vector<Triplet> &hloc );
  void BosonIntraHoppingPart( const size_t species_id,
    const std::vector< Node<Tnum>* > &lt,
    const Basis &bs, std::vector<Triplet> &hhop );
  /* ^^^^^^^ Boson Functions ^^^^^^^ */
  /* vvvvvvv Fermion Functions vvvvvvv */
  void FermionIntraLocalPart( const size_t species_id,
    const std::vector<Tnum> &Vloc, const Basis &bs, std::vector<Triplet> &hloc );
  void FermionInterLocalPart( const std::vector<int> species_id, const Tnum &Uloc,
      const std::vector<Basis> &bs, std::vector<Triplet> &hloc );
  void FermionInterLocalPart( const std::vector<int> species_id, const std::vector<Tnum> &Uloc,
      const std::vector<Basis> &bs, std::vector<Triplet> &hloc );
  void FermionIntraHoppingPart( const size_t species_id,
    const std::vector< Node<Tnum>* > &lt,
    const Basis &bs, std::vector<Triplet> &hhop );
  /* ^^^^^^^ Fermion Functions ^^^^^^^ */
  /* vvvvvvv Spin Functions vvvvvvv */
  void SpinOneHalfXXZ( const Tnum Delta, const std::vector< Node<Tnum>* > &lt,
    const Basis &bs, std::vector<Triplet> &hhop);
  void TIsing( const Tnum Jz, const std::vector< Node<Tnum>* > &lt,
    const Basis &bs, std::vector<Triplet> &hhop );
  /* ^^^^^^^ Spin Functions ^^^^^^^ */
  void eigh( std::vector<RealType> &Val, VectorType &Vec);
  void expH( const ComplexType Prefactor, ComplexVectorType &Vec, const size_t Kmax = 15 );
  void diag( RealVectorType &Val, RealMatrixType &Vec);
  void diag( RealVectorType &Val, ComplexMatrixType &Vec);
  RealVectorType expVals( const RealType Prefactor, const RealVectorType Vec);
  void mvprod(Tnum* x, Tnum* y, RealType alpha) const;
  inline SparseMatrixType getTotalHamiltonian()const{return H_total;};
  inline size_t DetermineTotalIndex( const std::vector<size_t> ids ){
    assert( ids.size() == HilbertSpaces.size() );
    size_t tidx = 0;
    size_t factor = 1;
    for (size_t cnt = 0; cnt < ids.size(); cnt++) {
      tidx += ids.at(cnt) * factor;
      factor *= HilbertSpaces.at(cnt);
    }
    return tidx;
  };
private:
  std::vector<size_t> HilbertSpaces;
  SparseMatrixType H_hop;
  SparseMatrixType H_local;
  SparseMatrixType H_total;
  void arpackDiagonalize(int n, Tnum* input_ptr, std::vector<RealType> &evals,
    int nev = 1, RealType tol = 0.0e0);
};
#endif//__HAMILTONIAN_HPP__
