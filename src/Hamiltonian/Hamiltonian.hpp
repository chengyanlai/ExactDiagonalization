#ifndef __HAMILTONIAN_HPP__
#define __HAMILTONIAN_HPP__
#include <vector>
#include <utility>
#include <tuple>
#include "src/EDType.hpp"
#include "src/ArmadilloMatrix.hpp"
#include "src/Node/Node.hpp"
#include "src/Basis/Basis.hpp"

template<typename Tnum = RealType>
class Hamiltonian{
public:
  typedef arma::Col<Tnum> VectorType;
  typedef arma::Mat<Tnum> MatrixType;
  typedef arma::SpMat<Tnum> SparseMatrixType;

  /* Constructors */
  Hamiltonian(){};

  Hamiltonian( const std::vector<Basis> &bs ){
    HilbertSpaces.clear();
    for ( auto &b : bs ){
      HilbertSpaces.push_back(b.GetHilbertSpace());
    }
    size_t TotalDim = GetTotalHilbertSpace();
  };

  Hamiltonian( const int& Fi, const std::vector<Basis> &bs ){
    HilbertSpaces.clear();
    HilbertSpaces.push_back(Fi);
    for ( auto &b : bs ){
      HilbertSpaces.push_back(b.GetHilbertSpace());
    }
    size_t TotalDim = GetTotalHilbertSpace();
  };

  virtual ~Hamiltonian(){};

  inline size_t GetTotalHilbertSpace()const{
    size_t tmp = 1;
    for (auto &j : HilbertSpaces){
      tmp *= j;
    }
    return tmp;
  };

  inline int CheckHermitian(){
    return arma::approx_equal(H_total, H_total.t(), "absdiff", 1.0e-5);
  };

  inline SparseMatrixType GetTotalHamiltonian()const{
    return H_total;
  };

  inline void ShiftEnergy(const RealType& Esf){
    SparseMatrixType EDiag = Esf * arma::speye<SparseMatrixType>(GetTotalHilbertSpace(), GetTotalHilbertSpace());
    H_total += EDiag;
  };

  inline size_t DetermineTotalIndex( const std::vector<size_t> ids )const{
    assert( ids.size() == HilbertSpaces.size() );
    size_t tidx = 0;
    size_t factor = 1;
    for (size_t cnt = 0; cnt < ids.size(); cnt++) {
      tidx += ids.at(cnt) * factor;
      factor *= HilbertSpaces.at(cnt);
    }
    return tidx;
  };

  /* vvvvvvv Linear Algebra vvvvvvv */
  void eigh( RealVectorType &Vals, MatrixType &Vecs, const int nev=4, const bool randomInitial=true, const std::string Target="SR")const;
  void diag( RealVectorType &Vals, MatrixType &Vec)const;
  void expH( const ComplexType Prefactor, ComplexVectorType &Vec, const size_t Kmax = 20 )const;
  void HKrylov( RealVectorType &Vals, ComplexMatrixType &Vecs, const ComplexVectorType& Vec, const size_t Kmax )const;
  // RealVectorType expVals( const RealType Prefactor, const RealVectorType Vec);

protected:
  SparseMatrixType H_total;
  std::vector<size_t> HilbertSpaces;

private:
  void mvprod(Tnum* x, Tnum* y, RealType alpha) const;
  void arpackDiagonalize(int n, Tnum* input_ptr, std::vector<RealType> &evals, int nev = 1, RealType tol = 0.0e0, std::string Target="SR")const;
};
#endif//__HAMILTONIAN_HPP__
