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

  inline SparseMatrixType BuildSparseHamiltonian( const std::vector<std::tuple<int, int, Tnum> >& MatElemts ){
    ULongMatrixType Locations(2, MatElemts.size());
    VectorType Values( MatElemts.size() );
    typename std::vector<std::tuple<int, int, Tnum> >::const_iterator it = MatElemts.begin();
    size_t cnt = 0;
    for (; it != MatElemts.end(); ++it ){
      int row, col;
      Tnum val;
      std::tie(row, col, val) = *it;
      // armadillo
      Locations(0,cnt) = row;
      Locations(1,cnt) = col;
      Values(cnt) = val;
      // Eigen3
      cnt++;
    }
    // First true allows repeated matrix elements
    SparseMatrixType out(true, Locations, Values, GetTotalHilbertSpace(), GetTotalHilbertSpace());//, sort_locations = true, check_for_zeros = true);
    return out;
  };

  /* vvvvvvv Linear Algebra vvvvvvv */
  void eigh( RealVectorType &Vals, MatrixType &Vecs, const int nev=4, const bool randomInitial=true);
  void diag( RealVectorType &Vals, MatrixType &Vec);
  void expH( const ComplexType Prefactor, ComplexVectorType &Vec, const size_t Kmax = 15 );
  // RealVectorType expVals( const RealType Prefactor, const RealVectorType Vec);

protected:
  SparseMatrixType H_total;
  std::vector<size_t> HilbertSpaces;

private:
  void mvprod(Tnum* x, Tnum* y, RealType alpha) const;
  void arpackDiagonalize(int n, Tnum* input_ptr, std::vector<RealType> &evals, int nev = 1, RealType tol = 0.0e0);
};
#endif//__HAMILTONIAN_HPP__
