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
  Hamiltonian(){};
  Hamiltonian( const std::vector<Basis> &bs );
  virtual ~Hamiltonian(){};
  inline size_t GetTotalHilbertSpace()const{
    size_t tmp = 1;
    for (auto &j : HilbertSpaces){
      tmp *= j;
    }
    return tmp;
  };
  inline int CheckHermitian(){return arma::approx_equal(H_total, H_total.t(), "absdiff", 1.0e-5);};// armadillo
  void BuildTotalHamiltonian( const std::vector<std::tuple<int, int, Tnum> >& MatElemts );
  /* vvvvvvv Bose Hubbard model vvvvvvv */
  void BoseHubbardModel( const std::vector<Basis>& bs, const std::vector< Node<Tnum>* >& lattice, const std::vector<Tnum>& Vloc, const std::vector<Tnum>& Uloc );
  void LocalTerms( const std::vector<Tnum> &Vloc, const std::vector<Tnum> &Uloc, const Basis &bs, std::vector<std::tuple<int, int, Tnum> > &MatElemts );
  void NNHopping( const std::vector< Node<Tnum>* > &lt, const Basis &bs, std::vector<std::tuple<int, int, Tnum> > &MatElemts );
  /* ^^^^^^^ Bose Hubbard model ^^^^^^^ */
  /* vvvvvvv Fermi Hubbard model vvvvvvv */
  void FermiHubbardModel( const std::vector<Basis>& bs, const std::vector< Node<Tnum>* >& lattice, const std::vector< std::vector<Tnum> >& Vloc, const std::vector<Tnum>& Uloc );
  void LocalPotential( const size_t species_id, const std::vector<Tnum> &Vloc, const Basis &bs, std::vector<std::tuple<int, int, Tnum> > &MatElemts );
  void HubbardInteraction( const std::vector<int> species_id, const std::vector<Tnum> &Uloc, const std::vector<Basis> &bs, std::vector<std::tuple<int, int, Tnum> > &MatElemts );
  void NNHopping( const size_t species_id, const std::vector< Node<Tnum>* > &lt, const Basis &bs, std::vector<std::tuple<int, int, Tnum> > &MatElemts );
  /* ^^^^^^^  Fermi Hubbard model ^^^^^^^ */
  // void FermionIntraNN( const int speciesId, const std::vector<std::tuple<int, int, Tnum> > betweenSitesVals, const Basis &bs, std::vector<std::tuple<int, int, Tnum> > &MatElemts );
  // void FermionDensityDensity(
  //   const std::vector<std::pair<int,int> > betweenSpecies, const std::vector<std::tuple<int, int, Tnum> > betweenSitesVals,
  //   const std::vector<Basis> &bs, std::vector<Triplet> &hloc );
  /* vvvvvvv Spin Functions vvvvvvv */
  // void SpinOneHalfXXZ( const Tnum Delta, const std::vector< Node<Tnum>* > &lt, const Basis &bs, std::vector<Triplet> &hhop);
  // void TIsing( const Tnum Jz, const std::vector< Node<Tnum>* > &lt, const Basis &bs, std::vector<Triplet> &hhop );
  /* ^^^^^^^ Spin Functions ^^^^^^^ */
  /* vvvvvvv Hybrid systems vvvvvvv */
  // void BuildHybridHamiltonian( const int species1, const int species2, const std::vector< std::tuple<int, int, Tnum> > &hybVals, const std::vector<Basis> &bs, const int maxLocalB = 0);
  // void Hybridization( const int species1, const int species2, const std::vector< std::tuple<int, int, Tnum> > &hybVals, const std::vector<Basis> &bs, std::vector<Triplet> &hhyb, const int maxLocalB);
  /* ^^^^^^^ Hybrid systems ^^^^^^^ */
  void eigh( RealVectorType &Vals, MatrixType &Vecs, const int nev=4, const bool randomInitial=true);
  void diag( RealVectorType &Vals, MatrixType &Vec);
  void expH( const ComplexType Prefactor, ComplexVectorType &Vec, const size_t Kmax = 15 );
  RealVectorType expVals( const RealType Prefactor, const RealVectorType Vec);
  void mvprod(Tnum* x, Tnum* y, RealType alpha) const;
  inline SparseMatrixType GetTotalHamiltonian()const{return H_total;};
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
private:
  std::vector<size_t> HilbertSpaces;
  SparseMatrixType H_total;
  void arpackDiagonalize(int n, Tnum* input_ptr, std::vector<RealType> &evals, int nev = 1, RealType tol = 0.0e0);
};
#endif//__HAMILTONIAN_HPP__
