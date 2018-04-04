#include "src/Hamiltonian/Hamiltonian.hpp"
#include "src/Lanczos/krylov.hpp"
#include "src/numeric/lapack.h"

template<typename Tnum>
void Hamiltonian<Tnum>::eigh( RealVectorType &Vals, MatrixType &Vecs, const int nev, const bool randomInitial, const std::string Target)const{
  size_t dim = GetTotalHilbertSpace();
  if ( randomInitial ) Vecs = MatrixType(dim, nev, arma::fill::randu);
  Tnum* input_ptr = Vecs.memptr();
  std::vector<RealType> Val;
  arpackDiagonalize(dim, input_ptr, Val, nev, /*tol*/0.0e0, Target);
  Vecs = MatrixType(input_ptr, dim, nev);
  Vals = RealVectorType(Val);
}

template<typename Tnum>
void Hamiltonian<Tnum>::diag( RealVectorType &Vals, MatrixType &Vecs)const{
  size_t dim = GetTotalHilbertSpace();
  // convert H_total to dense matrix
  MatrixType Mat(H_total);
  // working space
  Tnum* EigVec = (Tnum*)malloc( dim * dim * sizeof(Tnum) );
  RealType* Eig = (RealType*)malloc( dim * sizeof(RealType) );
  syDiag(Mat.memptr(), dim, Eig, EigVec);
  Vecs = MatrixType(EigVec, dim, dim);
  Vals = RealVectorType(Eig, dim);
}

template<>
void Hamiltonian<ComplexType>::expH( const ComplexType Prefactor, ComplexVectorType& Vec, const size_t Kmax )const{
  krylovEXP(H_total, Vec, Prefactor, Kmax);
}

template<>
void Hamiltonian<ComplexType>::HKrylov( RealVectorType &Vals, ComplexMatrixType &Vecs, const ComplexVectorType& Vec, const size_t Kmax )const{
  krylov(H_total, Vec, Vecs, Vals, Kmax);
}

/* Matrix vector product with MomHamiltonian: y = H_total * x + alpha * y
 * @param x the input vector
 * @param y the output vector
 * @param alpha the scaling value
 */
template<typename Tnum>
void Hamiltonian<Tnum>::mvprod(Tnum* x, Tnum* y, RealType alpha)const {
  size_t dim = GetTotalHilbertSpace();
  VectorType Vin(x, dim, true, false);//, copy_aux_mem = true, strict = false)// armadillo
  VectorType Vout(y, dim, true, false);//, copy_aux_mem = true, strict = false)// armadillo
  Vout = H_total * Vin + alpha * Vout;
  memcpy(y, Vout.memptr(), dim * sizeof(Tnum) );// armadillo
}

template class Hamiltonian<RealType>;
template class Hamiltonian<ComplexType>;


// template<typename Tnum>
// void Hamiltonian<Tnum>::BuildXXZHamiltonian(const Tnum Delta,
//   const std::vector<Basis> &bs, const std::vector< Node<Tnum>* > &lt )
// {
//   // std::vector<Triplet> hhop;// Eigen3
//   // hhop.clear();// Eigen3
//   for ( auto &b : bs ){
//     // SpinOneHalfXXZ( Delta, lt, b, hhop );// Eigen3
//   }
//   // H_hop.setFromTriplets(hhop.begin(), hhop.end());// Eigen3
//   // if (DEBUG) std::cout << "Non-zero matrix elements = " << hhop.size() << std::endl;// Eigen3
// }

// template<typename Tnum>
// void Hamiltonian<Tnum>::BuildTIsingHamiltonian(const Tnum hz,
//   const std::vector<Basis> &bs, const std::vector< Node<Tnum>* > &lt )
// {
//   // std::vector<Triplet> hhop;// Eigen3
//   // hhop.clear();// Eigen3
//   for ( auto &b : bs ){
//     // TIsing( hz, lt, b, hhop );// Eigen3
//   }
//   // H_hop.setFromTriplets(hhop.begin(), hhop.end());// Eigen3
//   // if (DEBUG) std::cout << "Non-zero matrix elements = " << hhop.size() << std::endl;// Eigen3
// }

// template<typename Tnum>
// void Hamiltonian<Tnum>::BuildHybridHamiltonian( const int species1, const int species2, const std::vector< std::tuple<int, int, Tnum> > &hybVals, const std::vector<Basis> &bs, const int maxLocalB ){
//   // std::vector<Triplet> hhyb;// Eigen3
//   // Hybridization( species1, species2, hybVals, bs, hhyb, maxLocalB );// Eigen3
//   // H_hybridization.setFromTriplets(hhyb.begin(), hhyb.end());// Eigen3
//   // if (DEBUG) std::cout << "Non-zero matrix elements = " << hhyb.size() << " from H_hyb!" << std::endl;// Eigen3
// }
