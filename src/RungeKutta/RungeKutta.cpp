#include "src/RungeKutta/RungeKutta.hpp"

/* NOTE: 4-th order Runge-Kutta integration.

d \Rhos = i[\Rhos, H] + \gamma ( Lmat1 - Lmat2 )
Lmat1 = a \Rhos a^\dagger
Lmat2 = 0.5 * {a^\dagger a, \Rhos}

f() is the Lindblad equation
x[i] is the Mat(\Rhos)

# k1 = dt * f( x[i], t )
# k2 = dt * f( x[i] + 0.5 * k1, t + 0.5 * dt )
# k3 = dt * f( x[i] + 0.5 * k2, t + 0.5 * dt )
# k4 = dt * f( x[i] + k3, t + dt )
# x[i+1] = x[i] + ( k1 + 2.0 * ( k2 + k3 ) + k4 ) / 6.0
*/
void Lindblad_RK4( const RealType &dt, const RealType &gamma,
  const size_t TBloc, const std::vector<Basis> &bas,
  const std::vector<Hamiltonian<ComplexType> > &ham,
  const std::vector<std::vector<size_t> > &CIdx,
  std::vector<ComplexMatrixType> &Rhos) {
  std::vector<ComplexMatrixType> k1 = Rhos;
  Lindblad_Newton( dt, gamma, TBloc, bas, ham, CIdx, k1);

  std::vector<ComplexMatrixType> k2;
  for (size_t cnt = 0; cnt < Rhos.size(); cnt++) {
    k2.push_back(Rhos.at(cnt) + 0.50e0* k1.at(cnt) );
  }
  Lindblad_Newton( dt, gamma, TBloc, bas, ham, CIdx, k2);

  std::vector<ComplexMatrixType> k3;
  for (size_t cnt = 0; cnt < Rhos.size(); cnt++) {
    k3.push_back(Rhos.at(cnt) + 0.50e0* k2.at(cnt) );
  }
  Lindblad_Newton( dt, gamma, TBloc, bas, ham, CIdx, k3);

  std::vector<ComplexMatrixType> k4;
  for (size_t cnt = 0; cnt < Rhos.size(); cnt++) {
    k4.push_back(Rhos.at(cnt) + k3.at(cnt) );
  }
  Lindblad_Newton( dt, gamma, TBloc, bas, ham, CIdx, k4);

  size_t cnt = 0;
  for ( auto &kf: Rhos ){
    kf += ( k1.at(cnt) + 2.0e0 * (k2.at(cnt) + k3.at(cnt)) + k4.at(cnt) ) / 6.0e0;
    cnt++;
  }
}

void Lindblad_Newton( const RealType &dt, const RealType &gamma,
  const size_t TBloc, const std::vector<Basis> &bas,
  const std::vector<Hamiltonian<ComplexType> > &ham,
  const std::vector<std::vector<size_t> > &CIdx,
  std::vector<ComplexMatrixType> &Rhos) {
  assert( ham.size() == Rhos.size() );
  assert( ham.size() == CIdx.size() + 1);
  const ComplexType imagI = ComplexType(0.0e0, 1.0e0);
  for (size_t cnt = 0; cnt < ham.size(); cnt++) {
    ComplexSparseMatrixType h = ham.at(cnt).getTotalHamiltonian();
    ComplexMatrixType LBmat1 = Nb_tbloc(TBloc, bas.at(cnt), Rhos.at(cnt));
    if ( cnt + 1 < ham.size() ) {
      assert( CIdx.at(cnt).size() == Rhos.at(cnt).cols() );
      ComplexMatrixType LBmat2 = Lindblad1(TBloc, Rhos.at(cnt+1), bas.at(cnt), CIdx.at(cnt));
      Rhos.at(cnt) = dt * ( imagI * (Rhos.at(cnt) * h - h * Rhos.at(cnt)) -
                            gamma * LBmat1 + gamma * LBmat2 );
    } else {
      Rhos.at(cnt) = dt * ( imagI * (Rhos.at(cnt) * h - h * Rhos.at(cnt)) -
                            gamma * LBmat1 );
    }
  }
}

ComplexMatrixType Lindblad1(const size_t TBloc, const ComplexMatrixType &MapMat,
  const Basis &bs, const std::vector<size_t> &CIdx){
  std::vector< std::vector<int> > b = bs.getBStates();
  assert( b.size() == CIdx.size() );
  ComplexMatrixType work = ComplexMatrixType::Zero(CIdx.size(), CIdx.size());
  for (size_t row = 0; row < CIdx.size(); row++) {
    for (size_t col = 0; col < CIdx.size(); col++) {
      RealType val = 0.0e0;
      if ( row == col ) {
        val = (RealType)(b.at(row).at(TBloc) + 1);
      } else {
        RealType val1 = (RealType)(b.at(row).at(TBloc) + 1);
        RealType val2 = (RealType)(b.at(col).at(TBloc) + 1);
        val = std::sqrt(val1*val2);
      }
      work(row,col) = val * MapMat(CIdx.at(row),CIdx.at(col));
    }
  }
  return work;
}

ComplexMatrixType Nb_tbloc( const size_t TBloc, const Basis &bs,
  const ComplexMatrixType &rho){
  ComplexMatrixType tmp1 = rho;
  std::vector< std::vector<int> > b = bs.getBStates();
  assert( b.size() == rho.cols() );
  assert( b.size() == rho.rows() );
  size_t coff = 0;
  for ( auto &nbi : b ){
    tmp1.col(coff) *= (RealType)nbi.at(TBloc);
    coff++;
  }
  return 0.50e0 * ( tmp1 + tmp1.adjoint() );
}
