#include "src/Dephasing/Dephasing.hpp"

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
  const std::vector<size_t> &sites, const std::vector<Basis> &bas,
  const Hamiltonian<ComplexType> &ham, ComplexMatrixType &Rhos){
  ComplexMatrixType k1 = Rhos;
  Lindblad_Newton( dt, gamma, sites, bas, ham, k1);

  ComplexMatrixType k2 = (Rhos + 0.50e0* k1 );
  Lindblad_Newton( dt, gamma, sites, bas, ham, k2);

  ComplexMatrixType k3 = (Rhos + 0.50e0* k2 );
  Lindblad_Newton( dt, gamma, sites, bas, ham, k3);

  ComplexMatrixType k4 = (Rhos + k3 );
  Lindblad_Newton( dt, gamma, sites, bas, ham, k4);

  Rhos += ( k1 + 2.0e0 * (k2 + k3) + k4 ) / 6.0e0;
}

void Lindblad_Newton( const RealType &dt, const RealType &gamma,
  const std::vector<size_t> &sites, const std::vector<Basis> &bas,
  const Hamiltonian<ComplexType> &ham, ComplexMatrixType &Rhos) {
  const ComplexType imagI = ComplexType(0.0e0, 1.0e0);
  ComplexSparseMatrixType h = ham.getTotalHamiltonian();
  ComplexMatrixType LBmat1 = Lindblad1(sites, bas.at(0), Rhos);
  ComplexMatrixType LBmat2 = Lindblad2(sites, bas.at(0), Rhos);
  Rhos = dt * ( imagI * (Rhos * h - h * Rhos) +
                        gamma * LBmat1 - gamma * LBmat2 );
}

ComplexMatrixType Lindblad1( const std::vector<size_t> &sites, const Basis &bs,
  const ComplexMatrixType &rho){
  std::vector< std::vector<int> > b = bs.getBStates();
  ComplexMatrixType tmp = ComplexMatrixType::Zero(rho.rows(), rho.cols());
  assert( b.size() == rho.cols() );
  assert( b.size() == rho.rows() );
  for ( auto &site : sites){
    size_t coff = 0;
    ComplexMatrixType work = rho;
    for ( auto &nbi : b ){
      RealType val = (RealType)nbi.at(site);
      work.col(coff) *= val;
      work.row(coff) *= val;
      coff++;
    }
    tmp += work;
  }
  return tmp;
}

ComplexMatrixType Lindblad2( const std::vector<size_t> &sites, const Basis &bs,
  const ComplexMatrixType &rho){
  ComplexMatrixType tmp = rho;
  std::vector< std::vector<int> > b = bs.getBStates();
  assert( b.size() == rho.cols() );
  assert( b.size() == rho.rows() );
  size_t coff = 0;
  for ( auto &nbi : b ){
    RealType val = 0.0e0;
    for ( auto &site : sites){
      RealType nb = (RealType)nbi.at(site);
      val += nb * nb;
    }
    tmp.col(coff) *= val;
    coff++;
  }
  return 0.50e0 * ( tmp + tmp.adjoint() );
}
