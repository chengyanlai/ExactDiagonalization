#include "src/Lindblad-SS/BEC.hpp"

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
void Lindblad_RK4( const RealType &dt,
	const std::vector< std::vector< std::tuple<ptrdiff_t,RealType> > > Idx1,
	const std::vector< std::tuple<int,int> > Idx2,
	const std::vector<RealType> Gammas, const std::vector<Basis> &bas,
	const Hamiltonian<ComplexType> &ham, ComplexMatrixType &Rhos){
		/* NOTE: This is only for quench */
	ComplexMatrixType k1 = Rhos;
	Lindblad_Newton( dt, Idx1, Idx2, Gammas, bas, ham, k1);

	ComplexMatrixType k2 = ( Rhos + 0.50e0 * k1 );
	Lindblad_Newton( dt, Idx1, Idx2, Gammas, bas, ham, k2);

	ComplexMatrixType k3 = ( Rhos + 0.50e0 * k2 );
	Lindblad_Newton( dt, Idx1, Idx2, Gammas, bas, ham, k3);

	ComplexMatrixType k4 = ( Rhos + k3 );
	Lindblad_Newton( dt, Idx1, Idx2, Gammas, bas, ham, k4);

	Rhos += ( k1 + 2.0e0 * (k2 + k3) + k4 ) / 6.0e0;
}

void Lindblad_Newton( const RealType &dt,
	const std::vector< std::vector< std::tuple<ptrdiff_t,RealType> > > Idx1,
	const std::vector< std::tuple<int,int> > Idx2,
	const std::vector<RealType> Gammas,
	const std::vector<Basis> &bas,
	const Hamiltonian<ComplexType> &ham,
	ComplexMatrixType &Rhos) {
	const ComplexType imagI = ComplexType(0.0e0, 1.0e0);
	ComplexSparseMatrixType h = ham.getTotalHamiltonian();
	ComplexMatrixType LBmat1 = Lindblad1(Idx1, Gammas, bas.at(0), Rhos);
	ComplexMatrixType LBmat2 = Lindblad2(Idx2, Gammas, bas.at(0), Rhos);
	ComplexMatrixType Commutator = Rhos * h - h * Rhos;
	// std::cout << "Commatator " << Commutator.trace() << std::endl;
	std::cout.precision(18);
	std::cout << "LBmat1 " << LBmat1.trace() << std::endl;
	std::cout << "LBmat2 " << LBmat2.trace() << std::endl;
	Rhos = dt * ( imagI * Commutator +	LBmat1 - LBmat2 );
}

ComplexMatrixType Lindblad1( const std::vector< std::vector< std::tuple<ptrdiff_t,RealType> > > Idx,
	const std::vector<RealType> Gammas,
	const Basis &bs, const ComplexMatrixType &rho) {
	int id1, id2;
	RealType prefactor1, prefactor2;
	ComplexMatrixType work = ComplexMatrixType::Zero(rho.rows(), rho.cols());
	std::vector< std::vector<int> > b = bs.getBStates();
	assert( b.size() == rho.cols() );
	assert( b.size() == rho.rows() );
	assert( Gammas.size() == Idx.size() );
	for (size_t cnt = 0; cnt < Gammas.size(); cnt++) {
		RealType gamma = Gammas.at(cnt);
		assert( b.size() == Idx.at(cnt).size() );
		size_t dim = b.size();
		for (size_t row = 0; row < dim; row++) {
			std::tie(id1,prefactor1) = Idx.at(cnt).at(row);
			if ( id1 < 0 ) continue;
			for (size_t col = 0; col < dim; col++) {
				std::tie(id2,prefactor2) = Idx.at(cnt).at(col);
				if ( id2 < 0 ) continue;
				work(row,col) += gamma * prefactor1 * prefactor2 * rho(id1,id2);
			}
		}
	}
	return 0.50e0 * ( work + work.adjoint() );
}

ComplexMatrixType Lindblad2( const std::vector< std::tuple<int,int> > Idx,
	const std::vector<RealType> Gammas,
	const Basis &bs, const ComplexMatrixType &rho) {
	/* NOTE: This implements L_n = \gamm_n c^\dagger_i c_j.
					 (i,j) pair and their corresponding gamma are stored in the tuples. */
	int idxi, idxj;
	ComplexMatrixType tmp1 = rho;
	std::vector< std::vector<int> > b = bs.getBStates();
	assert( b.size() == rho.cols() );
	assert( b.size() == rho.rows() );
	assert( Idx.size() == Gammas.size() );
	size_t coff = 0;
	for ( auto &nbi : b ){
		RealType val = 0.0e0;
		for (size_t cnt = 0; cnt < Idx.size(); cnt++) {
			std::tie(idxi,idxj) = Idx.at(cnt);
			RealType gamma = Gammas.at(cnt);
			if ( idxi == idxj ) {
				val += gamma * (RealType)(nbi.at(idxi) * nbi.at(idxj));
			} else {
				val += gamma * (RealType)(nbi.at(idxj) * (1 + nbi.at(idxi)));
			}
		}
		tmp1.col(coff) *= val;
		coff++;
	}
	return 0.50e0 * ( tmp1 + tmp1.adjoint() );
}
