#ifndef __HOLSTEIN_HPP__
#define __HOLSTEIN_HPP__

#include "src/Hamiltonian/Hamiltonian.hpp"

template<typename Tnum=ComplexType>
class Holstein: public Hamiltonian<Tnum> {
private:
  arma::SpMat<Tnum> H_Kinetic;
  arma::SpMat<Tnum> H_Phonon;
  arma::SpMat<Tnum> H_Couple;
public:
  Holstein ():Hamiltonian<Tnum>(){};
  Holstein ( const std::vector<Basis> &bs ):Hamiltonian<Tnum>(bs){};
  virtual ~Holstein (void){};

  void HolsteinModel( const std::vector<Basis>& bs, const RealType& k, const std::vector< Node<Tnum>* >& lattice, const std::vector<Tnum>& Wloc, const std::vector<Tnum>& Gloc );
  void PhononLocal( const std::vector<Tnum>& Wloc, const Basis& bs );
  void FermionPhononCoupling( const std::vector<Tnum>& Gloc, const Basis& bs );
  void FermionNNHopping( const RealType& k, const std::vector< Node<Tnum>* > &lt, const Basis& bs );
};

#endif	/* end of include guard: __HOLSTEIN_HPP__ */
