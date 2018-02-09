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
  void FermionNNHopping( const RealType& k, const std::vector< Node<Tnum>* > &lt, const Basis& bs );

  void HolsteinModel( const std::vector<Basis>& bs, const RealType& k, const RealType& J, const std::vector<Tnum>& Wloc, const std::vector<Tnum>& Gloc );
  void FermionNNHoppingInfinite( const RealType& k, const RealType& J, const Basis& bs, const RealType& Phi=0.0e0);

  void FermionPhononCoupling( const std::vector<Tnum>& Gloc, const Basis& bs );
  void PhononLocal( const std::vector<Tnum>& Wloc, const Basis& bs );

  inline arma::SpMat<Tnum> GetHKinetic()const{
    return H_Kinetic;
  };
};

#endif	/* end of include guard: __HOLSTEIN_HPP__ */
