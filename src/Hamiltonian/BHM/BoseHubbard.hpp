#ifndef __BOSE_HUBBARD_HPP__
#define __BOSE_HUBBARD_HPP__

#include "src/Hamiltonian/Hamiltonian.hpp"

template<typename Tnum>
class BHM: public Hamiltonian<Tnum> {
private:
  arma::SpMat<Tnum> H_J;
  arma::SpMat<Tnum> H_V;
  arma::SpMat<Tnum> H_U;
public:
  BHM ():Hamiltonian<Tnum>(){};
  BHM (const std::vector<Basis> &bs):Hamiltonian<Tnum>(bs){};
  virtual ~BHM (void){};

  /* vvvvvvv Bose Hubbard model vvvvvvv */
  void BoseHubbardModel( const std::vector<Basis>& bs, const std::vector< Node<Tnum>* >& lattice, const std::vector<Tnum>& Vloc, const std::vector<Tnum>& Uloc );
  void LocalTerms( const std::vector<Tnum> &Vloc, const std::vector<Tnum> &Uloc, const Basis &bs );
  void NNHopping( const std::vector< Node<Tnum>* > &lt, const Basis &bs );
};

#endif	/* end of include guard: __BOSE_HUBBARD_HPP__ */
