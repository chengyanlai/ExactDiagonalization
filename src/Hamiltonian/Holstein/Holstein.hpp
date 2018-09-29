#ifndef __HOLSTEIN_HPP__
#define __HOLSTEIN_HPP__
#include <algorithm>//* find, distance

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
  Holstein ( const int Fi, const std::vector<Basis> &bs ):Hamiltonian<Tnum>(Fi, bs){};
  virtual ~Holstein (void){};

  void HolsteinModel( const std::vector<Basis>& bs, const RealType& k, const RealType& J, const RealType& Wloc, const RealType& Gloc );
  void PhononLocal( const RealType& Wloc, const Basis& bs );
  void FermionNNHopping( const RealType& k, const RealType& J, const Basis& bs, const RealType& Phi=0.0e0);
  void FermionPhononCoupling( const RealType& Gloc, const Basis& bs );

  void HolsteinModelR( const int& KPoints, const std::vector<Basis>& bs, const RealType& W, const RealType& G, const RealType& J );
  void PhononR( const int& KPoints, const RealType& W, const Basis& bs );
  void FermionR( const int& KPoints, const RealType& J, const Basis& bs );
  void FermionPhononR( const int& KPoints, const RealType& G, const Basis& bs);

  void HolsteinModelK( const std::vector<int>& KPoints, const std::vector<Basis>& bs, const RealType& W, const RealType& G, const RealType& J, const int& Nmax, const int WithoutK0Phonon );
  void PhononK( const Basis& bs, const RealType& W );
  void FermionK( const std::vector<int>& KPoints, const Basis& bs, const RealType& J, const RealType& Phi=0.0e0);
  void FermionPhononK( const std::vector<int>& KPoints, const Basis& bs, const RealType& G, const int& Nmax, const int WithoutK0Phonon);

  inline arma::SpMat<Tnum> GetHKinetic()const{
    return H_Kinetic;
  };
  inline arma::SpMat<Tnum> GetHPhonon()const{
    return H_Phonon;
  };
  inline arma::SpMat<Tnum> GetHCouple()const{
    return H_Couple;
  };
};

#endif	/* end of include guard: __HOLSTEIN_HPP__ */
