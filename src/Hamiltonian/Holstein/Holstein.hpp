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
  Holstein ( const int& KPoints, const std::vector<Basis> &bs ):Hamiltonian<Tnum>(KPoints, bs){};
  virtual ~Holstein (void){};

  void HolsteinModel( const std::vector<Basis>& bs, const RealType& k, const RealType& J, const std::vector<Tnum>& Wloc, const std::vector<Tnum>& Gloc );
  void PhononLocal( const std::vector<Tnum>& Wloc, const Basis& bs );
  void FermionNNHopping( const RealType& k, const RealType& J, const Basis& bs, const RealType& Phi=0.0e0);
  void FermionPhononCoupling( const std::vector<Tnum>& Gloc, const Basis& bs );

  void HolsteinModelK( const std::vector<int>& KPoints, const std::vector<Basis>& bs, const RealType& W, const RealType& G, const RealType& J, const int& Nmax );
  void PhononK( const int& KPoints, const RealType& W, const Basis& bs );
  void FermionK( const int& KPoints, const RealType& J, const Basis& bs, const RealType& Phi=0.0e0);
  void FermionPhononK( const std::vector<int>& KPoints, const RealType& G, const Basis& bs, const int& Nmax);

  inline arma::SpMat<Tnum> GetHKinetic()const{
    return H_Kinetic;
  };
  inline arma::SpMat<Tnum> GetHPhonon()const{
    return H_Phonon;
  };
  inline arma::SpMat<Tnum> GetHCouple()const{
    return H_Couple;
  };

  inline int KnToIndex( const int Kn )const{
    if ( Kn < 0 ) return 2 * std::abs(Kn);
    else if ( Kn > 0 ) return 2 * Kn - 1;
    else return Kn;
  };

  inline int IndexToKn( const int Index )const{
    if ( Index % 2 == 0 ){
      return -1 * Index / 2;
    }else{
      return (Index + 1) / 2;
    }
  };

  inline void FirstBZ( const int& L,  int& Kn )const{
    assert( L % 2 == 0 );
    while ( Kn > 0 && Kn > L/2 ) Kn -= L;
    while ( Kn < 0 && Kn <= -L/2 ) Kn += L;
    if ( Kn == - L / 2 ) Kn = L / 2;
  };
};

#endif	/* end of include guard: __HOLSTEIN_HPP__ */
