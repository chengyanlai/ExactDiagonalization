#ifndef __FERMI_HUBBARD_HPP__
#define __FERMI_HUBBARD_HPP__

#include "src/Hamiltonian/Hamiltonian.hpp"

template<typename Tnum>
class FHM: public Hamiltonian<Tnum> {
private:
  arma::SpMat<Tnum> H_Jup, H_Jdn;
  arma::SpMat<Tnum> H_Vup, H_Vdn;
  arma::SpMat<Tnum> H_U, H_W;
public:
  FHM ():Hamiltonian<Tnum>(){};
  FHM (const std::vector<Basis> &bs):Hamiltonian<Tnum>(bs){};
  virtual ~FHM (void){};

  void FermiHubbardModel( const std::vector<Basis>& bs, const std::vector< Node<Tnum>* >& lattice, const std::vector< std::vector<Tnum> >& Vloc, const std::vector<Tnum>& Uloc );
  void ExtendedFermiHubbardModel( const std::vector<Basis>& bs, const std::vector< Node<Tnum>* >& lattice, const std::vector< std::vector<Tnum> >& Vloc, const std::vector<Tnum>& Uloc, const std::vector<Tnum>& Wloc );
  void LocalPotential( const size_t species_id, const std::vector<Tnum> &Vloc, const Basis &bs );
  void ExtendedHubbardInteraction( const std::vector<Tnum> &Wloc, const std::vector<Basis> &bs);
  void HubbardInteraction( const std::vector<Tnum> &Uloc, const std::vector<Basis> &bs );
  void NNHopping( const size_t species_id, const std::vector< Node<Tnum>* > &lt, const Basis &bs );

  inline arma::SpMat<Tnum> GetHoppingHamiltonian(const size_t species_id){
    if( species_id == 0 ) return H_Jup;
    else if( species_id == 1 ) return H_Jdn;
  };
};

#endif	/* end of include guard: __FERMI_HUBBARD_HPP__ */
