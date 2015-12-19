#include <cassert>
#include <cmath>//sqrt
#include <map>
#include "src/Hamiltonian/Hamiltonian.hpp"

template<typename Tnum>
void Hamiltonian<Tnum>::BosonIntraLocalPart( const size_t species_id,
  const std::vector<Tnum> &Vloc, const std::vector<Tnum> &Uloc,
  const Basis &bs, std::vector<Triplet> &hloc )
{
  int state_id = 0;
  for ( std::vector<int> b : bs.BStates ){
    int loc = 0;
    Tnum val = (Tnum)0.0;
    for ( int &p : b ){
      val += Vloc.at(loc) * (Tnum)p +
        (Tnum)0.50e0 * Uloc.at(loc) * (Tnum)(p * (p - 1));
      loc++;
    }
    //FIXME: Need to loop over all species_id!
    //use while
    // std::vector<size_t> idrun(HilbertSpaces.size(), 0);
    // idrun.at(species_id) = HilbertSpaces.at(HilbertSpaces);
    // bool DONE = false;
    // while ( !DONE ){
    //   for (size_t cnt = 0; cnt < HilbertSpaces.size(); cnt++) {
    //   }
    // }
    size_t idx = DetermineTotalIndex( std::vector<size_t>(1, state_id) );
    hloc.push_back(Triplet(idx, idx, val));
    state_id++;
  }
}

template<typename Tnum>
void Hamiltonian<Tnum>::BosonIntraHoppingPart( const size_t species_id,
const std::vector< Node<Tnum>* > &lt,
const Basis &bs, std::vector<Triplet> &hhop )
{
  int state_id1 = 0;
  for ( std::vector<int> b : bs.BStates ){
    for ( Node<Tnum>* l : lt ) {
      size_t site_i = l->data;
      std::vector< Node<Tnum>* > nn = l->getNeighbors();
      std::vector< Tnum > nnJ = l->getJval();
      for (size_t cnt = 0; cnt < l->NumNeighbors(); cnt++) {
        std::vector<int> b_copy = b;
        size_t site_j = nn.at(cnt)->data;
        if ( b_copy.at(site_i) > 0 ){
          Tnum val = (Tnum)(-1.0e0) * nnJ.at(cnt) *
            sqrt( (Tnum)((b.at(site_j) + 1) * b.at(site_i)) );
          b_copy.at(site_j) = b.at(site_j) + 1;
          b_copy.at(site_i) = b.at(site_i) - 1;
          size_t state_id2 = bs.getIndexFromTag( BosonBasisTag(b_copy) );
          if ( state_id2 < bs.getHilbertSpace() ){
            //FIXME: Need to loop over all species_id!
            // size_t idx1 = DetermineTotalIndex(species_id, state_id1);
            size_t idx1 = DetermineTotalIndex(std::vector<size_t>(1, state_id1));
            // size_t idx2 = DetermineTotalIndex(species_id, state_id2);
            size_t idx2 = DetermineTotalIndex(std::vector<size_t>(1, state_id2));
            hhop.push_back(Triplet(idx2, idx1, val));
            // INFO(idx2 << "(" << BosonBasisTag(b_copy) << " " << bs.BTags.at(idx2) << ") " << idx1 << " " << val);
          }
        }
      }
    }
    state_id1++;
  }
}

template class Hamiltonian<RealType>;
template class Hamiltonian<ComplexType>;
