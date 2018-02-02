#include <cassert>
#include <cmath>//sqrt
#include <map>
#include "src/bitwise.h"
#include "src/Hamiltonian/Hamiltonian.hpp"

template<typename Tnum>
void Hamiltonian<Tnum>::SpinOneHalfXXZ( const Tnum Delta,
  const std::vector< Node<Tnum>* > &lt, const Basis &bs, std::vector<Triplet> &hhop )
{
  size_t rid, cid;
  size_t bid = 0, pid = 0;//l and p's index
  for ( int b : bs.FStates ){
    for ( Node<Tnum>* l : lt ) {
      size_t site_i = l->data;
      std::vector< Node<Tnum>* > nn = l->getNeighbors();
      std::vector< Tnum > Jxy = l->getJval();
      for (size_t cnt = 0; cnt < l->NumNeighbors(); cnt++) {
        size_t site_j = nn.at(cnt)->data;// spin flipping from site-j to site-i
        /* see if flipping exist */
        if ( btest(b, site_j) && !(btest(b, site_i)) ) {
          /* if yes, no particle in i and one particle in j. */
          size_t p = ibset(b,site_i);
          p = ibclr(p,site_j);
          bid = bs.FTags.at(b);// Find their indices
          pid = bs.FTags.at(p);// Find their indices
          Tnum value = (Tnum)(0.50e0 * Jxy.at(cnt) - 0.250e0 * Delta * Jxy.at(cnt));
          hhop.push_back(Triplet(bid, pid, value));
        } else{
          // they have the same spin
          bid = bs.FTags.at(b);
          Tnum value = (Tnum)(0.250e0 * Delta * Jxy.at(cnt));
          hhop.push_back(Triplet(bid, bid, value));
        }
      }
    }
  }
}

template<typename Tnum>
void Hamiltonian<Tnum>::TIsing( const Tnum hz,
  const std::vector< Node<Tnum>* > &lt, const Basis &bs, std::vector<Triplet> &hhop )
{
  size_t bid = 0, pid = 0;//l and p's index
  for ( int b : bs.FStates ){
    bid = bs.FTags.at(b);// Find their indices
    Tnum valueLocal = (Tnum)(-1.0e0 * hz * (Tnum)bs.getSzTotal(b));
    hhop.push_back(Triplet(bid, bid, valueLocal));
    for ( Node<Tnum>* l : lt ) {
      size_t site_i = l->data;
      std::vector< Node<Tnum>* > nn = l->getNeighbors();
      std::vector< Tnum > Jxx = l->getJval();
      for (size_t cnt = 0; cnt < l->NumNeighbors(); cnt++) {
        size_t site_j = nn.at(cnt)->data;
        int p = b;
        if ( btest(b, site_i) ){
          p = ibclr(p,site_i);
        }else{
          p = ibset(p,site_i);
        }
        if ( btest(b, site_j) ){
          p = ibclr(p,site_j);
        }else{
          p = ibset(p,site_j);
        }
        pid = bs.FTags.at(p);// Find their indices
        Tnum value = (Tnum)(-0.1250e0 * Jxx.at(cnt));
        hhop.push_back(Triplet(bid, pid, value));
      }
    }
  }
}

template class Hamiltonian<RealType>;
template class Hamiltonian<ComplexType>;
