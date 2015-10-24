#include <cassert>
#include <cmath>//sqrt
#include <map>
#include "src/bitwise.h"
#include "src/Hamiltonian/Hamiltonian.hpp"

template<typename Tnum, typename Tlabel>
void Hamiltonian<Tnum, Tlabel>::FermionIntraLocalPart( const size_t species_id,
  const std::vector<Tnum> &Vloc, const Basis &bs, std::vector<Triplet> &hloc )
{
  int state_id = 0;
  for ( int b : bs.FStates ){
    Tnum val = (Tnum)0.0;
    for (size_t loc = 0; loc < bs.getL(); loc++) {
      if ( btest(b, loc) ){
        val += Vloc.at(loc);
      }
    }
    std::vector<size_t> ids(2, state_id);
    if (species_id == 0) {
      for (size_t id2 = 0; id2 < HilbertSpaces.at(1); id2++) {
        ids.at(1) = id2;
        size_t idx = DetermineTotalIndex( ids );
        hloc.push_back(Triplet(idx, idx, val));
      }
    } else if (species_id == 1){
      for (size_t id2 = 0; id2 < HilbertSpaces.at(0); id2++) {
        ids.at(0) = id2;
        size_t idx = DetermineTotalIndex( ids );
        hloc.push_back(Triplet(idx, idx, val));
      }
    } else{
      RUNTIME_ERROR("Not support more than 2 species fermion yet!");
    }
    state_id++;
  }
}

template<typename Tnum, typename Tlabel>
void Hamiltonian<Tnum, Tlabel>::FermionInterLocalPart(
  const std::vector<int> species_id, const Tnum &Uloc,
  const std::vector<Basis> &bs, std::vector<Triplet> &hloc )
{
  /*NOTE: Calculate interaction between any two input bases.
  Assume two input bases has the same L.
  Assume the interaction is uniform in space.*/
  assert( bs.size() == 2 );
  assert( bs.at(0).getL() == bs.at(1).getL() );
  std::vector< std::vector<size_t> > IndexU1;
  std::vector< std::vector<size_t> > IndexU2;
  size_t point1, point2;
  int TotalN = bs.at(0).getN() + bs.at(1).getN();
  /*NOTE: Build the IndexU1 and IndexU2. */
  for (size_t i = 0; i < bs.at(0).getL(); i++) {
    std::vector<size_t> work;
    point1 = 0;
    for (size_t cnt_up = 0; cnt_up < HilbertSpaces.at(species_id.at(0)); cnt_up++) {
      if ( TotalN <= bs.at(0).getL() ) {
        if ( btest(bs.at(0).FStates.at(cnt_up), i) ) {
          point1++;
          work.push_back(cnt_up);
        }
      }
      else{
        if ( !(btest(bs.at(0).FStates.at(cnt_up), i)) ) {
          point1++;
          work.push_back(cnt_up);
        }
      }
    }
    IndexU1.push_back(work);
    work.clear();
  }
  if ( bs.at(0).getN() == bs.at(1).getN() ) {
    point2 = point1;
    IndexU2 = IndexU1;
  }
  else{
    for (size_t i = 0; i < bs.at(0).getL(); i++) {
      std::vector<size_t> work;
      point2 = 0;
      for (size_t cnt_dn = 0; cnt_dn < HilbertSpaces.at(species_id.at(1)); cnt_dn++) {
        if ( TotalN <= bs.at(1).getL() ) {
          if ( btest(bs.at(1).FStates.at(cnt_dn), i) ) {
            point2++;
            work.push_back(cnt_dn);
          }
        }
        else{
          if ( !(btest(bs.at(1).FStates.at(cnt_dn), i)) ) {
            point2++;
            work.push_back(cnt_dn);
          }
        }
      }
      IndexU2.push_back(work);
      work.clear();
    }
  }
  /*NOTE: Add the overall interaction. */
  Tnum g = (Tnum)0.0e0;
  if ( TotalN > bs.at(1).getL() ) {
    g = Uloc * (Tnum)( TotalN - bs.at(0).getL() );
    for (size_t cnt = 0; cnt < getTotalHilbertSpace(); cnt++) {
      hloc.push_back(Triplet(cnt, cnt, g));
    }
  }
  if ( bs.at(0).getN() == bs.at(1).getN() ) {
    if ( TotalN <= bs.at(0).getL() ) {
      g = Uloc * (Tnum)(bs.at(0).getN());
    }
    else{
      g = Uloc * (Tnum)( bs.at(0).getL() - bs.at(0).getN() );
    }
    for (size_t cnt_up = 0; cnt_up < HilbertSpaces.at(species_id.at(1)); cnt_up++) {
      size_t id = DetermineTotalIndex( std::vector<size_t>(2, cnt_up) );
      hloc.push_back(Triplet(id, id, g));
    }
  }
  /*NOTE: for detatil difference*/
  for (size_t i = 0; i < bs.at(0).getL(); i++) {
    for (size_t up = 0; up < point1; up++) {
      for (size_t dn = 0; dn < point2; dn++) {
        if ( bs.at(0).getN() == bs.at(1).getN() ) {
          if ( up != dn ) {
            std::vector<size_t> ids(2,0);
            ids.at(0) = IndexU1.at(i).at(up);
            ids.at(1) = IndexU2.at(i).at(dn);
            size_t id = DetermineTotalIndex( ids );
            hloc.push_back(Triplet(id, id, Uloc));
          }
        }
        else {
          std::vector<size_t> ids(2,0);
          ids.at(0) = IndexU1.at(i).at(up);
          ids.at(1) = IndexU2.at(i).at(dn);
          size_t id = DetermineTotalIndex( ids );
          hloc.push_back(Triplet(id, id, Uloc));
        }
      }
    }
  }
}

template<typename Tnum, typename Tlabel>
void Hamiltonian<Tnum, Tlabel>::FermionIntraHoppingPart( const size_t species_id,
  const std::vector< Node<Tnum, Tlabel>* > &lt,
  const Basis &bs, std::vector<Triplet> &hhop )
{
  size_t rid, cid;
  size_t bid = 0, pid = 0;//l and p's index
  for ( int b : bs.FStates ){
    for ( Node<Tnum, Tlabel>* l : lt ) {
      Tlabel site_i = l->data;
      std::vector< Node<Tnum, Tlabel>* > nn = l->getNeighbors();
      std::vector< Tnum > nnJ = l->getJval();
      for (size_t cnt = 0; cnt < l->NumNeighbors(); cnt++) {
        Tlabel site_j = nn.at(cnt)->data;// hopping from site-j to site-i
        /* see if hopping exist */
        if ( btest(b, site_j) && !(btest(b, site_i)) ) {
          /* if yes, no particle in i and one particle in j. */
          size_t p = ibset(b,site_i);
          p = ibclr(p,site_j);
          bid = bs.FTags.at(b);
          pid = bs.FTags.at(p);
          // INFO(lid << " " << pid);
          size_t count;
          if ( species_id == 0 ) {
            count = HilbertSpaces.at(1);
          } else {
            count = HilbertSpaces.at(0);
          }
          std::vector<size_t> rids(2, bid);
          std::vector<size_t> cids(2, pid);
          for (size_t loop_id = 0; loop_id < count; loop_id++) {
            if ( species_id == 0 ) {
              rids.at(1) = loop_id;
              cids.at(1) = loop_id;
            } else {
              rids.at(0) = loop_id;
              cids.at(0) = loop_id;
            }
            rid = DetermineTotalIndex( rids );
            cid = DetermineTotalIndex( cids );
            Tnum value = (Tnum)(-1.0e0) * nnJ.at(cnt);
            hhop.push_back(Triplet(rid, cid, value));
            // hhop.push_back(Triplet(cid, rid, value));
            // INFO( rid << " " <<  cid << " " <<  value );
          }
        }
      }
    }
  }
}

template class Hamiltonian<RealType, int>;
template class Hamiltonian<ComplexType, int>;
