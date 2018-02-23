#include <cassert>
#include <cmath>//sqrt
#include <map>
#include <tuple>
#include "src/bitwise.h"
#include "src/Hamiltonian/FHM/FermiHubbard.hpp"

template<typename Tnum>
void FHM<Tnum>::LocalPotential( const size_t spin1, const std::vector<Tnum> &Vloc, const Basis &bs){
  std::vector<std::tuple<int, int, Tnum> > MatElemts;
  MatElemts.clear();
  int spin2 = (spin1 == 0)?  1 : 0;
  int state_id = 0;
  for ( int b : bs.FStates ){
    Tnum val = (Tnum)0.0;
    for (size_t loc = 0; loc < bs.GetL(); loc++) {
      if ( btest(b, loc) ){
        val += Vloc.at(loc);
      }
    }
    if ( std::abs(val) > 1.0e-12 ){
      std::vector<size_t> ids(2, state_id);
      for (size_t id2 = 0; id2 < this->HilbertSpaces.at(spin2); id2++) {
        ids.at(spin2) = id2;
        size_t idx = this->DetermineTotalIndex( ids );
        MatElemts.push_back(std::make_tuple(idx, idx, val));
      }
    }
    state_id++;
  }
  if (spin1 == 0) {
    H_Vup = this->BuildSparseHamiltonian( MatElemts );
  } else if (spin1 == 1){
    H_Vdn = this->BuildSparseHamiltonian( MatElemts );
  } else{
    RUNTIME_ERROR("Not support more than 2 species fermion yet!");
  }
}

template<typename Tnum>
void FHM<Tnum>::HubbardInteraction( const std::vector<Tnum> &Uloc, const std::vector<Basis> &bs){
  std::vector<std::tuple<int, int, Tnum> > MatElemts;
  MatElemts.clear();
  /*NOTE: Calculate interaction between any two input bases.
          Assume two input bases has the same L. */
  assert( bs.size() == 2 );
  assert( bs.at(0).GetL() == bs.at(1).GetL() );
  std::vector< std::vector<size_t> > IndexU1;
  std::vector< std::vector<size_t> > IndexU2;
  size_t point1, point2;
  /*NOTE: Build the IndexU1 and IndexU2.
          IndexU has all index of species-1 which has occupied particle at site-i
  */
  for (size_t i = 0; i < bs.at(0).GetL(); i++) {
    std::vector<size_t> work;
    point1 = 0;// point1 will be the same for all i
    for (size_t cnt_up = 0; cnt_up < this->HilbertSpaces.at(0); cnt_up++) {
      if ( btest(bs.at(0).FStates.at(cnt_up), i) ) {
        point1++;
        work.push_back(cnt_up);
      }
    }
    IndexU1.push_back(work);
    work.clear();
  }
  if ( bs.at(0).GetN() == bs.at(1).GetN() ) {
    point2 = point1;
    IndexU2 = IndexU1;
  }
  else{
    for (size_t i = 0; i < bs.at(1).GetL(); i++) {
      std::vector<size_t> work;
      point2 = 0;// point2 will be the same for all i
      for (size_t cnt_dn = 0; cnt_dn < this->HilbertSpaces.at(1); cnt_dn++) {
        if ( btest(bs.at(1).FStates.at(cnt_dn), i) ) {
          point2++;
          work.push_back(cnt_dn);
        }
      }
      IndexU2.push_back(work);
      work.clear();
    }
  }
  // Calculate the interactions
  for (size_t i = 0; i < bs.at(0).GetL(); i++) {
    for (size_t up = 0; up < point1; up++) {
      for (size_t dn = 0; dn < point2; dn++) {
        std::vector<size_t> ids(2,0);
        ids.at(0) = IndexU1.at(i).at(up);
        ids.at(1) = IndexU2.at(i).at(dn);
        size_t id = this->DetermineTotalIndex( ids );
        MatElemts.push_back( std::make_tuple(id, id, Uloc.at(i)) );
        // std::cout << "2 " << id << " " << id << " " << Uloc << std::endl;
      }
    }
  }
  H_U = this->BuildSparseHamiltonian( MatElemts );
}

template<typename Tnum>
void FHM<Tnum>::NNHopping( const size_t spin1, const std::vector< Node<Tnum>* > &lt, const Basis &bs){
  std::vector<std::tuple<int, int, Tnum> > MatElemts;
  MatElemts.clear();
  int spin2 = (spin1 == 0)?  1 : 0;
  size_t bid = 0, pid = 0;//l and p's index
  for ( int b : bs.FStates ){
    for ( Node<Tnum>* l : lt ) {
      int site_i = l->data;
      std::vector< Node<Tnum>* > nn = l->GetNeighbors();
      std::vector< Tnum > nnJ = l->GetJval();
      for (size_t cnt = 0; cnt < l->NumNeighbors(); cnt++) {
        int site_j = nn.at(cnt)->data;// hopping from site-j to site-i
        /* see if hopping exist */
        if ( btest(b, site_j) && !(btest(b, site_i)) ) {
          /* if yes, no particle in i and one particle in j. */
          int CrossFermionNumber = 0;
          Tnum tsign = (Tnum)(1.0e0);
          if ( std::labs(site_i - site_j) > 1 ){
            // possible cross fermions and sign change.
            if ( site_i > site_j ){
              for ( int k = site_j+1; k < site_i; k++){
                CrossFermionNumber += btest(b, k);
              }
            }else if ( site_i < site_j ){
              for ( int k = site_i+1; k < site_j; k++){
                CrossFermionNumber += btest(b, k);
              }
            }else{
              std::cout << "Something wrong in hopping sign......" << std::endl;
            }
            // std::cout << "Sign change " << site_i << " " << site_j << " " << CrossFermionNumber << std::endl;
          }
          if (CrossFermionNumber % 2 == 1) tsign = (Tnum)(-1.0e0);
          size_t p = ibset(b,site_i);
          p = ibclr(p,site_j);
          // std::cout << site_i << " " << site_j << " " << b << " " << p << " " << bs.FTags.size() << std::endl;
          // bs.printFermionBasis(b);
          // std::cout << "->" << std::flush;
          // bs.printFermionBasis(p);
          // std::cout << "." << std::endl;
          bid = bs.FTags.at(b);// Find their indices
          pid = bs.FTags.at(p);// Find their indices
          // INFO(lid << " " << pid);
          std::vector<size_t> rids(2, bid);
          std::vector<size_t> cids(2, pid);
          for (size_t loop_id = 0; loop_id < this->HilbertSpaces.at(spin2); loop_id++) {
            rids.at(spin2) = loop_id;
            cids.at(spin2) = loop_id;
            size_t rid = this->DetermineTotalIndex( rids );
            size_t cid = this->DetermineTotalIndex( cids );
            Tnum value = (Tnum)(-1.0e0) * tsign * nnJ.at(cnt);
            MatElemts.push_back( std::make_tuple(rid, cid, value) );
          }
        }
      }
    }
  }
  if (spin1 == 0) {
    H_Jup = this->BuildSparseHamiltonian( MatElemts );
  } else if (spin1 == 1){
    H_Jdn = this->BuildSparseHamiltonian( MatElemts );
  }
}

template<typename Tnum>
void FHM<Tnum>::FermiHubbardModel( const std::vector<Basis>& bs, const std::vector< Node<Tnum>* >& lattice, const std::vector< std::vector<Tnum> >& Vloc, const std::vector<Tnum>& Uloc ){
  // Both species live in the same lattice and the same happing amplitudes
  assert( bs.size() == 2 );//NOTE: Only support two or one species right now.
  assert( Vloc.size() == bs.size() );
  int cnt = 0;
  for ( auto &b : bs ){
    /* For intra-species local terms: Potential */
    assert( b.GetL() == Vloc.at(cnt).size() );
    assert( b.GetL() == Uloc.size() );
    LocalPotential( cnt, Vloc.at(cnt), b );
    /* For intra-species N-N hopping */
    assert( b.GetL() == lattice.size() );
    NNHopping( cnt, lattice, b );
    cnt++;
  }
  /* For inter-species local terms: Hubbard U */
  HubbardInteraction(Uloc, bs);
  /* Build H_total */
  this->H_total = H_Jup + H_Jdn + H_Vup + H_Vdn + H_U;
}

// template<typename Tnum>// Extended Hubbard
// void Hamiltonian<Tnum>::FermionIntraNN( const int speciesId, const std::vector<std::tuple<int, int, Tnum> > betweenSitesVals, const Basis &bs, std::vector<std::tuple<int, int, Tnum> > &MatElemts ){
//   size_t count;
//   if ( speciesId == 0 ) {
//     count = HilbertSpaces.at(1);
//   } else {
//     count = HilbertSpaces.at(0);
//   }
//   int stateId = 0;
//   for ( int b : bs.FStates ){
//     Tnum FinalVal = (Tnum)(0.0e0);
//     // for ( auto obj : betweenSitesVals){
//     for ( typename std::vector<std::tuple<int, int, Tnum> >::const_iterator obj=betweenSitesVals.begin(); obj != betweenSitesVals.end(); ++obj ){
//       int site1, site2;
//       Tnum val;
//       // std::tie(site1, site2, val) = obj;
//       std::tie(site1, site2, val) = *obj;
//       if ( btest(b, site1) && btest(b, site2) ){
//         FinalVal += val;
//       }
//     }
//     if ( std::abs(FinalVal) > 1.0e-12 ){
//       std::vector<size_t> ids(2,stateId);
//       for ( size_t p=0; p < count; p++){
//         if ( speciesId == 0 ) {
//           ids.at(1) = p;
//         } else {
//           ids.at(0) = p;
//         }
//         size_t id = DetermineTotalIndex( ids );
//         MatElemts.push_back( std::make_tuple(id, id, FinalVal) );
//       }
//     }
//     stateId++;
//   }
// }

template class FHM<RealType>;
template class FHM<ComplexType>;
