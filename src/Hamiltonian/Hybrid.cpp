#include <cassert>
#include <cmath>//sqrt
#include <map>
#include "src/bitwise.h"
#include "src/Hamiltonian/Hamiltonian.hpp"

template<typename Tnum>
void Hamiltonian<Tnum>::Hybridization( const int species1, const int species2, const std::vector< std::tuple<int, int, Tnum> > &hybVals, const std::vector<Basis> &bs, std::vector<Triplet> &hhyb, const int maxLocalB ){
  /* assuming only two species now */
  assert( bs.size() > 1 );
  Basis b1 = bs.at(species1);
  assert( !(b1.getType()) );
  Basis f2 = bs.at(species2);
  assert( f2.getType() );
  for ( int f : f2.FStates ){
    for ( typename std::vector< std::tuple<int, int, Tnum> >::const_iterator obj = hybVals.begin(); obj != hybVals.end(); ++obj ){
      int site1, site2;
      Tnum val;
      std::tie(site1, site2, val) = *obj;
      int bid1 = 0;
      for ( std::vector<int> b : b1.BStates ){
        size_t rid = 0, cid = 0;
        Tnum value = (Tnum)0.0e0;
        if ( btest(f, site2) && ( (b.at(site1) < maxLocalB && maxLocalB) || (b.at(site1) < b1.getN() && !(maxLocalB)) ) ){// hopping happens
          size_t fid1 = f2.FTags.at(f);
          rid = DetermineTotalIndex( vec<size_t>(bid1, fid1) );
          int CrossFermionNumber = 0;
          Tnum tsign = (Tnum)(1.0e0);
          if ( site2 > 0 ){
            for( int k = 0; k < site2; k++){
              CrossFermionNumber += btest(f,k);
            }
          }
          if (CrossFermionNumber % 2 == 1) tsign = (Tnum)(-1.0e0);
          // Find new basis for fermion
          size_t p = ibclr(f, site2);
          size_t fid2 = f2.FTags.at(p);
          // Find new basis for boson
          std::vector<int> q = b;
          q.at(site1) += 1;//creation
          size_t bid2 = b1.getIndexFromTag( BosonBasisTag(q) );
          if ( bid2 == b1.getHilbertSpace() ){
            continue;
          }
          cid = DetermineTotalIndex( vec<size_t>(bid2, fid2) );
          value = (Tnum)( val * sqrt( (Tnum)(q.at(site1)) ) ) * tsign;
          // std::cout << site1 << " " << site2 << ":" << rid << "(" << bid1 << "," << fid1 << ") " << cid << "(" << bid2 << "," << fid2 << ") " << value << std::endl;
        }else if ( !(btest(f, site2)) && b.at(site1) > 0 ){
          size_t fid1 = f2.FTags.at(f);
          rid = DetermineTotalIndex( vec<size_t>(bid1, fid1) );
          int CrossFermionNumber = 0;
          Tnum tsign = (Tnum)(1.0e0);
          if ( site2 > 0 ){
            for( int k = 0; k < site2; k++){
              CrossFermionNumber += btest(f,k);
            }
          }
          if (CrossFermionNumber % 2 == 1) tsign = (Tnum)(-1.0e0);
          // Find new basis for fermion
          size_t p = ibset(f, site2);
          size_t fid2 = f2.FTags.at(p);
          // Find new basis for boson
          std::vector<int> q = b;
          q.at(site1) -= 1;//annihilation
          size_t bid2 = b1.getIndexFromTag( BosonBasisTag(q) );
          if ( bid2 == b1.getHilbertSpace() ){
            continue;
          }
          cid = DetermineTotalIndex( vec<size_t>(bid2, fid2) );
          value = (Tnum)( val * sqrt( (Tnum)(b.at(site1)) ) ) * tsign;
        }
        if ( std::abs(value) > 1.0e-12 ){
          hhyb.push_back(Triplet(rid, cid, value));
        }
        bid1++;
      }
    }
  }
}

template class Hamiltonian<RealType>;
template class Hamiltonian<ComplexType>;
