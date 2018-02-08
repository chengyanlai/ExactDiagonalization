#include "src/Lindblad/ParticleExchange/lindblad.hpp"

/* NOTE: 4-th order Runge-Kutta integration.

d \Rhos = i[\Rhos, H] + \gamma ( Lmat1 - Lmat2 )
Lmat1 = L \Rhos L^\dagger
Lmat2 = 0.5 * {L^\dagger L, \Rhos}

f() is the Lindblad equation
x[i] is the Mat(\Rhos)

# k1 = dt * f( x[i], t )
# k2 = dt * f( x[i] + 0.5 * k1, t + 0.5 * dt )
# k3 = dt * f( x[i] + 0.5 * k2, t + 0.5 * dt )
# k4 = dt * f( x[i] + k3, t + dt )
# x[i+1] = x[i] + ( k1 + 2.0 * ( k2 + k3 ) + k4 ) / 6.0
*/

void FRK4( const RealType& dt, const std::vector<RealType>& gammas, const std::vector<std::tuple<int,int,int> >& SiteTypesSpin, const std::vector<std::vector<std::pair<int,int> > >& BasisIds, const std::vector<std::vector<std::vector<std::pair<size_t, size_t> > > >& CollapseIds, const std::vector<std::vector<Basis> >& bas, const std::vector<FHM<ComplexType> >& hams, std::vector<ComplexMatrixType>& Rhos ) {
  /* NOTE: This is only for quench */
  std::vector<ComplexMatrixType> k1 = Rhos;
  FNewton( dt, gammas, SiteTypesSpin, BasisIds, CollapseIds, bas, hams, k1);

  std::vector<ComplexMatrixType> k2;
  for (size_t cnt = 0; cnt < Rhos.size(); cnt++) {
    k2.push_back(Rhos.at(cnt) + 0.50e0* k1.at(cnt) );
  }
  FNewton( dt, gammas, SiteTypesSpin, BasisIds, CollapseIds, bas, hams, k2);

  std::vector<ComplexMatrixType> k3;
  for (size_t cnt = 0; cnt < Rhos.size(); cnt++) {
    k3.push_back(Rhos.at(cnt) + 0.50e0* k2.at(cnt) );
  }
  FNewton( dt, gammas, SiteTypesSpin, BasisIds, CollapseIds, bas, hams, k3);

  std::vector<ComplexMatrixType> k4;
  for (size_t cnt = 0; cnt < Rhos.size(); cnt++) {
    k4.push_back(Rhos.at(cnt) + k3.at(cnt) );
  }
  FNewton( dt, gammas, SiteTypesSpin, BasisIds, CollapseIds, bas, hams, k4);

  size_t cnt = 0;
  for ( auto &kf: Rhos ){
    kf += ( k1.at(cnt) + 2.0e0 * (k2.at(cnt) + k3.at(cnt)) + k4.at(cnt) ) / 6.0e0;
    cnt++;
  }
}

void FNewton( const RealType& dt, const std::vector<RealType>& gammas, const std::vector<std::tuple<int,int,int> >& SiteTypesSpin, const std::vector<std::vector<std::pair<int,int> > >& BasisIds, const std::vector<std::vector<std::vector<std::pair<size_t, size_t> > > >& CollapseIds, const std::vector<std::vector<Basis> >& bas, const std::vector<FHM<ComplexType> >& hams, std::vector<ComplexMatrixType>& Rhos ){
  assert( hams.size() == Rhos.size() );
  assert( gammas.size()  == SiteTypesSpin.size() );
  assert( SiteTypesSpin.size()  == BasisIds.size() );
  assert( BasisIds.size() == CollapseIds.size() );
  const ComplexType imagI = ComplexType(0.0e0, 1.0e0);
  std::vector<ComplexMatrixType> Commutator = Rhos;
  std::vector<ComplexMatrixType> LindbladMatrix;
  for ( size_t i = 0; i < hams.size(); i++){
    Commutator.at(i) = imagI * (Rhos.at(i) * hams.at(i).GetTotalHamiltonian() - hams.at(i).GetTotalHamiltonian() * Rhos.at(i));
    LindbladMatrix.push_back( ComplexMatrixType(Rhos.at(i).n_rows, Rhos.at(i).n_cols, arma::fill::zeros) );
  }
  for (size_t cnt = 0; cnt < gammas.size(); cnt++) {
    int site, type, spin;
    std::tie(site, type, spin) = SiteTypesSpin.at(cnt);
    RealType gamma = gammas.at(cnt);
    std::vector<std::pair<int, int> > BasisIdx = BasisIds.at(cnt);
    for ( size_t j = 0; j < BasisIdx.size(); j++ ){
      int NewBsIdx = BasisIdx.at(j).first;
      int OldBsIdx = BasisIdx.at(j).second;
      std::vector<std::pair<size_t, size_t> > CIdx = CollapseIds.at(cnt).at(j);
      LindbladMatrix.at(NewBsIdx) += gamma * LindbladTerm(LindbladMatrix.at(NewBsIdx).n_cols, CIdx, Rhos.at(OldBsIdx)) ;
    }
    for ( size_t i = 0; i < hams.size(); i++){
      LindbladMatrix.at(i) -= gamma * AntiCommutatorF( type, site, spin, bas.at(i), hams.at(i), Rhos.at(i));
    }
  }
  for ( size_t cnt = 0; cnt < Rhos.size(); cnt++ ){
    Rhos.at(cnt) = dt * ( Commutator.at(cnt) + LindbladMatrix.at(cnt) );
  }
}

ComplexMatrixType LindbladTerm( const size_t dim, const std::vector<std::pair<size_t, size_t> > &CIds, const ComplexMatrixType &MapMat) {
  // std::cout << "LB" << std::endl;
  ComplexMatrixType work(dim, dim, arma::fill::zeros);
  for ( auto row : CIds ){
    size_t newr = row.first;
    size_t oldr = row.second;
    for ( auto col : CIds ){
      size_t newc = col.first;
      size_t oldc = col.second;
      // std::cout << dim << " " << MapMat.rows() << " " << newr << " " << newc << " <- " << oldc << " " << newc << std::endl;
      work(newr,newc) = MapMat(oldr, oldc);
    }
  }
  // std::cout << "LB done" << std::endl;
  return work;
}

ComplexMatrixType AntiCommutatorF( const int type, const size_t Site, const size_t Spin, const std::vector<Basis> Bs, FHM<ComplexType> ham, const ComplexMatrixType &rho){
  /* NOTE: This returns the desity matrix multiply the particle number on Site
  respectively to each basis.
  type -  1 : c^\dagger c
         -1 : c c^\dagger
  */
  // std::cout << "AC" << std::endl;
  ComplexMatrixType tmp1 = rho;
  // std::cout << tmp1.rows() << " " << tmp1.cols() << " " << Bs.at(0).GetHilbertSpace() << " " <<  Bs.at(1).GetHilbertSpace() << std::endl;
  std::vector<int> f = Bs.at(Spin).GetFStates();
  int Spin2;
  if ( Spin ) Spin2 = 0;
  else Spin2 = 1;
  size_t coff = 0;
  for ( auto &fstate : f ){
    if ( (!(btest(fstate, Site)) && type == -1) || (btest(fstate, Site) && type == 1) ){
      // Without particle and destroy OR with particle and create.
      for ( size_t p=0; p < Bs.at(Spin2).GetHilbertSpace(); p++){
        size_t idx;
        if ( Spin ){
          std::vector<size_t> ids = vec<size_t>(p, coff);
          idx = ham.DetermineTotalIndex( ids );
        }else{
          std::vector<size_t> ids = vec<size_t>(coff, p);
          idx = ham.DetermineTotalIndex( ids );
        }
        tmp1.col(idx) *= (RealType)(0.0e0);
      }
    }
    coff++;
  }
  // std::cout << "AC done" << std::endl;
  return 0.50e0 * ( tmp1 + tmp1.t() );
}

/* This is Lindblad operator set to C which belongs to type -1
    d \rho ~ c \rho c^dagger
*/
void Cf(const size_t Site, const int Spin, const std::vector<std::vector<Basis> > &Bs, std::vector<FHM<ComplexType> > &hams, const std::map<std::pair<int,int>, int > PairIndex1, const std::map<int, std::pair<int,int> > PairIndex2, std::vector<std::pair<int,int> > &BasisIdx, std::vector<std::vector<std::pair<size_t, size_t> > > &CollapseIdx) {
  BasisIdx.clear();
  CollapseIdx.clear();
  int spin1, spin2;
  if ( Spin == 0 ){
    spin1 = 0;
    spin2 = 1;
  }else if ( Spin == 1 ){
    spin1 = 1;
    spin2 = 0;
  }
  for (size_t cnt = 0; cnt < Bs.size(); cnt++) {
    int Nup_Bfr = PairIndex2.at(cnt).first;
    int Ndn_Bfr = PairIndex2.at(cnt).second;
    int Nup_Afr, Ndn_Afr;
    // goes to a higher U(1) sector
    if ( Spin == 0 ){
      Nup_Afr = Nup_Bfr + 1;
      Ndn_Afr = Ndn_Bfr;
      /* NOTE: Hot fix */
      if ( Nup_Afr > Bs.at(cnt).at(Spin).GetL() ) continue;
    }else if ( Spin == 1 ){
      Nup_Afr = Nup_Bfr;
      Ndn_Afr = Ndn_Bfr + 1;
      /* NOTE: Hot fix */
      if ( Ndn_Afr > Bs.at(cnt).at(Spin).GetL() ) continue;
    }
    // std::cout << Nup_Bfr << " " << Ndn_Bfr << " -> " << Nup_Afr << " " << Ndn_Afr << std::endl;
    int BsIndex = PairIndex1.at(std::make_pair(Nup_Afr, Ndn_Afr));
    BasisIdx.push_back(std::make_pair(cnt, BsIndex));
    std::vector<std::pair<size_t, size_t> > tmp_idx;
    size_t oid1 = 0;
    for ( auto &ef: Bs.at(cnt).at(spin1).GetFStates() ) {
      if ( !(btest(ef, Site)) ){
        int nef = ibset(ef, Site);
        size_t id1 = Bs.at(BsIndex).at(spin1).GetFTags().at(nef);
        for ( size_t id2 = 0; id2 < Bs.at(cnt).at(spin2).GetHilbertSpace(); id2++ ){
          std::vector<size_t> oids, ids;
          if ( Spin == 0 ){
            oids.push_back(oid1);
            oids.push_back(id2);
            ids.push_back(id1);
            ids.push_back(id2);
          }else if ( Spin == 1 ){
            oids.push_back(id2);
            oids.push_back(oid1);
            ids.push_back(id2);
            ids.push_back(id1);
          }
          tmp_idx.push_back( std::make_pair(hams.at(cnt).DetermineTotalIndex(oids), hams.at(BsIndex).DetermineTotalIndex(ids)) );
        }
      }
      oid1 += 1;
    }
    CollapseIdx.push_back(tmp_idx);
  }
}

/* This is Lindblad operator set to C^\dagger which belongs to type 1
    d \rho ~ c^dagger \rho c
*/
void Cfdagger( const size_t Site, const int Spin, const std::vector<std::vector<Basis> > &Bs, std::vector<FHM<ComplexType> > &hams, const std::map<std::pair<int,int>, int > PairIndex1, const std::map<int, std::pair<int,int> > PairIndex2, std::vector<std::pair<int, int> > &BasisIdx, std::vector<std::vector<std::pair<size_t, size_t> > > &CollapseIdx) {
  BasisIdx.clear();
  CollapseIdx.clear();
  int spin1, spin2;
  if ( Spin == 0 ){
    spin1 = 0;
    spin2 = 1;
  }else if ( Spin == 1 ){
    spin1 = 1;
    spin2 = 0;
  }
  for (size_t cnt = 0; cnt < Bs.size(); cnt++) {
  // for (size_t cnt = Bs.size() - 1; cnt > 0; cnt--) {
    int Nup_Bfr = PairIndex2.at(cnt).first;
    int Ndn_Bfr = PairIndex2.at(cnt).second;
    int Nup_Afr, Ndn_Afr;
    if ( Spin == 0 && Nup_Bfr > 0 ){
      Nup_Afr = Nup_Bfr - 1;
      Ndn_Afr = Ndn_Bfr;
    }else if ( Spin == 1 && Ndn_Bfr > 0){
      Nup_Afr = Nup_Bfr;
      Ndn_Afr = Ndn_Bfr - 1;
    }else{
      continue;
    }
    // std::cout << Nup_Bfr << " " << Ndn_Bfr << " -> " << Nup_Afr << " " << Ndn_Afr << std::endl;
    int BsIndex = PairIndex1.at(std::make_pair(Nup_Afr, Ndn_Afr));
    // First is the current basis index and the later is the index collapse from.
    BasisIdx.push_back(std::make_pair(cnt, BsIndex));
    std::vector<std::pair<size_t, size_t> > tmp_idx;
    size_t oid1 = 0;
    for ( auto &ef: Bs.at(cnt).at(spin1).GetFStates() ) {
      if ( btest(ef, Site) ){
        int nef = ibclr(ef, Site);
        size_t id1 = Bs.at(BsIndex).at(spin1).GetFTags().at(nef);
        for ( size_t id2 = 0; id2 < Bs.at(cnt).at(spin2).GetHilbertSpace(); id2++ ){
          std::vector<size_t> oids, ids;
          if ( Spin == 0 ){
            oids.push_back(oid1);
            oids.push_back(id2);
            ids.push_back(id1);
            ids.push_back(id2);
          }else if ( Spin == 1 ){
            oids.push_back(id2);
            oids.push_back(oid1);
            ids.push_back(id2);
            ids.push_back(id1);
          }
          tmp_idx.push_back( std::make_pair(hams.at(cnt).DetermineTotalIndex(oids), hams.at(BsIndex).DetermineTotalIndex(ids)) );
        }
      }
      oid1 += 1;
    }
    CollapseIdx.push_back(tmp_idx);
  }
}
