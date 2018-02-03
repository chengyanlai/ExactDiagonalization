#include "src/Basis/Basis.hpp"
#include "src/hdf5io/hdf5io.hpp"

Basis::Basis(const bool _isFermion):
isFermion(_isFermion){
  BStates.clear();
  BTags.clear();
  FStates.clear();
  FTags.clear();
}

Basis::Basis(const size_t _L, const size_t _N, const bool _isFermion):
L(_L), N(_N), isFermion(_isFermion){
  BStates.clear();
  BTags.clear();
  FStates.clear();
  FTags.clear();
}

Basis::~Basis(){
  BStates.clear();
  BTags.clear();
  FStates.clear();
  FTags.clear();
}

void Basis::Save( const std::string filename, const std::string gname ){
  HDF5IO file(filename);
  if ( isFermion ) {
    file.SaveStdVector(gname, "Tags", FTags);
    file.SaveStdVector(gname, "State", FStates);
  }else {
    file.SaveStdVector(gname, "Tags", BTags);
    size_t cnt = 0;
    std::string sub_gname = gname;
    sub_gname.append("/Fock/");
    for ( auto &st : BStates){
      std::string dset_name = "State-";
      dset_name.append(std::to_string((unsigned long long)cnt));
      file.SaveStdVector(sub_gname, dset_name, st);
      cnt++;
    }
  }
}

void Basis::Load( const std::string filename, const std::string gname ){
  HDF5IO file(filename);
  if ( isFermion ) {
    file.LoadStdVector(gname, "Tags", FTags);
    file.LoadStdVector(gname, "State", FStates);
  } else {
    file.LoadStdVector(gname, "Tags", BTags);
    BStates.clear();
    std::string sub_gname = gname;
    sub_gname.append("/Fock/");
    for (size_t cnt = 0; cnt < BTags.size(); cnt++) {
      std::string dset_name = "State-";
      dset_name.append(std::to_string((unsigned long long)cnt));
      std::vector<int> tmp;
      file.LoadStdVector(sub_gname, dset_name, tmp);
      BStates.push_back(tmp);
    }
    assert( BStates.size() == BTags.size() );
  }
}
