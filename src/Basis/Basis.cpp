#include "src/Basis/Basis.hpp"

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
