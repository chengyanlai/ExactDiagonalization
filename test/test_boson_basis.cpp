#include <iostream>
#include <vector>
#include "src/EDType.hpp"
#include "src/Basis/Basis.hpp"

int main(int argc, char const *argv[]) {
  int L = 4;
  int N = 4;
  Basis B1(L, N);
  B1.Boson();
  std::vector< std::vector<int> > st = B1.getBStates();
  std::vector< double > tg = B1.getBTags();
  for (size_t cnt = 0; cnt < tg.size(); cnt++) {
    INFO_NONEWLINE(cnt << " ");
    for (auto &j : st.at(cnt)){
      INFO_NONEWLINE(j << " ");
    }
    INFO(tg.at(cnt));
  }
  return 0;
}
