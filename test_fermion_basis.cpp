#include <iostream>
#include <vector>
#include "src/EDType.hpp"
#include "src/Basis/Basis.hpp"

int main(int argc, char const *argv[]) {
  int L = 4;
  int N = 2;
  Basis F1(L, N, true);
  F1.Fermion();
  std::vector<int> st = F1.getFStates();
  std::vector<size_t> tg = F1.getFTags();
  for (size_t cnt = 0; cnt < st.size(); cnt++) {
    INFO(cnt << " " << st.at(cnt) << " " << tg.at(st.at(cnt)));
  }
  return 0;
}
