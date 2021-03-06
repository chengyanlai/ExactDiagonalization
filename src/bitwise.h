#ifndef __BITWISE_H__
#define __BITWISE_H__

inline bool btest(const int i, const int j)
{
  /* This is equavalent to n_j = c^\dagger_jc_j */
  return ( i & (1 << j) );
}

inline int ibset(const int i, const int j)
{
  /* This is equivalent to c^\dagger_j */
  int number = i;
  number |= 1 << j;
  return number;
}

inline int ibclr(const int i, const int j)
{
  /* This is equivalent to c_j */
  int number = i;
  number &= ~(1 << j);
  return number;
}

#endif//__BITWISE_H__
