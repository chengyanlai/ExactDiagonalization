#ifndef __EDTYPE_H__
#define __EDTYPE_H__
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <complex>
#include <cmath>
#include <vector>

#define INFO(MSG)             std::cout << MSG << std::endl
#define INFO_NONEWLINE(MSG)   std::cout << MSG << std::flush

#define _S(x) #x
#define S_(x) _S(x)
#define S__LINE__ S_(__LINE__)
/* use S__LINE__ instead of __LINE__ */
#define RUNTIME_ERROR(MSG)    throw std::runtime_error(MSG __FILE__ ": " S__LINE__ ": ")
#define LOGIC_ERROR(MSG)      throw std::logic_error(MSG __FILE__ ": " S__LINE__ ": ");
#define OVERFLOW_ERROR(MSG)   throw std::overflow_error(MSG __FILE__ ": " S__LINE__ ": ");

/** Real floating point type. */
typedef double RealType;
/** Complex type. */
typedef std::complex<double> ComplexType;

inline RealType Conjg( RealType number ){return number;};
inline ComplexType Conjg( ComplexType number ){return std::conj(number);};
inline RealType RealPart( RealType number ){return number;};
inline RealType RealPart( ComplexType number ){return number.real();};
inline RealType ImaginaryPart( RealType number ){return 0.0e0;};
inline RealType ImaginaryPart( ComplexType number ){return number.imag();};

const RealType PI = std::atan(1.0)*4.0;

template <typename Tnum> std::vector<Tnum> vec(const Tnum& a0) ;
template <typename Tnum> std::vector<Tnum> vec(const Tnum& a0, const Tnum& a1) ;
template <typename Tnum> std::vector<Tnum> vec(const Tnum& a0, const Tnum& a1, const Tnum& a2) ;
template <typename Tnum> std::vector<Tnum> vec(const Tnum& a0, const Tnum& a1, const Tnum& a2, const Tnum& a3) ;
template <typename Tnum> std::vector<Tnum> vec(const Tnum& a0, const Tnum& a1, const Tnum& a2, const Tnum& a3, const Tnum& a4) ;
template <typename Tnum> std::vector<Tnum> vec(const Tnum& a0, const Tnum& a1, const Tnum& a2, const Tnum& a3, const Tnum& a4, const Tnum& a5) ;
template <typename Tnum> std::vector<Tnum> vec(const Tnum& a0, const Tnum& a1, const Tnum& a2, const Tnum& a3, const Tnum& a4, const Tnum& a5, const Tnum& a6) ;
template <typename Tnum> std::vector<Tnum> vec(const Tnum& a0, const Tnum& a1, const Tnum& a2, const Tnum& a3, const Tnum& a4, const Tnum& a5, const Tnum& a6, const Tnum& a7) ;
template <typename Tnum> std::vector<Tnum> vec(const Tnum& a0, const Tnum& a1, const Tnum& a2, const Tnum& a3, const Tnum& a4, const Tnum& a5, const Tnum& a6, const Tnum& a7, const Tnum& a8) ;

template <typename Tnum>
std::vector<Tnum> vec(const Tnum& a0){
  std::vector<Tnum> A ;
  A.push_back(a0) ;
  return A ;
}

template <typename Tnum>
std::vector<Tnum> vec(const Tnum& a0, const Tnum& a1){
  std::vector<Tnum> A ;
  A.push_back(a0) ;
  A.push_back(a1) ;
  return A ;
}

template <typename Tnum>
std::vector<Tnum> vec(const Tnum& a0, const Tnum& a1, const Tnum& a2){
  std::vector<Tnum> A ;
  A.push_back(a0) ;
  A.push_back(a1) ;
  A.push_back(a2) ;
  return A ;
}

template <typename Tnum>
std::vector<Tnum> vec(const Tnum& a0, const Tnum& a1, const Tnum& a2,
  const Tnum& a3){
  std::vector<Tnum> A ;
  A.push_back(a0) ;
  A.push_back(a1) ;
  A.push_back(a2) ;
  A.push_back(a3) ;
  return A ;
}

template <typename Tnum>
std::vector<Tnum> vec(const Tnum& a0, const Tnum& a1, const Tnum& a2,
  const Tnum& a3, const Tnum& a4){
  std::vector<Tnum> A ;
  A.push_back(a0) ;
  A.push_back(a1) ;
  A.push_back(a2) ;
  A.push_back(a3) ;
  A.push_back(a4) ;
  return A ;
}

template <typename Tnum>
std::vector<Tnum> vec(const Tnum& a0, const Tnum& a1, const Tnum& a2,
  const Tnum& a3, const Tnum& a4, const Tnum& a5){
  std::vector<Tnum> A ;
  A.push_back(a0) ;
  A.push_back(a1) ;
  A.push_back(a2) ;
  A.push_back(a3) ;
  A.push_back(a4) ;
  A.push_back(a5) ;
  return A ;
}

template <typename Tnum>
std::vector<Tnum> vec(const Tnum& a0, const Tnum& a1, const Tnum& a2,
  const Tnum& a3, const Tnum& a4, const Tnum& a5, const Tnum& a6){
  std::vector<Tnum> A ;
  A.push_back(a0) ;
  A.push_back(a1) ;
  A.push_back(a2) ;
  A.push_back(a3) ;
  A.push_back(a4) ;
  A.push_back(a5) ;
  A.push_back(a6) ;
  return A ;
}

template <typename Tnum>
std::vector<Tnum> vec(const Tnum& a0, const Tnum& a1, const Tnum& a2,
  const Tnum& a3, const Tnum& a4, const Tnum& a5, const Tnum& a6, const Tnum& a7){
  std::vector<Tnum> A ;
  A.push_back(a0) ;
  A.push_back(a1) ;
  A.push_back(a2) ;
  A.push_back(a3) ;
  A.push_back(a4) ;
  A.push_back(a5) ;
  A.push_back(a6) ;
  A.push_back(a7) ;
  return A ;
}

template <typename Tnum>
std::vector<Tnum> vec(const Tnum& a0, const Tnum& a1, const Tnum& a2,
  const Tnum& a3, const Tnum& a4, const Tnum& a5, const Tnum& a6, const Tnum& a7,
  const Tnum& a8){
  std::vector<Tnum> A ;
  A.push_back(a0) ;
  A.push_back(a1) ;
  A.push_back(a2) ;
  A.push_back(a3) ;
  A.push_back(a4) ;
  A.push_back(a5) ;
  A.push_back(a6) ;
  A.push_back(a7) ;
  A.push_back(a8) ;
  return A ;
}

// template <typename T>
// void PrintVector(const std::vector<T> Vin, const int digit=9, const std::string separation="\t"){
//   for ( size_t i = 0; i < Vin.size(); i++){
//     std::cout << "\t" << std::scientific << std::setprecision(digit) << Vin.at(i) << separation << std::flush;
//   }
//   std::cout << std::endl;
// }

// template <typename T>
// void PrintVector(std::ofstream& file, const std::vector<T> Vin, const int digit=9, const std::string separation="\n"){
//   for ( size_t i = 0; i < Vin.size(); i++){
//     file << "\t" << std::setw(2) << i << " " << std::scientific << std::setprecision(digit) << Vin.at(i) << separation << std::flush;
//   }
//   file << "\n";
// }

#endif//__EDTYPE_H__
