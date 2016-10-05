#ifndef __HDF5IO_HPP__
#define __HDF5IO_HPP__
#include <string>
#include <vector>
#include <H5Cpp.h>
#include "src/EDType.hpp"
#include "src/EigenMatrix.hpp"

class HDF5IO: public H5::H5File {
private:
  std::string FileName;
  H5::CompType ComplexDataType;
  const H5::CompType initCompexDataType(void);
  H5::EnumType PtrLREnumType;
  const H5::EnumType initPtrLREnumType(void);
  bool fileExists(const std::string &FileName);

public:
  HDF5IO (const std::string &FileName);
  HDF5IO (const std::string &FileName, const bool force);
  virtual ~HDF5IO ();
  inline H5::H5File getFile(){return *this;};
  H5::Group getGroup(const std::string &GroupName);

  void saveNumber(const std::string& GroupName, const std::string& Name, int x);
  void saveNumber(const std::string& GroupName, const std::string& Name, unsigned long x);
  void saveNumber(const std::string& GroupName, const std::string& Name, RealType x);
  void saveNumber(const std::string& GroupName, const std::string& Name, ComplexType C);

  void saveStdVector(const std::string& GroupName, const std::string& Name, const std::vector<int>& V);
  void saveStdVector(const std::string& GroupName, const std::string& Name, const std::vector<size_t>& V);
  void saveStdVector(const std::string& GroupName, const std::string& Name, const std::vector<double>& V);
  void saveStdVector(const std::string& GroupName, const std::string& Name, const std::vector<std::complex<double> >& V);
  // only for Eigen3 matrix
  void saveVector(const std::string& GroupName, const std::string& Name, const RealVectorType& V);
  void saveVector(const std::string& GroupName, const std::string& Name, const ComplexVectorType& V);
  void saveMatrix(const std::string& GroupName, const std::string& Name, const RealMatrixType& M);
  void saveMatrix(const std::string& GroupName, const std::string& Name, const ComplexMatrixType& M);

  int loadInt(const std::string& GroupName, const std::string& Name);
  size_t loadUlong(const std::string& GroupName, const std::string& Name);
  RealType loadReal(const std::string& GroupName, const std::string& Name);
  ComplexType loadComplex(const std::string& GroupName, const std::string& Name);
  void loadStdVector(const std::string& GroupName, const std::string& Name, std::vector<int>& V);
  void loadStdVector(const std::string& GroupName, const std::string& Name, std::vector<size_t>& V);
  void loadStdVector(const std::string& GroupName, const std::string& Name, std::vector<double>& V);
  void loadStdVector(const std::string& GroupName, const std::string& Name, std::vector<std::complex<double> >& V);
  // only for Eigen3 matrix
  void loadVector(const std::string& GroupName, const std::string& Name, RealVectorType& V);
  void loadVector(const std::string& GroupName, const std::string& Name, ComplexVectorType& V);
  void loadMatrix(const std::string& GroupName, const std::string& Name, RealMatrixType& M);
  void loadMatrix(const std::string& GroupName, const std::string& Name, ComplexMatrixType& M);
};
#endif//__HDF5IO_HPP__
