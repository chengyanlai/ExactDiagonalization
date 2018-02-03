#ifndef __HDF5IO_HPP__
#define __HDF5IO_HPP__

#include <string>
#include <complex>
#include <vector>
#include <cassert>
#include <H5Cpp.h>
#include "src/ArmadilloMatrix.hpp"


class HDF5IO: public H5::H5File {
  private:
    std::string FileName;
    H5::CompType ComplexDataType;
    const H5::CompType InitComplexDataType(void);
    bool FileExists(const std::string &FileName);

    template<typename T>
    inline H5::DataType DetermineDataType(const T input){
      if ( sizeof(T) == 4){//int
        return H5::PredType::NATIVE_INT;
      }else if ( sizeof(T) == 8 ){//double
        return H5::PredType::NATIVE_DOUBLE;
      }else if ( sizeof(T) == 16 ){//complex<double>
        return ComplexDataType;
      }
    };

  public:
    HDF5IO (const std::string &FileName);
    HDF5IO (const std::string &FileName, const bool force);
    virtual ~HDF5IO (void){};
    inline H5::H5File GetFile(){return *this;};
    H5::Group GetGroup(const std::string &GroupName);

    template<typename T>
      void SaveRawBuffer(const std::string& GroupName, const std::string& SetName, const size_t dim, const T* x);
    template<typename T>
      void LoadRawBuffer(const std::string& GroupName, const std::string& SetName, size_t& dim, T*& x);

    template<typename T>
      inline void SaveNumber(const std::string& GroupName, const std::string& SetName, const T& x){
        size_t dim = 1;
        this->SaveRawBuffer(GroupName, SetName, dim, &x);
      };
    template<typename T>
      inline void LoadNumber(const std::string& GroupName, const std::string& SetName, T& x){
        T* val;
        size_t dim = 1;
        this->LoadRawBuffer(GroupName, SetName, dim, val);
        x = val[0];
      };

    template<typename T>
      inline void SaveStdVector(const std::string& GroupName, const std::string& SetName, const std::vector<T>& V){
        size_t dim = V.size();
        this->SaveRawBuffer(GroupName, SetName, dim, V.data());
      };
    template<typename T>
      inline void LoadStdVector(const std::string& GroupName, const std::string& SetName, std::vector<T>& V){
        T* val;
        size_t dim;
        this->LoadRawBuffer(GroupName, SetName, dim, val);
        V.clear();
        V.assign(val, val+dim);
      };

    template<typename T>
      void SaveVector(const std::string& GroupName, const arma::Col<T> Ten);
    template<typename T>
      void LoadVector(const std::string& GroupName, arma::Col<T>& Ten);
    template<typename T>
      void SaveMatrix(const std::string& GroupName, const arma::Mat<T> Ten);
    template<typename T>
      void LoadMaxtrix(const std::string& GroupName, arma::Mat<T>& Ten);
};

#endif  /* end of include guard: __HDF5IO_HPP__ */
