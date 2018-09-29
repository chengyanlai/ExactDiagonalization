#include <exception>
#include <sstream>
#include <cstdlib>
#include "src/hdf5io/hdf5io.hpp"

#define DEBUG 0

/* For the HDF5IO class */
const H5::CompType HDF5IO::InitComplexDataType(){
  H5::CompType Type(sizeof(std::complex<double>));
  // NOTE: New HDF5 with native complex datatypes ??
  Type.insertMember("real", 0, H5::PredType::NATIVE_DOUBLE);
  Type.insertMember("imag", sizeof(double), H5::PredType::NATIVE_DOUBLE);
  try {
    Type.commit(*this, "complex");
  } catch(H5::DataTypeIException){};
  return Type;
}

bool HDF5IO::FileExists(const std::string& FileName){
  try {
    H5::Exception::dontPrint();
    return isHdf5(FileName.c_str());
  } catch(H5::FileIException){
    return false;
  }
}

HDF5IO::HDF5IO (const std::string& fn) :
H5File(fn.c_str(), FileExists(fn) ? H5F_ACC_RDWR : H5F_ACC_TRUNC),
FileName(fn),
ComplexDataType(InitComplexDataType()){}

HDF5IO::HDF5IO (const std::string& fn, const bool force):
H5File(fn.c_str(), H5F_ACC_TRUNC),
FileName(fn),
ComplexDataType(InitComplexDataType()){}

H5::Group HDF5IO::GetGroup(const std::string &GroupName){
  H5::Group group;
  try{
    H5::Exception::dontPrint();
    group = H5::Group( this->openGroup( GroupName.c_str() ) );
  } catch( H5::FileIException not_found_error ){
    group = H5::Group( this->createGroup( GroupName.c_str() ) );
  } catch( const H5::Exception err ){
    std::cout << "In Group - " << GroupName << std::endl;
    throw std::runtime_error(" HDF5IO::GetGroup fail. ");
  }
  return group;
}


/* +------------------------------+
   | Load/Save C/C++ raw buffer   | Everything works from here!
   +------------------------------+  */
template<typename T>
void HDF5IO::SaveRawBuffer(const std::string& GroupName, const std::string& SetName, const size_t dim, const T* x){
  H5::DataType InputDataType = DetermineDataType( T(0) );
  try{
    H5::Exception::dontPrint();
    H5::DataSpace dataspace;// Default to H5S_SCALAR
    if ( dim != 1 ){
      hsize_t Dim[1] = {hsize_t(dim)};
      dataspace = H5::DataSpace(1,Dim);// Will change to H5S_SIMPLE
    }
    H5::Group FG = GetGroup( GroupName );
    try{
      H5::DataSet dset = FG.openDataSet(SetName.c_str());
      dset.write(x, InputDataType, dataspace);
    } catch ( const H5::GroupIException not_found_error ){
      H5::DataSet dset = FG.createDataSet(SetName.c_str(), InputDataType, dataspace);
      dset.write(x, InputDataType);
    } catch( const H5::FileIException error){
      error.printErrorStack();
    } catch( const H5::DataSetIException error){
      error.printErrorStack();
    }
    FG.close();
  }catch ( const H5::Exception ) {
    std::cout << "In Group - " << GroupName << ", and SetName is " << SetName << std::endl;
    throw std::runtime_error(" HDF5IO::SaveRawBuffer<T> ");
  }
}
template void HDF5IO::SaveRawBuffer(const std::string& GroupName, const std::string& SetName, const size_t dim, const int* x);
template void HDF5IO::SaveRawBuffer(const std::string& GroupName, const std::string& SetName, const size_t dim, const double* x);
template void HDF5IO::SaveRawBuffer(const std::string& GroupName, const std::string& SetName, const size_t dim, const std::complex<double>* x);
template<>
void HDF5IO::SaveRawBuffer(const std::string& GroupName, const std::string& SetName, const size_t dim, const unsigned long* x){
  H5::DataType InputDataType = H5::PredType::NATIVE_ULLONG;
  try{
    H5::Exception::dontPrint();
    H5::DataSpace dataspace;// Default to H5S_SCALAR
    if ( dim != 1 ){
      hsize_t Dim[1] = {hsize_t(dim)};
      dataspace = H5::DataSpace(1,Dim);// Will change to H5S_SIMPLE
    }
    H5::Group FG = GetGroup( GroupName );
    try{
      H5::DataSet dset = FG.openDataSet(SetName.c_str());
      dset.write(x, InputDataType, dataspace);
    } catch ( const H5::GroupIException not_found_error ){
      H5::DataSet dset = FG.createDataSet(SetName.c_str(), InputDataType, dataspace);
      dset.write(x, InputDataType);
    } catch( const H5::FileIException error){
      error.printErrorStack();
    } catch( const H5::DataSetIException error){
      error.printErrorStack();
    }
    FG.close();
  }catch ( const H5::Exception ) {
    std::cout << "In Group - " << GroupName << ", and SetName is " << SetName << std::endl;
    throw std::runtime_error(" HDF5IO::SaveRawBuffer<ulong> ");
  }
}

template<typename T>
void HDF5IO::LoadRawBuffer(const std::string& GroupName, const std::string& SetName, size_t& dim, T*& x){
  H5::DataType OutputDataType = DetermineDataType( T(0) );
  try{
    H5::Group FG = GetGroup( GroupName );
    H5::DataSet DataSet = FG.openDataSet(SetName.c_str());
    try{
      H5::DataSpace DataSpace = DataSet.getSpace();
      hsize_t Dims[1];
      int NDims = DataSpace.getSimpleExtentDims(Dims);
      dim = ( NDims == 0 )? 1 : Dims[0];
      x = (T*)malloc( dim * sizeof(T) );
      DataSet.read(x, OutputDataType);
    } catch( const H5::DataSpaceIException error ){
      error.printErrorStack();
    } catch( const H5::DataSetIException error){
      error.printErrorStack();
    } catch( const H5::FileIException error){
      error.printErrorStack();
    }
    FG.close();
  }catch ( const H5::Exception error ) {
    error.printErrorStack();
    std::cout << "In Group - " << GroupName << ", and SetName is " << SetName << std::endl;
    throw std::runtime_error(" HDF5IO::LoadRawBuffer<T>. ");
  }
}
template void HDF5IO::LoadRawBuffer(const std::string& GroupName, const std::string& SetName, size_t& dim, int*& x);
template void HDF5IO::LoadRawBuffer(const std::string& GroupName, const std::string& SetName, size_t& dim, double*& x);
template void HDF5IO::LoadRawBuffer(const std::string& GroupName, const std::string& SetName, size_t& dim, std::complex<double>*& x);
template<>
void HDF5IO::LoadRawBuffer(const std::string& GroupName, const std::string& SetName, size_t& dim, unsigned long*& x){
  H5::DataType OutputDataType = H5::PredType::NATIVE_ULONG;
  try{
    H5::Group FG = GetGroup( GroupName );
    H5::DataSet DataSet = FG.openDataSet(SetName.c_str());
    try{
      H5::DataSpace DataSpace = DataSet.getSpace();
      hsize_t Dims[1];
      int NDims = DataSpace.getSimpleExtentDims(Dims);
      dim = ( NDims == 0 )? 1 : Dims[0];
      x = (unsigned long*)malloc( dim * sizeof(unsigned long) );
      DataSet.read(x, OutputDataType);
    } catch( const H5::DataSpaceIException error ){
        error.printErrorStack();
    } catch( const H5::DataSetIException error){
      error.printErrorStack();
    } catch( const H5::FileIException error){
      error.printErrorStack();
    }
    FG.close();
  }catch ( const H5::Exception ) {
    std::cout << "In Group - " << GroupName << ", and SetName is " << SetName << std::endl;
    throw std::runtime_error(" HDF5IO::LoadRawBuffer<ulong>. ");
  }
}






/* +-----------------------------------+
   | Load/Save armadillo matrix/vector |
   +-----------------------------------+ */
template<typename T>
void HDF5IO::SaveVector(const std::string& GroupName, const std::string& SetName, const arma::Col<T> Vec){
  std::string gname = GroupName;
  H5::Group FG1 = GetGroup( GroupName );
  gname.append("/");
  gname.append(SetName);
  H5::Group FG2 = GetGroup( gname );
  int rows = Vec.n_rows;
  SaveNumber(gname, "col", rows);
  SaveRawBuffer(gname, "Elem", rows, Vec.memptr());
}
template void HDF5IO::SaveVector(const std::string& GroupName, const std::string& SetName, const arma::Col<double> Vec);
template void HDF5IO::SaveVector(const std::string& GroupName, const std::string& SetName, const arma::Col<std::complex<double> > Vec);

template<typename T>
void HDF5IO::LoadVector(const std::string& GroupName, const std::string& SetName, arma::Col<T>& Vec){
  std::string gname = GroupName;
  gname.append("/");
  gname.append(SetName);
  T* Elem;
  size_t ElemNum;
  LoadNumber(gname, "col", ElemNum);
  LoadRawBuffer(gname, "Elem", ElemNum, Elem);
  Vec = arma::Col<T>(Elem, ElemNum);//, copy_aux_mem = true, strict = false);
}
template void HDF5IO::LoadVector(const std::string& GroupName, const std::string& SetName, arma::Col<double>& Ten);
template void HDF5IO::LoadVector(const std::string& GroupName, const std::string& SetName, arma::Col<std::complex<double> >& Ten);

template<typename T>
void HDF5IO::SaveMatrix(const std::string& GroupName, const std::string& SetName, const arma::Mat<T> Mat){
  std::string gname = GroupName;
  H5::Group FG1 = GetGroup( GroupName );
  gname.append("/");
  gname.append(SetName);
  H5::Group FG2 = GetGroup( gname );
  int rows = Mat.n_rows;
  SaveNumber(gname, "row", rows);
  int cols = Mat.n_cols;
  SaveNumber(gname, "col", cols);
  SaveRawBuffer(gname, "Elem", Mat.n_rows*Mat.n_cols, Mat.memptr());
}
template void HDF5IO::SaveMatrix(const std::string& GroupName, const std::string& SetName, const arma::Mat<double> Mat);
template void HDF5IO::SaveMatrix(const std::string& GroupName, const std::string& SetName, const arma::Mat<std::complex<double> > Mat);

template<typename T>
void HDF5IO::LoadMatrix(const std::string& GroupName, const std::string& SetName, arma::Mat<T>& Mat){
  std::string gname = GroupName;
  gname.append("/");
  gname.append(SetName);
  T* Elem;
  size_t rows, cols;
  LoadNumber(gname, "col", cols);
  LoadNumber(gname, "row", rows);
  size_t NumElem = rows * cols;
  LoadRawBuffer(gname, "Elem", NumElem, Elem);
  Mat = arma::Mat<T>(Elem, rows, cols);//, copy_aux_mem = true, strict = false);
}
template void HDF5IO::LoadMatrix(const std::string& GroupName, const std::string& SetName, arma::Mat<double>& Mat);
template void HDF5IO::LoadMatrix(const std::string& GroupName, const std::string& SetName, arma::Mat<std::complex<double> >& Mat);
