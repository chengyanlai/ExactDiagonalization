#include <exception>
#include <stdexcept>
#include <sstream>
#include <cassert>
#include <cstring>
#include "src/hdf5io/hdf5io.hpp"

/* For the HDF5IO class */
const H5::CompType HDF5IO::initCompexDataType()
{
  H5::CompType Type(sizeof(ComplexType));
  // TODO: New HDF5 with native complex datatypes ??
  Type.insertMember("real",0,H5::PredType::NATIVE_DOUBLE);
  Type.insertMember("imag",sizeof(RealType),H5::PredType::NATIVE_DOUBLE);
  try {
      Type.commit(*this,"complex");
  } catch(H5::DataTypeIException){};
  return Type;
}

bool HDF5IO::fileExists(const std::string& FileName)
{
  try {
    H5::Exception::dontPrint();
    return isHdf5(FileName.c_str());
  } catch(H5::FileIException){
    return false;
  }
}

HDF5IO::HDF5IO (const std::string& fn) :
FileName(fn),
H5File(fn.c_str(), fileExists(fn) ? H5F_ACC_RDWR : H5F_ACC_TRUNC),
ComplexDataType(initCompexDataType())
{}

HDF5IO::HDF5IO (const std::string& fn, const bool force):
FileName(fn),
H5File(fn.c_str(), H5F_ACC_TRUNC),
ComplexDataType(initCompexDataType())
{}

HDF5IO::~HDF5IO(void)
{
  this->close();
}

H5::Group HDF5IO::getGroup(const std::string &GroupName)
{
  H5::Group group;
  try{
    group = H5::Group( this->openGroup( GroupName.c_str() ) );
    // INFO("Open old group");
  } catch( H5::FileIException not_found_error ){
    // INFO("Create new group");
    group = H5::Group( this->createGroup( GroupName.c_str() ) );
  } catch( const H5::Exception err ){
    RUNTIME_ERROR("HDF5IO::getGroup");
  }
  return group;
}

void HDF5IO::saveNumber(const std::string& GroupName, const std::string& Name,
    int x)
{
    H5::Group FG = getGroup( GroupName );
    try{
      H5::Exception::dontPrint();
      H5::DataSet dataset = FG.openDataSet( Name.c_str() );
      dataset.write(&x, H5::PredType::NATIVE_INT);
    } catch ( const H5::GroupIException not_found_error ){
      H5::DataSet dataset = FG.createDataSet( Name.c_str(), H5::PredType::NATIVE_INT, H5::DataSpace());
      dataset.write(&x, H5::PredType::NATIVE_INT);
    }
    FG.close();
}

int HDF5IO::loadInt(const std::string& GroupName, const std::string& Name)
{
    try{
      H5::Group FG = getGroup( GroupName );
      H5::DataSet DataSet = FG.openDataSet( Name.c_str());
      int x;
      DataSet.read(&x,H5::PredType::NATIVE_INT);
      FG.close();
      return x;
    }catch( H5::GroupIException not_found_error ){
      RUNTIME_ERROR("No dataset found in loadInt. ");
    }
}

void HDF5IO::saveNumber(const std::string& GroupName, const std::string& Name,
    unsigned long x)
{
    H5::Group FG = getGroup( GroupName );
    try{
      H5::Exception::dontPrint();
      H5::DataSet dataset = FG.openDataSet( Name.c_str() );
      dataset.write(&x, H5::PredType::NATIVE_ULONG);
    } catch ( const H5::GroupIException not_found_error ){
      H5::DataSet dataset = FG.createDataSet( Name.c_str(), H5::PredType::NATIVE_ULONG, H5::DataSpace());
      dataset.write(&x, H5::PredType::NATIVE_ULONG);
    }
    FG.close();
}

size_t HDF5IO::loadUlong(const std::string& GroupName, const std::string& Name)
{
  try{
    H5::Group FG = getGroup( GroupName );
    H5::DataSet DataSet = FG.openDataSet( Name.c_str() );
    size_t x;
    DataSet.read(&x, H5::PredType::NATIVE_ULONG);
    return x;
    FG.close();
  }catch( H5::GroupIException not_found_error ){
    INFO("In Group - " << GroupName << ", and Name is " << Name);
    RUNTIME_ERROR("No dataset found in loadUlong. ");
  }
}

// Save/load RealType
void HDF5IO::saveNumber(const std::string& GroupName, const std::string& Name,
    RealType x)
{
    H5::Group FG = getGroup( GroupName );
    try{
      H5::Exception::dontPrint();
      H5::DataSet dataset = FG.openDataSet(Name.c_str());
      dataset.write(&x,H5::PredType::NATIVE_DOUBLE);
    } catch ( const H5::GroupIException not_found_error ) {
      H5::DataSet dataset = FG.createDataSet(Name.c_str(),H5::PredType::NATIVE_DOUBLE,H5::DataSpace());
      dataset.write(&x,H5::PredType::NATIVE_DOUBLE);
    }
    FG.close();
}

RealType HDF5IO::loadReal(const std::string& GroupName, const std::string& Name)
{
  try{
    H5::Group FG = getGroup( GroupName );
    H5::DataSet DataSet = FG.openDataSet(Name.c_str());
    RealType x;
    DataSet.read(&x,H5::PredType::NATIVE_DOUBLE);
    FG.close();
    return x;
  }catch( H5::GroupIException not_found_error ){
    RUNTIME_ERROR("No dataset found in loadReal. ");
  }
}

void HDF5IO::saveNumber(const std::string& GroupName, const std::string& Name,
    ComplexType C)
{
    H5::CompType ComplexDataType = openCompType("complex");
    H5::Group FG = getGroup( GroupName );
    try{
      H5::Exception::dontPrint();
      H5::DataSet dataset = FG.openDataSet(Name.c_str());
      RealType RealImag[2] = {real(C),imag(C)};
      dataset.write(RealImag, ComplexDataType);
    } catch ( const H5::GroupIException not_found_error ){
      H5::DataSet dataset = FG.createDataSet(Name.c_str(), ComplexDataType, H5::DataSpace());
      RealType RealImag[2] = {real(C),imag(C)};
      dataset.write(RealImag, ComplexDataType);
    }
    FG.close();
}

ComplexType HDF5IO::loadComplex(const std::string& GroupName, const std::string& Name)
{
  try{
    H5::CompType ComplexDataType = this->openCompType("complex");
    H5::Group FG = getGroup( GroupName );
    H5::DataSet DataSet = FG.openDataSet(Name.c_str());
    ComplexType C;
    RealType RealImag[2];
    DataSet.read(RealImag, ComplexDataType);
    FG.close();
    return ComplexType(RealImag[0],RealImag[1]);
  }catch( H5::GroupIException not_found_error ){
    RUNTIME_ERROR("No dataset found in loadComplex. ");
  }
}

void HDF5IO::saveStdVector(const std::string& GroupName, const std::string& Name,
    const std::vector<int>& V)
{
  try{
    hsize_t Dim[1] = {hsize_t(V.size())};
    H5::DataSpace dspace(1,Dim);
    H5::Group FG = getGroup( GroupName );
    try{
      H5::Exception::dontPrint();
      H5::DataSet DataSet = FG.openDataSet(Name.c_str());
      DataSet.write(V.data(),H5::PredType::NATIVE_INT, dspace);
    } catch ( const H5::GroupIException not_found_error ){
      H5::DataSet DataSet = FG.createDataSet(Name.c_str(),
        H5::PredType::NATIVE_INT, dspace);
      DataSet.write(V.data(),H5::PredType::NATIVE_INT);
    }
    FG.close();
  } catch( const H5::Exception err ){
    RUNTIME_ERROR("HDF5IO::saveUlongStdVector");
  }
}

void HDF5IO::saveStdVector(const std::string& GroupName, const std::string& Name,
    const std::vector<size_t>& V)
{
  try{
    hsize_t Dim[1] = {hsize_t(V.size())};
    H5::DataSpace dspace(1,Dim);
    H5::Group FG = getGroup( GroupName );
    try{
      H5::Exception::dontPrint();
      H5::DataSet DataSet = FG.openDataSet(Name.c_str());
      DataSet.write(V.data(),H5::PredType::NATIVE_ULONG, dspace);
    } catch ( const H5::GroupIException not_found_error ){
      H5::DataSet DataSet = FG.createDataSet(Name.c_str(),
        H5::PredType::NATIVE_ULONG, dspace);
      DataSet.write(V.data(), H5::PredType::NATIVE_ULONG);
    }
    FG.close();
  } catch( const H5::Exception err ){
    RUNTIME_ERROR("HDF5IO::saveUlongStdVector ");
  }
}

void HDF5IO::saveStdVector(const std::string& GroupName, const std::string& Name,
    const std::vector<double>& V)
{
  try{
    hsize_t Dim[1] = {hsize_t(V.size())};
    H5::DataSpace dataspace(1,Dim);
    H5::Group FG = getGroup( GroupName );
    try{
      H5::Exception::dontPrint();
      H5::DataSet dataset = FG.openDataSet(Name.c_str());
      dataset.write(V.data(),H5::PredType::NATIVE_DOUBLE, dataspace);
    } catch ( const H5::GroupIException not_found_error ){
      H5::DataSet dataset = FG.createDataSet(Name.c_str(),
        H5::PredType::NATIVE_DOUBLE,dataspace);
      dataset.write(V.data(),H5::PredType::NATIVE_DOUBLE);
    }
    FG.close();
  } catch( const H5::Exception err ){
    RUNTIME_ERROR("HDF5IO::saveRealStdVector");
  }
}

void HDF5IO::saveStdVector(const std::string& GroupName, const std::string& Name,
    const std::vector<std::complex<double> >& V)
{
  try{
    H5::CompType ComplexDataType = openCompType("complex");
    hsize_t Dim[1] = {hsize_t(V.size())};
    H5::DataSpace dataspace(1,Dim);
    H5::Group FG = getGroup( GroupName.c_str() );
    try{
      H5::Exception::dontPrint();
      H5::DataSet dataset = FG.openDataSet(Name.c_str());
      dataset.write(V.data(), ComplexDataType, dataspace);
    } catch( const H5::GroupIException not_found_error ){
      H5::DataSet dataset = FG.createDataSet(Name.c_str(), ComplexDataType,
        dataspace);
      dataset.write(V.data(), ComplexDataType);
    } catch( const H5::FileIException error){
      error.printError();
    } catch( const H5::DataSetIException error){
      error.printError();
    }
    FG.close();
  } catch( const H5::Exception err ){
    err.printError();
    RUNTIME_ERROR("HDF5IO::saveComplexStdVector. ");
  }
}

void HDF5IO::loadStdVector(const std::string& GroupName, const std::string& Name,
    std::vector<size_t>& V)
{
  try{
    H5::Group FG = getGroup( GroupName );
    H5::DataSet DataSet = FG.openDataSet(Name.c_str());
    H5::DataSpace DataSpace = DataSet.getSpace();
    if(DataSpace.getSimpleExtentNdims() != 1)
      throw(H5::DataSpaceIException("HDF5IO::loadRealVector()",
        "Unexpected multidimentional dataspace."));
    V.resize(DataSpace.getSimpleExtentNpoints());
    DataSet.read(V.data(),H5::PredType::NATIVE_ULONG);
    FG.close();
  } catch( const H5::Exception err ){
    RUNTIME_ERROR("HDF5IO::loadUlongStdVector");
  }
}

void HDF5IO::loadStdVector(const std::string& GroupName, const std::string& Name,
    std::vector<int>& V)
{
  try{
    H5::Group FG = getGroup( GroupName );
    H5::DataSet DataSet = FG.openDataSet(Name.c_str());
    H5::DataSpace DataSpace = DataSet.getSpace();
    if(DataSpace.getSimpleExtentNdims() != 1)
      throw(H5::DataSpaceIException("HDF5IO::loadRealVector()",
        "Unexpected multidimentional dataspace."));
    V.resize(DataSpace.getSimpleExtentNpoints());
    DataSet.read(V.data(),H5::PredType::NATIVE_INT);
    FG.close();
  } catch( const H5::Exception err ){
    RUNTIME_ERROR("HDF5IO::loadUlongStdVector");
  }
}

void HDF5IO::loadStdVector(const std::string& GroupName, const std::string& Name,
    std::vector<RealType>& V)
{
  try{
    H5::Group FG = getGroup( GroupName );
    H5::DataSet DataSet = FG.openDataSet(Name.c_str());
    H5::DataSpace DataSpace = DataSet.getSpace();
    if(DataSpace.getSimpleExtentNdims() != 1)
      throw(H5::DataSpaceIException("HDF5IO::loadRealVector()","Unexpected multidimentional dataspace."));
    V.resize(DataSpace.getSimpleExtentNpoints());
    DataSet.read(V.data(),H5::PredType::NATIVE_DOUBLE);
    FG.close();
  } catch( const H5::Exception err ){
    RUNTIME_ERROR("HDF5IO::loadRealStdVector");
  }
}

void HDF5IO::loadStdVector(const std::string& GroupName, const std::string& Name,
    std::vector<ComplexType>& V)
{
  try{
    H5::CompType ComplexDataType = this->openCompType("complex");
    H5::Group FG = getGroup( GroupName );
    H5::DataSet DataSet = FG.openDataSet(Name.c_str());
    H5::DataSpace DataSpace = DataSet.getSpace();
    if(DataSpace.getSimpleExtentNdims() != 1)
      throw(H5::DataSpaceIException("HDF5IO::loadComplexVector()","Unexpected multidimentional dataspace."));
    V.resize(DataSpace.getSimpleExtentNpoints());
    DataSet.read(V.data(),ComplexDataType);
    FG.close();
  } catch( const H5::Exception err ){
    RUNTIME_ERROR("HDF5IO::loadComplexStdVector");
  }
}

/* only for Eigen3 matrix/vector */
void HDF5IO::saveVector(const std::string& GroupName,
  const std::string& Name, const RealVectorType& V)
{
    hsize_t Dim[1] = {hsize_t(V.size())};
    H5::DataSpace dspace(1,Dim);
    H5::Group FG = getGroup( GroupName );
    try{
      H5::Exception::dontPrint();
      H5::DataSet DataSet = FG.openDataSet(Name.c_str());
      DataSet.write(V.data(),H5::PredType::NATIVE_DOUBLE, dspace);
    } catch ( const H5::GroupIException not_found_error ){
      H5::DataSet DataSet = FG.createDataSet(Name.c_str(),
        H5::PredType::NATIVE_DOUBLE,dspace);
      DataSet.write(V.data(),H5::PredType::NATIVE_DOUBLE);
    }
    FG.close();
}

void HDF5IO::loadVector(const std::string& GroupName, const std::string& Name,
  RealVectorType& V)
{
  H5::Group FG = getGroup( GroupName );
  H5::DataSet DataSet = FG.openDataSet(Name.c_str());
  H5::DataSpace DataSpace = DataSet.getSpace();
  if(DataSpace.getSimpleExtentNdims() != 1)
	throw(H5::DataSpaceIException("HDF5IO::loadRealVector()",
    "Unexpected multidimentional dataspace."));
  V.resize(DataSpace.getSimpleExtentNpoints());
  try{
    DataSet.read(V.data(),H5::PredType::NATIVE_DOUBLE);
  }catch( H5::GroupIException not_found_error ){
    RUNTIME_ERROR("No dataset found in loadRealVector. ");
  }
  FG.close();
}

void HDF5IO::saveVector(const std::string& GroupName,
  const std::string& Name, const ComplexVectorType& V)
{
  try{
    H5::CompType ComplexDataType = this->openCompType("complex");
    hsize_t Dim[1] = {hsize_t(V.size())};
    H5::DataSpace dspace(1,Dim);
    H5::Group FG = getGroup( GroupName );
    try{
      H5::Exception::dontPrint();
      H5::DataSet DataSet = FG.openDataSet(Name.c_str());
      DataSet.write(V.data(), ComplexDataType, dspace);
    } catch ( const H5::GroupIException not_found_error ){
      H5::DataSet DataSet = FG.createDataSet(Name.c_str(), ComplexDataType,
        dspace);
      DataSet.write(V.data(), ComplexDataType);
    }
    FG.close();
  } catch ( const H5::DataSetIException error ){
    error.printError();
    RUNTIME_ERROR("HDF5IO::saveComplexVector at ");
  }
}

void HDF5IO::loadVector(const std::string& GroupName, const std::string& Name,
  ComplexVectorType& V)
{
  try{
    H5::CompType ComplexDataType = this->openCompType("complex");
    H5::Group FG = getGroup( GroupName );
    H5::DataSet DataSet = FG.openDataSet(Name.c_str());
    H5::DataSpace DataSpace = DataSet.getSpace();
    if(DataSpace.getSimpleExtentNdims() != 1)
      throw(H5::DataSpaceIException("HDF5IO::loadRealVector()",
        "Unexpected multidimentional dataspace."));
    V.resize(DataSpace.getSimpleExtentNpoints());
    try{
      DataSet.read(V.data(), ComplexDataType);
    }catch( H5::GroupIException not_found_error ){
      RUNTIME_ERROR("No dataset found in loadRealVector. ");
    }
    FG.close();
  } catch( const H5::Exception err ){
    RUNTIME_ERROR("HDF5IO::loadComplexVector at ");
  }
}

void HDF5IO::saveMatrix(const std::string& GroupName, const std::string& Name,
    const RealMatrixType& M)
{
  try{
    hsize_t Dims[2] = {hsize_t(M.rows()),hsize_t(M.cols())};
    H5::DataSpace dataspace(2,Dims);
    H5::Group FG = getGroup( GroupName );
    try{
      H5::Exception::dontPrint();
      H5::DataSet DataSet = FG.openDataSet(Name.c_str());
      DataSet.write(M.data(),H5::PredType::NATIVE_DOUBLE,dataspace);
    } catch ( const H5::GroupIException not_found_error ){
      H5::DataSet DataSet = FG.createDataSet(Name.c_str(),H5::PredType::NATIVE_DOUBLE,dataspace);
      DataSet.write(M.data(),H5::PredType::NATIVE_DOUBLE);
    }
    FG.close();
  } catch( const H5::Exception err ){
    RUNTIME_ERROR("HDF5IO::saveRealMatrix");
  }
}

void HDF5IO::loadMatrix(const std::string& GroupName, const std::string& Name,
    RealMatrixType& M)
{
  try{
    H5::Group FG = getGroup( GroupName );
    H5::DataSet DataSet = FG.openDataSet(Name.c_str());
    H5::DataSpace DataSpace = DataSet.getSpace();
    if( DataSpace.getSimpleExtentNdims() != 2 )
      throw(H5::DataSpaceIException("HDF5IO::loadRealMatrix()",
        "A dataspace must be precisely two-dimensional."));
    hsize_t Dims[2];
    DataSpace.getSimpleExtentDims(Dims);
    M.resize(Dims[0], Dims[1]);
    DataSet.read(M.data(), H5::PredType::NATIVE_DOUBLE);
    FG.close();
  } catch( const H5::Exception err ){
    RUNTIME_ERROR("HDF5IO::loadRealMatrix");
  }
}

void HDF5IO::saveMatrix(const std::string& GroupName, const std::string& Name,
    const ComplexMatrixType& M)
{
  try{
    H5::CompType ComplexDataType = this->openCompType("complex");
    hsize_t Dims[2] = {hsize_t(M.rows()),hsize_t(M.cols())};
    H5::DataSpace dataspace(2,Dims);
    H5::Group FG = getGroup( GroupName );
    try{
      H5::Exception::dontPrint();
      H5::DataSet dset = FG.openDataSet(Name.c_str());
      // dset.extend( Dims );not working
      dset.write(M.data(), ComplexDataType);
    } catch ( const H5::GroupIException not_found_error ){
      H5::DataSet dset = FG.createDataSet(Name.c_str(), ComplexDataType, dataspace);
      dset.write(M.data(), ComplexDataType);
    } catch ( const H5::DataSetIException error ){
      error.printError();
      RUNTIME_ERROR("HDF5IO::saveComplexMatrix at ");
    }
    FG.close();
  } catch( const H5::Exception error ){
    error.printError();
    RUNTIME_ERROR("HDF5IO::saveComplexMatrix at ");
  }
}

void HDF5IO::loadMatrix(const std::string& GroupName, const std::string& Name,
    ComplexMatrixType& M)
{
  try{
    H5::CompType ComplexDataType = this->openCompType("complex");
    H5::Group FG = getGroup( GroupName );
    H5::DataSet DataSet = FG.openDataSet(Name.c_str());
    H5::DataSpace DataSpace = DataSet.getSpace();
    if(DataSpace.getSimpleExtentNdims() != 2)
	throw(H5::DataSpaceIException("HDF5IO::loadMatrix()","A dataspace must be precisely two-dimensional."));
    hsize_t Dims[2];
    DataSpace.getSimpleExtentDims(Dims);
    M.resize(Dims[0],Dims[1]);
    DataSet.read(M.data(), ComplexDataType);
    FG.close();
  } catch( const H5::Exception err ){
    RUNTIME_ERROR("HDF5IO::loadComplexMatrix at ");
  }
}
