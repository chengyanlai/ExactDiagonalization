ifeq "$(ARCH)" ""
ARCH = x86_64
endif

ifeq "$(OS)" ""
OS = $(shell uname -s)
endif

ifeq "$(NODE)" ""
NODE = $(shell uname -n)
endif

ifeq ("$(NODE)", "CHEN-YENs-MacBook-Pro.local")
	EIGENINC = /Volumes/Files/GitHub/eigen/
	HDF5ROOT = /usr/local/Cellar/hdf5/1.8.15/
	MKLROOT =
else ifeq ("$(NODE)", "kagome.rcc.ucmerced.edu")
	EIGENINC = /usr/local/include
	HDF5ROOT = /data2/hdf5/HDF_Group/HDF5/1.8.15-patch1
	MKLROOT = /data2/intel/composer_xe_2015.2.164/mkl
else ifeq ("$(NODE)", "edgestate.rcc.ucmerced.edu")
	EIGENINC = /usr/local/include
	HDF5ROOT = /opt/HDF_Group/HDF5/1.8.15-patch1
	MKLROOT = /opt/intel/composer_xe_2015.2.164/mkl
else ifeq ("$(NODE)", "atomtronics.ucmerced.edu")
	EIGENINC = /usr/local/include
	HDF5ROOT = /opt/HDF_Group/HDF5/1.8.15-patch1
	MKLROOT = /opt/intel/composer_xe_2015.2.164/mkl
else ifneq (, $(filter "$(NODE)", "comet-ln1.sdsc.edu" "comet-ln2.sdsc.edu" "comet-ln3.sdsc.edu" "comet-ln4.sdsc.edu"))
	EIGENINC = /opt/eigen/include
	HDF5ROOT = /home/chenyen/hdf5/HDF_Group/HDF5/1.8.15-patch1
	MKLROOT = /opt/intel/composer_xe_2015.2.164/mkl
else ifneq (, $(filter "$(NODE)", "login1.stampede.tacc.utexas.edu" "login2.stampede.tacc.utexas.edu" "login3.stampede.tacc.utexas.edu" "login4.stampede.tacc.utexas.edu"))
	EIGENINC = /home1/03731/chenyen/eigen
	HDF5ROOT = /home1/03731/chenyen/HDF_Group/HDF5/1.8.15-patch1
	MKLROOT = /opt/apps/intel/15/composer_xe_2015.2.164/mkl
endif


ifeq ("$(OS)", "Darwin")
	OPENMP =#-fopenmp
	HDF5LIB = -L$(HDF5ROOT)/lib -lhdf5 -lhdf5_cpp
	LAPACK = -lblas -llapack -lm
	LAPACK_OMP = $(LAPACK)
	CC = clang++ $(OPENMP) -O3 -m64 -std=c++11 -stdlib=libc++ -I$(HDF5ROOT)/include -I$(EIGENINC)
else ifeq ("$(OS)", "Linux")
	OPENMP =#-openmp
	#NOTE: the order of linker matters!
	HDF5LIB = $(HDF5ROOT)/lib/libhdf5_cpp.a $(HDF5ROOT)/lib/libhdf5.a $(HDF5ROOT)/lib/libz.a
	LAPACK = $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a \
	$(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group \
	$(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
	$(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a \
	-Wl,--end-group -lpthread -lm
	LAPACK_OMP = $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a \
	$(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group \
	${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
	${MKLROOT}/lib/intel64/libmkl_core.a \
	${MKLROOT}/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -lpthread -lm
	CC = icpc $(OPENMP) -O3 -Wall -std=c++11 -I./ -I$(HDF5ROOT)/include -I$(EIGENINC) -DMKL
endif

MODULES   := Node Lattice Basis Hamiltonian Lanczos hdf5io
SRC_DIR   := $(addprefix src/,$(MODULES))
BUILD_DIR := $(addprefix build/,$(MODULES))

SRC       := $(foreach sdir,$(SRC_DIR),$(wildcard $(sdir)/*.cpp))
OBJ       := $(patsubst src/%.cpp,build/%.o,$(SRC))
# INCLUDES  := $(addprefix -I,$(SRC_DIR))
INCLUDES  := -I./

vpath %.cpp $(SRC_DIR)

define make-goal
$1/%.o: %.cpp
	$(CC) $(INCLUDES) -c $$< -o $$@
endef

.PHONY: all checkdirs clean

all: checkdirs build/SSH.f build/SSH.b build/test.fermion build/test.boson build/test.lattice build/test.node

build/%.o: %.cpp
	$(CC) $(INCLUDES) -c $< -o $@

build/SSH.b: build/SSH_boson.o $(OBJ)
	$(CC) $^ -o $@ $(LAPACK) $(HDF5LIB)

build/SSH.f: build/SSH_fermion.o $(OBJ)
	$(CC) $^ -o $@ $(LAPACK) $(HDF5LIB)

# build/test.fermion: build/test_fermion_hamiltonian.o $(OBJ)
build/test.fermion: build/test_fermion_complex.o $(OBJ)
	$(CC) $^ -o $@ $(LAPACK) $(HDF5LIB)

build/test.boson: build/test_boson_hamiltonian.o $(OBJ)
	$(CC) $^ -o $@ $(LAPACK) $(HDF5LIB)

build/test.lattice: build/test_lattice.o $(OBJ)
	$(CC) $^ -o $@ $(LAPACK) $(HDF5LIB)

build/test.node: build/test_node.o $(OBJ)
	$(CC) $^ -o $@ $(LAPACK) $(HDF5LIB)

checkdirs: $(BUILD_DIR)

$(BUILD_DIR):
	@mkdir -p $@

clean:
	@rm -rf $(BUILD_DIR) build/*.o

$(foreach bdir,$(BUILD_DIR),$(eval $(call make-goal,$(bdir))))
