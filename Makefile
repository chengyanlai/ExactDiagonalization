ifeq "$(ARCH)" ""
ARCH = x86_64
endif

ifeq "$(OS)" ""
OS = $(shell uname -s)
endif

ifeq "$(NODE)" ""
NODE = $(shell uname -n)
endif

ifeq ("$(OS)", "Darwin")
	EIGENINC = /Volumes/Files/GitHub/eigen
	HDF5ROOT = /usr/local/opt/hdf5
	MKLROOT =
else ifeq ("$(NODE)", "kagome.rcc.ucmerced.edu")
	EIGENINC = /usr/local/include
	HDF5ROOT = /condensate1/hdf5/HDF_Group/HDF5/1.8.15-patch1
	MKLROOT = /condensate1/intel/composer_xe_2015.2.164/mkl
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
	LAPACK = -lblas -llapack -lm -larpack
	LAPACK_OMP = $(LAPACK)
	CC = clang++ $(OPENMP) -O3 -m64 -std=c++11 -stdlib=libc++ -I$(HDF5ROOT)/include -I$(EIGENINC)
else ifeq ("$(OS)", "Linux")
	OPENMP = -openmp
	#NOTE: the order of linker matters!
	HDF5LIB = $(HDF5ROOT)/lib/libhdf5_cpp.a $(HDF5ROOT)/lib/libhdf5.a $(HDF5ROOT)/lib/libz.a
	LAPACK = $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a \
	$(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group \
	$(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
	$(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a \
	-Wl,--end-group -lpthread -lm -larpack
	LAPACK_OMP = $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a \
	$(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group \
	${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
	${MKLROOT}/lib/intel64/libmkl_core.a \
	${MKLROOT}/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -lpthread -lm -larpack
	CC = icpc $(OPENMP) -O3 -Wall -std=c++11 -I./ -I$(HDF5ROOT)/include -I$(EIGENINC) -DMKL
endif

MODULES   := Node Lattice Basis Hamiltonian Lanczos hdf5io
SRC_DIR   := $(addprefix src/,$(MODULES))
SRC       := $(foreach sdir,$(SRC_DIR),$(wildcard $(sdir)/*.cpp))
OBJ       := $(patsubst src/%.cpp,build/%.o,$(SRC))

TB_MODULES   := TerminatorBeam
TB_SRC_DIR   := $(addprefix src/,$(MODULES) $(TB_MODULES))
TB_SRC       := $(foreach sdir,$(TB_SRC_DIR),$(wildcard $(sdir)/*.cpp))
TB_OBJ       := $(patsubst src/%.cpp,build/%.o,$(TB_SRC))

DP_MODULES   := Dephasing
DP_SRC_DIR   := $(addprefix src/,$(MODULES) $(DP_MODULES))
DP_SRC       := $(foreach sdir,$(DP_SRC_DIR),$(wildcard $(sdir)/*.cpp))
DP_OBJ       := $(patsubst src/%.cpp,build/%.o,$(DP_SRC))

OP_MODULES   := Lindblad
OP_SRC_DIR   := $(addprefix src/,$(MODULES) $(OP_MODULES))
OP_SRC       := $(foreach sdir,$(OP_SRC_DIR),$(wildcard $(sdir)/*.cpp))
OP_OBJ       := $(patsubst src/%.cpp,build/%.o,$(OP_SRC))

BUILD_DIR := $(addprefix build/,$(MODULES) $(TB_MODULES) $(DP_MODULES) $(OP_MODULES))

# INCLUDES  := $(addprefix -I,$(SRC_DIR))
INCLUDES  := -I./

vpath %.cpp $(SRC_DIR)
vpath %.cpp $(TB_SRC_DIR)
vpath %.cpp $(DP_SRC_DIR)
vpath %.cpp $(OP_SRC_DIR)

define make-goal
$1/%.o: %.cpp
	$(CC) $(INCLUDES) -c $$< -o $$@
endef

.PHONY: all checkdirs clean

all: checkdirs build/1D.b build/SSH.f build/SSH.b build/SSWF.b build/TBWF.b build/TBLB.b build/SSLB.b build/SSOP.b
# all: checkdirs build/SSH.f build/SSH.b build/SourceSinkDyn.b build/TB.b build/test.fermion build/test.boson build/test.lattice build/test.node

build/%.o: %.cpp
	$(CC) $(INCLUDES) -c $< -o $@

build/1D.b: build/1D_boson.o $(OBJ)
	$(CC) $^ -o $@ $(LAPACK_OMP) $(HDF5LIB)

build/SSH.b: build/SSH_boson.o $(OBJ)
	$(CC) $^ -o $@ $(LAPACK_OMP) $(HDF5LIB)

build/SSH.f: build/SSH_fermion.o $(OBJ)
	$(CC) $^ -o $@ $(LAPACK_OMP) $(HDF5LIB)

build/TBWF.b: build/TBWF_boson.o $(OBJ)
	$(CC) $^ -o $@ $(LAPACK_OMP) $(HDF5LIB)

build/TBLB.b: build/TBLB_boson.o $(TB_OBJ)
	$(CC) $^ -o $@ $(LAPACK_OMP) $(HDF5LIB)

build/SSWF.b: build/SSWF_boson.o $(OBJ)
	$(CC) $^ -o $@ $(LAPACK_OMP) $(HDF5LIB)

build/SSLB.b: build/SSLB_boson.o $(DP_OBJ)
	$(CC) $^ -o $@ $(LAPACK_OMP) $(HDF5LIB)

build/SSOP.b: build/SSOP_boson.o $(OP_OBJ)
	$(CC) $^ -o $@ $(LAPACK_OMP) $(HDF5LIB)

checkdirs: $(BUILD_DIR)

$(BUILD_DIR):
	@mkdir -p $@

clean:
	@rm -rf $(BUILD_DIR) build/*.o

$(foreach bdir,$(BUILD_DIR),$(eval $(call make-goal,$(bdir))))
