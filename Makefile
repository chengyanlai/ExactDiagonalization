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
	EIGENHOME = /Volumes/Files/PublicRepo/eigen
	HDF5HOME = /usr/local/opt/hdf5
	NumCores = 1
	HDF5LIB = -L$(HDF5HOME)/lib -lhdf5 -lhdf5_cpp
	LAPACK = -lblas -llapack -lm -larpack
	CC = clang++ -O3 -m64 -std=c++11 -stdlib=libc++ -Isrc/ -I$(HDF5HOME)/include -I$(EIGENHOME)
else
	ifeq ("$(NODE)", "kagome.ucmerced.edu")
		OPENMP = -qopenmp
		NumCores = 16
		ARPACK = -L$(ARPACKHOME)/lib -larpack
	else ifeq ("$(NODE)", "braid.cnsi.ucsb.edu")
		OPENMP = -openmp
		NumCores = 20
		ARPACK = -larpack
	else ifeq ("$(NODE)", "atomtronics.ucmerced.edu")
		OPENMP = -openmp
		NumCores = 4
		ARPACK = /home/chengyanlai/Downloads/arpack-ng/build/lib/libarpack.so
	else ifeq ("$(NODE)", "merced.cluster")
		OPENMP = -openmp
		NumCores = 20
		ARPACK = -L/home/clai24/apps/arpack-ng/lib -larpack
	else ifneq (, $(filter "$(NODE)", "comet-ln1.sdsc.edu" "comet-ln2.sdsc.edu" "comet-ln3.sdsc.edu" "comet-ln4.sdsc.edu"))
		OPENMP = -openmp
		NumCores = 24
		ARPACK = -larpack
	else ifneq (, $(filter "$(NODE)", "login1.stampede.tacc.utexas.edu" "login2.stampede.tacc.utexas.edu" "login3.stampede.tacc.utexas.edu" "login4.stampede.tacc.utexas.edu"))
		OPENMP = -openmp
		NumCores = 16
		ARPACK = -larpack
	endif
	# NOTE: the order of linker matters!
	HDF5LIB = $(HDF5HOME)/lib/libhdf5_cpp.a $(HDF5HOME)/lib/libhdf5.a $(HDF5HOME)/lib/libz.a $(HDF5HOME)/lib/libszip.a
	LAPACK = $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a \
	$(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group \
	${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
	${MKLROOT}/lib/intel64/libmkl_core.a \
	${MKLROOT}/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -lpthread -lm $(ARPACK)
	CC = icpc $(OPENMP) -DNumCores=$(NumCores) -O3 -Wall -std=c++11 -Isrc/ -I$(HDF5HOME)/include -I$(EIGENHOME) -DMKL
endif

MODULES   := Node Lattice Basis Hamiltonian Lanczos hdf5io
SRC_DIR   := $(addprefix src/,$(MODULES))
SRC       := $(foreach sdir,$(SRC_DIR),$(wildcard $(sdir)/*.cpp))
OBJ       := $(patsubst src/%.cpp,build/%.o,$(SRC))

TB_MODULES   := Lindblad-TB
TB_SRC_DIR   := $(addprefix src/,$(MODULES) $(TB_MODULES))
TB_SRC       := $(foreach sdir,$(TB_SRC_DIR),$(wildcard $(sdir)/*.cpp))
TB_OBJ       := $(patsubst src/%.cpp,build/%.o,$(TB_SRC))

DP_MODULES   := Lindblad-DP
DP_SRC_DIR   := $(addprefix src/,$(MODULES) $(DP_MODULES))
DP_SRC       := $(foreach sdir,$(DP_SRC_DIR),$(wildcard $(sdir)/*.cpp))
DP_OBJ       := $(patsubst src/%.cpp,build/%.o,$(DP_SRC))

OP_MODULES   := Lindblad-OP
OP_SRC_DIR   := $(addprefix src/,$(MODULES) $(OP_MODULES))
OP_SRC       := $(foreach sdir,$(OP_SRC_DIR),$(wildcard $(sdir)/*.cpp))
OP_OBJ       := $(patsubst src/%.cpp,build/%.o,$(OP_SRC))

SS_MODULES   := Lindblad-SS
SS_SRC_DIR   := $(addprefix src/,$(MODULES) $(SS_MODULES))
SS_SRC       := $(foreach sdir,$(SS_SRC_DIR),$(wildcard $(sdir)/*.cpp))
SS_OBJ       := $(patsubst src/%.cpp,build/%.o,$(SS_SRC))

BUILD_DIR := $(addprefix build/,$(MODULES) $(TB_MODULES) $(DP_MODULES) $(OP_MODULES) $(SS_MODULES) apps)

# INCLUDES  := $(addprefix -I,$(SRC_DIR))
INCLUDES  := -I./

vpath %.cpp $(SRC_DIR)
vpath %.cpp $(TB_SRC_DIR)
vpath %.cpp $(DP_SRC_DIR)
vpath %.cpp $(OP_SRC_DIR)
vpath %.cpp $(SS_SRC_DIR)
vpath %.cpp apps/

define make-goal
$1/%.o: %.cpp
	$(CC) $(INCLUDES) -c $$< -o $$@
endef

.PHONY: all checkdirs clean

all: checkdirs build/1D.b build/SSH.f build/SSH.b# build/SSWF.b build/TBLB.b build/SSLB.b build/SSOP.b build/SSt.b

build/apps/%.o: apps/%.cpp
	$(CC) $(INCLUDES) -c $< -o $@

build/1D.b: build/apps/1D_boson.o $(OBJ)
	$(CC) $^ -o $@ $(LAPACK) $(HDF5LIB)

build/SSH.b: build/apps/SSH_boson.o $(OBJ)
	$(CC) $^ -o $@ $(LAPACK) $(HDF5LIB)

build/SSH.f: build/apps/SSH_fermion.o $(OBJ)
	$(CC) $^ -o $@ $(LAPACK) $(HDF5LIB)

build/TBLB.b: build/apps/TBLB_boson.o $(TB_OBJ)
	$(CC) $^ -o $@ $(LAPACK) $(HDF5LIB)

build/SSWF.b: build/apps/SSWF_boson.o $(OBJ)
	$(CC) $^ -o $@ $(LAPACK) $(HDF5LIB)

build/SSLB.b: build/apps/SSLB_boson.o $(DP_OBJ)
	$(CC) $^ -o $@ $(LAPACK) $(HDF5LIB)

build/SSOP.b: build/apps/SSOP_boson.o $(OP_OBJ)
	$(CC) $^ -o $@ $(LAPACK) $(HDF5LIB)

build/SSt.b: build/apps/SSt_boson.o $(SS_OBJ)
	$(CC) $^ -o $@ $(LAPACK) $(HDF5LIB)

checkdirs: $(BUILD_DIR)

$(BUILD_DIR):
	@mkdir -p $@

clean:
	@rm -rf $(BUILD_DIR) build/*.o

$(foreach bdir,$(BUILD_DIR),$(eval $(call make-goal,$(bdir))))
