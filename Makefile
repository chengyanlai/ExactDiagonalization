include env.in

MODULES   := Node Lattice Basis Hamiltonian Lanczos hdf5io numeric
SRC_DIR   := $(addprefix src/,$(MODULES))
SRC       := $(foreach sdir,$(SRC_DIR),$(wildcard $(sdir)/*.cpp))
OBJ       := $(patsubst src/%.cpp,build/%.o,$(SRC))

PE_MODULES   := Lindblad-PE
PE_SRC_DIR   := $(addprefix src/,$(MODULES) $(PE_MODULES))
PE_SRC       := $(foreach sdir,$(PE_SRC_DIR),$(wildcard $(sdir)/*.cpp))
PE_OBJ       := $(patsubst src/%.cpp,build/%.o,$(PE_SRC))

BUILD_DIR := $(addprefix build/,$(MODULES) $(PE_MODULES) apps)

INCLUDES  := -I./

vpath %.cpp $(SRC_DIR)
vpath %.cpp $(PE_SRC_DIR)
vpath %.cpp apps/

define make-goal
$1/%.o: %.cpp
	$(CC) $(INCLUDES) -c $$< -o $$@
endef

.PHONY: all checkdirs clean

all: checkdirs build/loop.f build/xas.f build/plex build/1D.f build/rixs

mpi: checkdirs build/plex.mpi build/rixs.mpi

build/apps/%.o: apps/%.cpp
	$(CC) $(INCLUDES) -c $< -o $@

build/1D.f: build/apps/1D_fermion.o $(OBJ)
	$(CC) $^ -o $@ $(LAPACK) $(HDF5LIB)

build/xas.f: build/apps/TRXAS.o $(OBJ)
# build/xas.f: build/apps/XAS2.o $(OBJ)
	$(CC) $^ -o $@ $(LAPACK) $(HDF5LIB)

build/loop.f: build/apps/LoopFermion.o $(PE_OBJ)
	$(CC) $^ -o $@ $(LAPACK) $(HDF5LIB)

build/plex: build/apps/FermionBosonMix.o $(OBJ)
	$(CC) $^ -o $@ $(LAPACK) $(HDF5LIB)

build/rixs: build/apps/RIXS.o $(OBJ)
	$(CC) $^ -o $@ $(LAPACK) $(HDF5LIB)

build/apps/FermionBosonMixMPI.o: apps/FermionBosonMix.cpp
	$(MPICC) $(INCLUDES) -c $< -o $@

build/plex.mpi: build/apps/FermionBosonMixMPI.o $(OBJ)
	$(MPICC) $^ -o $@ $(LAPACK) $(HDF5LIB)

build/apps/RIXSMPI.o: apps/RIXS.cpp
	$(MPICC) $(INCLUDES) -c $< -o $@

build/rixs.mpi: build/apps/RIXSMPI.o $(OBJ)
	$(MPICC) $^ -o $@ $(LAPACK) $(HDF5LIB)

checkdirs: $(BUILD_DIR)

$(BUILD_DIR):
	@mkdir -p $@

clean:
	@rm -rf $(BUILD_DIR) build/*.o

$(foreach bdir,$(BUILD_DIR),$(eval $(call make-goal,$(bdir))))
