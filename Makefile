include env.in

MODULES   := Node Lattice Basis Hamiltonian Lanczos hdf5io numeric
SRC_DIR   := $(addprefix src/,$(MODULES))
SRC       := $(foreach sdir,$(SRC_DIR),$(wildcard $(sdir)/*.cpp))
OBJ       := $(patsubst src/%.cpp,build/%.o,$(SRC))

FHM_MODULES   := Hamiltonian/FHM
FHM_SRC_DIR   := $(addprefix src/,$(MODULES) $(FHM_MODULES))
FHM_SRC       := $(foreach sdir,$(FHM_SRC_DIR), $(wildcard $(sdir)/*.cpp))
FHM_OBJ       := $(patsubst src/%.cpp, build/%.o,$(FHM_SRC))

BHM_MODULES   := Hamiltonian/BHM
BHM_SRC_DIR   := $(addprefix src/,$(MODULES) $(BHM_MODULES))
BHM_SRC       := $(foreach sdir,$(BHM_SRC_DIR), $(wildcard $(sdir)/*.cpp))
BHM_OBJ       := $(patsubst src/%.cpp, build/%.o,$(BHM_SRC))

Holstein_MODULES   := Hamiltonian/Holstein
Holstein_SRC_DIR   := $(addprefix src/,$(MODULES) $(Holstein_MODULES))
Holstein_SRC       := $(foreach sdir,$(Holstein_SRC_DIR), $(wildcard $(sdir)/*.cpp))
Holstein_OBJ       := $(patsubst src/%.cpp, build/%.o,$(Holstein_SRC))

PE_MODULES   := Lindblad/ParticleExchange
PE_SRC_DIR   := $(addprefix src/,$(MODULES) $(PE_MODULES))
PE_SRC       := $(foreach sdir,$(PE_SRC_DIR), $(wildcard $(sdir)/*.cpp))
PE_OBJ       := $(patsubst src/%.cpp, build/%.o,$(PE_SRC))

BUILD_DIR := $(addprefix build/,$(MODULES) $(FHM_MODULES) $(BHM_MODULES) $(Holstein_MODULES) $(PE_MODULES) apps)

INCLUDES  := -I./

vpath %.cpp $(SRC_DIR)
vpath %.cpp $(PE_SRC_DIR)
vpath %.cpp apps/

define make-goal
$1/%.o: %.cpp
	$(CC) $(INCLUDES) -c $$< -o $$@
endef

.PHONY: all checkdirs clean

all: checkdirs build/fhm.1d build/bhm.1d build/tqdm build/holstein.1dInfty

mpi: checkdirs build/plex.mpi build/rixs.mpi

build/apps/%.o: apps/%.cpp
	$(CC) $(INCLUDES) -c $< -o $@

build/fhm.1d: build/apps/FHM1D.o $(OBJ) $(FHM_OBJ)
	$(CC) $^ -o $@ $(ARMADILLO) $(LAPACK) $(ARPACK) $(HDF5LIB)

build/bhm.1d: build/apps/BHM1D.o $(OBJ) $(BHM_OBJ)
	$(CC) $^ -o $@ $(ARMADILLO) $(LAPACK) $(ARPACK) $(HDF5LIB)

build/tqdm: build/apps/TQDM.o $(OBJ) $(FHM_OBJ) $(PE_OBJ)
	$(CC) $^ -o $@ $(ARMADILLO) $(LAPACK) $(ARPACK) $(HDF5LIB)

build/holstein.1dInfty: build/apps/Holstein1DInfty.o $(OBJ) $(Holstein_OBJ)
	$(CC) $^ -o $@ $(ARMADILLO) $(LAPACK) $(ARPACK) $(HDF5LIB)

checkdirs: $(BUILD_DIR)

$(BUILD_DIR):
	@mkdir -p $@

clean:
	@rm -rf $(BUILD_DIR) build/*.o

$(foreach bdir,$(BUILD_DIR),$(eval $(call make-goal,$(bdir))))
