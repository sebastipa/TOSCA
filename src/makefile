# Required list of packages needed to compile TOSCA
# ==================================================================================================================================== #
# 1. A working compiler (only tested with gcc)
# 2. OpenMPI   (we suggest using the default e.g. from environment-modules)
# 3. PETSc     (we suggest compiling from source, see User-Guide at https://sebastipa.github.io/TOSCA/tosca_installation.html)
# 4. HDF5      (we suggest using the PETSc automatic installation)
# 5. HYPRE     (we suggest using the PETSc automatic installation)
# 6. Catalyst2 (optional, see User-Guide at https://sebastipa.github.io/TOSCA/tosca_installation.html#paraview-catalyst2-os-rendering)
# ==================================================================================================================================== #

# User defined options 
# ==================================================================================================================================== #
USE_HDF5     ?= 1
USE_CATALYST ?= 0    
# ==================================================================================================================================== #

# Environment variables should be defined in the shell before running make
HDF5       = ${HDF5_DIR}
CATALYST   = ${CATALYST_DIR}

# Compiler and flags
CC         = mpicxx
CFLAGS     = -O3 -g
CPPFLAGS   = $(INC_FLAGS) -DPARAVIEW_IMPL_DIR=\"\"
LDFLAGS    = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib 
LIBS       = -lpthread -lrt -ldl -lstdc++ -lpetsc -lf2cblas -lf2clapack -lHYPRE 

# HDF5 support
ifeq ($(USE_HDF5), 1)
CPPFLAGS  += -DUSE_HDF5=1
LDFLAGS   += -L$(HDF5)/lib
INC_FLAGS += -I$(HDF5)/include
LIBS      += -lhdf5 
else
CPPFLAGS  += -DUSE_HDF5=0
endif

# Catalyst support
ifeq ($(USE_CATALYST), 1)
CPPFLAGS  += -DUSE_CATALYST=1
LDFLAGS   += -L$(CATALYST)/lib
INC_FLAGS += -I$(CATALYST)/include/catalyst-2.0
LIBS      += -lcatalyst
else
CPPFLAGS  += -DUSE_CATALYST=0
endif

# Include paths
INC_FLAGS  = -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include 

# Directories
BUILD_DIR  = build
SRC_DIR    = .
OBJ_DIR    = $(BUILD_DIR)

# Object files
OBJS_TOSCA = $(OBJ_DIR)/main.o $(OBJ_DIR)/initialization.o $(OBJ_DIR)/io.o $(OBJ_DIR)/overset.o \
             $(OBJ_DIR)/wallfunctions.o $(OBJ_DIR)/boundary.o $(OBJ_DIR)/mesh.o $(OBJ_DIR)/inflow.o \
             $(OBJ_DIR)/abl.o $(OBJ_DIR)/ueqn.o $(OBJ_DIR)/peqn.o $(OBJ_DIR)/les.o $(OBJ_DIR)/teqn.o \
             $(OBJ_DIR)/turbines.o $(OBJ_DIR)/initialField.o $(OBJ_DIR)/acquisition.o $(OBJ_DIR)/ibm.o \
             $(OBJ_DIR)/clock.o $(OBJ_DIR)/precursor.o $(OBJ_DIR)/ibmInput.o

OBJS_PP    = $(OBJ_DIR)/tosca2PV.o $(OBJ_DIR)/initialization.o $(OBJ_DIR)/io.o $(OBJ_DIR)/wallfunctions.o \
             $(OBJ_DIR)/overset.o $(OBJ_DIR)/boundary.o $(OBJ_DIR)/mesh.o $(OBJ_DIR)/inflow.o \
             $(OBJ_DIR)/abl.o $(OBJ_DIR)/ueqn.o $(OBJ_DIR)/peqn.o $(OBJ_DIR)/les.o $(OBJ_DIR)/teqn.o \
             $(OBJ_DIR)/turbines.o $(OBJ_DIR)/initialField.o $(OBJ_DIR)/acquisition.o $(OBJ_DIR)/ibm.o \
             $(OBJ_DIR)/clock.o $(OBJ_DIR)/precursor.o $(OBJ_DIR)/ibmInput.o

# Targets
all: $(BUILD_DIR) tosca tosca2PV

tosca: $(OBJS_TOSCA)
	$(CC) -o $@ $^ $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $(LIBS)

tosca2PV: $(OBJS_PP)
	$(CC) -o $@ $^ $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $(LIBS)

# Pattern rule for object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(OBJ_DIR)
	$(CC) -c $< -o $@ $(CFLAGS) $(CPPFLAGS)

# Clean targets
clean:
	rm -rf $(BUILD_DIR) tosca tosca2PV

# Help target
help:
	@echo "Available targets:"
	@echo "  all         - Build all targets (default, HDF5 required)"
	@echo "  tosca       - Build the TOSCA executable (HDF5 not required)"
	@echo "  tosca2PV    - Build the TOSCA2PV executable (HDF5 required)"
	@echo "  clean       - Remove object files and executables"
	@echo "  help        - Show this help message"

# Build directory
$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)