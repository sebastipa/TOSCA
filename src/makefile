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
USE_CATALYST ?= 0
# ==================================================================================================================================== #

# mpixx compiler should be available with OpenMPI (required)
CC         = mpicxx

# to be defined if PETSc is compiled from source (suggested)
PETSC      = ${PETSC_DIR}

# not required if hdf5 is loaded with e.g. environent-modules (suggested)
HDF5       = ${HDF5_DIR}

# required to build with paraview-catalyst off-screen rendering capabilities
CATALYST   = ${CATALYST_DIR}

# Update LIBS 
LIBS       = -lpthread \
             -lrt \
             -ldl \
             -lstdc++ \
             -lpetsc \
             -lf2cblas \
             -lf2clapack \
             -lHYPRE \
             -lhdf5

ifeq ($(USE_CATALYST), 1)
LIBS      += -lcatalyst
endif

# Update LDFLAGS 
LDFLAGS    = -L$(PETSC)/$(PETSC_ARCH)/lib \
             -L$(HDF5)/lib

ifeq ($(USE_CATALYST), 1)
LDFLAGS   += -L$(CATALYST)/lib
endif

# Update INC_FLAGS 
INC_FLAGS  = -I$(PETSC)/include \
             -I$(PETSC)/$(PETSC_ARCH)/include \
             -I$(HDF5)/include

ifeq ($(USE_CATALYST), 1)
INC_FLAGS += -I$(CATALYST)/include/catalyst-2.0
endif

# Update CPPFLAGS 
CPPFLAGS   = $(INC_FLAGS) -DPARAVIEW_IMPL_DIR=\"\"

ifeq ($(USE_CATALYST), 1)
CPPFLAGS  += -DUSE_CATALYST=1
else
CPPFLAGS  += -DUSE_CATALYST=0
endif

CFLAGS     = -O3 -g 

OBJS_TOSCA =  main.o initialization.o io.o overset.o wallfunctions.o boundary.o \
    mesh.o inflow.o abl.o ueqn.o peqn.o les.o teqn.o turbines.o \
    initialField.o acquisition.o ibm.o clock.o precursor.o ibmInput.o

OBJS_PP   =  tosca2PV.o initialization.o io.o wallfunctions.o overset.o boundary.o \
    mesh.o inflow.o abl.o ueqn.o peqn.o les.o teqn.o turbines.o \
    initialField.o acquisition.o ibm.o clock.o precursor.o ibmInput.o

tosca:  $(OBJS_TOSCA)
	$(CC) -o $@ $(OBJS_TOSCA) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $(LIBS)
	
tosca2PV:   $(OBJS_PP)
	$(CC) -o $@ $(OBJS_PP) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $(LIBS)
	
all:    tosca tosca2PV

clean:
	rm *.o tosca tosca2PV
