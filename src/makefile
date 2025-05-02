# Required list of packages needed to compile TOSCA
# ========================================================================== #
# 1. A working compiler (only tested with gcc)
# 2. OpenMPI   (we suggest using the default e.g. from environment-modules)
# 3. PETSc     (we suggest compiling from source, see README file)
# 4. HDF5      (we suggest using the default e.g. from environment-modules)
# 5. Catalyst2 (optional, see README file)
# ========================================================================== #

# mpixx compiler should be available with OpenMPI (required)
CC         = mpicxx

# to be defined if PETSc is compiled from source (suggested)
PETSC      = ${PETSC_DIR}

# not required if hdf5 is loaded with e.g. environent-modules (suggested)
HDF5       = ${HDF5_DIR}

# required to build with paraview-catalyst off-screen rendering capabilities
CATALYST   = ${CATALYST_DIR}

# NOTEs:
# 1)  -lhdf5     is only needed to compile tosca2PV (post-processor). Remove if only 'tosca' target is invoked
# 2)  -lcatalyst is only needed if paraview-catalyst off-screen rendering capabilities are enabled with -DUSE_CATALYST=1                .

LIBS       = -lpthread \
             -lrt \
             -ldl \
             -lstdc++ \
		     -lpetsc \
		     -lf2cblas \
		     -lf2clapack \
		     -lHYPRE \
		     -lhdf5 \
#		     -lcatalyst

LDFLAGS    = -L$(PETSC)/$(PETSC_ARCH)/lib\
             -L$(HDF5)/lib \
             -L$(CATALYST)/lib \
#			 -L$(SLEPC)/lib \
			 -L$(SLEPC)/arch-opt/lib

INC_FLAGS  = -I$(PETSC)/include \
             -I$(PETSC)/$(PETSC_ARCH)/include \
             -I$(HDF5)/include \
             -I$(CATALYST)/include/catalyst-2.0 \
	         -I$(EIGEN_DIR)
			#-I/usr/include/eigen3
#			 -I$(SLEPC)/include \
			 -I$(SLEPC)/arch-opt/include

CFLAGS     = -O3 -g

CPPFLAGS   =  $(INC_FLAGS) -DPARAVIEW_IMPL_DIR=\"\" -DUSE_CATALYST=0

OBJS_TOSCA =  main.o initialization.o io.o overset.o wallfunctions.o boundary.o \
    mesh.o inflow.o abl.o ueqn.o peqn.o les.o teqn.o turbines.o \
    initialField.o acquisition.o ibm.o clock.o precursor.o ibmInput.o vents.o sm.o

OBJS_PP   =  tosca2PV.o initialization.o io.o wallfunctions.o overset.o boundary.o \
    mesh.o inflow.o abl.o ueqn.o peqn.o les.o teqn.o turbines.o \
    initialField.o acquisition.o ibm.o clock.o precursor.o ibmInput.o vents.o sm.o

tosva:  $(OBJS_TOSCA)
	$(CC) -o $@ $(OBJS_TOSCA) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $(LIBS)

tosva2PV:   $(OBJS_PP)
	$(CC) -o $@ $(OBJS_PP) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $(LIBS)

all:    tosva tosva2PV

clean:
	rm *.o tosva tosva2PV
