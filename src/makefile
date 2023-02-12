# Set PETSC_DIR and PETSC_ARCH variables in the bash file 

CC        = mpicxx
PETSC     = ${PETSC_DIR}
PETSCLIB  = $(PETSC)/$(PETSC_ARCH)/lib/
PETSCINC1 = $(PETSC)/include/
PETSCINC2 = $(PETSC)/$(PETSC_ARCH)/include/

LIBFLAG   = -lpthread -lrt -ldl -lstdc++ \
		    -lpetsc -lf2cblas -lf2clapack -lHYPRE #-lcatalyst

LIBDIR    = -L$(PETSCLIB) -L$(HDF5LIB) -L$(CATALYST)/lib

CPPFLAGS  = -I$(PETSCINC1) -I$(PETSCINC2) -I$(CATALYST)/include/catalyst-2.0/ -O3 -g -lhdf5\
            -DPARAVIEW_IMPL_DIR=\"\" -DUSE_CATALYST=0


SOURCEC =       main.c initialization.c io.c overset.c wallfunctions.c boundary.c \
		mesh.c inflow.c abl.c ueqn.c peqn.c les.c teqn.c turbines.c \
		initialField.c acquisition.c ibm.c clock.c precursor.c ibmInput.c\

OBJSM   =       main.o initialization.o io.o overset.o wallfunctions.o boundary.o \
		mesh.o inflow.o abl.o ueqn.o peqn.o les.o teqn.o turbines.o \
		initialField.o acquisition.o ibm.o clock.o precursor.o ibmInput.o\

OBJSW   =       windToPW.o initialization.o io.o wallfunctions.o overset.o boundary.o \
                mesh.o inflow.o abl.o ueqn.o peqn.o les.o teqn.o turbines.o \
                initialField.o acquisition.o ibm.o clock.o precursor.o ibmInput.o\

tosca:   $(OBJSM)
		$(CC) -o tosca $(OBJSM) $(CPPFLAGS) $(LIBDIR) $(LIBFLAG)

windToPW:   $(OBJSW)
		$(CC) -o windToPW $(OBJSW) $(CPPFLAGS) $(LIBDIR) $(LIBFLAG)

clean:
	rm *.o tosca windToPW
