CC      = mpicxx
PETSC_DIR=${HOME}/petsc-3.14.6
PETSC_ARCH=arch-linux-c-opt
PETSC=${PETSC_DIR}
PETSCLIB=$(PETSC)/$(PETSC_ARCH)/lib/
PETSCINC1=$(PETSC)/include/
PETSCINC2=$(PETSC)/$(PETSC_ARCH)/include/

#Serial HDF5
#HDF5INC=/usr/lib/x86_64-linux-gnu/hdf5/serial/include/
#HDF5LIB=/usr/lib/x86_64-linux-gnu/hdf5/serial/

#Parallel HDF5
HDF5INC=/usr/lib/x86_64-linux-gnu/hdf5/openmpi/include/
HDF5LIB=/usr/lib/x86_64-linux-gnu/hdf5/openmpi/

LIBFLAG =       -lpthread -lrt -ldl -lstdc++ \
		        -lpetsc -lf2cblas -lf2clapack -lhdf5 -lHYPRE

LIBDIR =        -L$(PETSCLIB) -L$(HDF5LIB)

CPPFLAGS =      -DNDEBUG -I$(PETSCINC1) -I$(PETSCINC2) -I$(HDF5INC) -O -g


SOURCEC =       main.c initialization.c io.c overset.c wallfunctions.c boundary.c \
		mesh.c inflow.c abl.c ueqn.c peqn.c les.c teqn.c turbines.c \
		initialField.c acquisition.c ibm.c clock.c precursor.c\

OBJSM   =       main.o initialization.o io.o overset.o wallfunctions.o boundary.o \
		mesh.o inflow.o abl.o ueqn.o peqn.o les.o teqn.o turbines.o \
		initialField.o acquisition.o ibm.o clock.o precursor.o\

OBJSW   =       windToPW.o initialization.o io.o wallfunctions.o overset.o boundary.o \
                mesh.o inflow.o abl.o ueqn.o peqn.o les.o teqn.o turbines.o \
                initialField.o acquisition.o ibm.o clock.o precursor.o\

tosca:   $(OBJSM)
		$(CC) -o ../bin/tosca $(OBJSM) $(CPPFLAGS) $(LIBDIR) $(LIBFLAG)

windToPW:   $(OBJSW)
		$(CC) -o ../bin/windToPW $(OBJSW) $(CPPFLAGS) $(LIBDIR) $(LIBFLAG)

clean:
	rm *.o ../bin/tosca ../bin/windToPW
