include ../../make.inc
F90_OPTS=-O3 -fast -fastsse -Mpreprocess -Munroll -Minfo=all -D_HDF5_ ${HDF5_FLIB}

.SUFFIXES: .o .f90
.f90.o:
	$(F90) $(F90_OPTS) -c -o $(*F).o $<

all: pp_u

pp_u: mod_hdf5.o pp_u.o
	$(F90) -o $@ mod_hdf5.o pp_u.o $(LIB_SYS)
		
clean:
	rm *.o *.mod