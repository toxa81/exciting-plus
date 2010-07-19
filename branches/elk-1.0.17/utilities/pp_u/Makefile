include ../../make.inc
F90_OPTS=-O2 -Wall -cpp -D_HDF5_ -I/Users/anton/local/include

.SUFFIXES: .o .f90
.f90.o:
	$(F90) $(F90_OPTS) -c -o $(*F).o $<

all: pp_u

pp_u: mod_hdf5.o pp_u.o
	$(F90) -o $@ mod_hdf5.o pp_u.o $(LIB_SYS)
		
clean:
	rm *.o *.mod