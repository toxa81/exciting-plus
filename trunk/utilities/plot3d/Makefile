include ../../make.inc
F90_OPTS_SERIAL=$(subst -D_MPI_,,$(F90_OPTS))

.SUFFIXES: .o .f90
.f90.o:
	$(F90) $(F90_OPTS_SERIAL) -c -o $(*F).o $<

all: sicvlm plot_wan_dens

mod_hdf5.o:
	$(F90) -c $(F90_OPTS_SERIAL) ../../src/addons/mod_hdf5.f90 -o mod_hdf5.o

genrlm.o:
	$(F90) -c $(F90_OPTS_SERIAL) ../../src/genrlm.f90 -o genrlm.o

genylm.o:
	$(F90) -c $(F90_OPTS_SERIAL) ../../src/genylm.f90 -o genylm.o

sphcrd.o:
	$(F90) -c $(F90_OPTS_SERIAL) ../../src/sphcrd.f90 -o sphcrd.o

sicvlm: mod_hdf5.o sicvlm.o genrlm.o genylm.o sphcrd.o
	$(F90) -o $@ mod_hdf5.o sicvlm.o genrlm.o genylm.o sphcrd.o $(LIBS)

plot_wan_dens: mod_hdf5.o plot_wan_dens.o genrlm.o genylm.o sphcrd.o
	$(F90) -o $@ mod_hdf5.o plot_wan_dens.o genrlm.o genylm.o sphcrd.o $(LIBS)
		
clean:
	rm *.o *.mod