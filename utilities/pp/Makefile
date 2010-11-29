include ../../make.inc

.SUFFIXES: .o .f90
.f90.o:
	$(F90) $(F90_OPTS) -c -o $(*F).o $<

all: pdos bndchr

pdos: addorb.o pdos.o
	$(F90) -o $@ addorb.o pdos.o
bndchr: addorb.o bndchr.o
	$(F90) -o $@ addorb.o bndchr.o
		
clean:
	rm *.o *.mod