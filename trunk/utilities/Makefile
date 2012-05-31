include ../make.inc

.SUFFIXES: .o .f90
.f90.o:
	$(F90) $(F90_OPTS) -I../src -c -o $(*F).o $<

all: qdir 

qdir: qdir.o
	$(F90) -o $@ qdir.o  $(LIBS) ../src/libelk.a
		
clean:
	rm *.o *.mod