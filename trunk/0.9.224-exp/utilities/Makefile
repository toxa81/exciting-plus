include ../make.inc

.SUFFIXES: .o .f90
.f90.o:
	$(F90) $(F77_OPTS) -c -o $(*F).o $<

all: writeinput qdir

writeinput: writeinput.o
	$(F90) -o $@ writeinput.o exciting.a
	
qdir: qdir.o
	$(F90) -o $@ qdir.o exciting.a fftlib.a -llapack -lblas

clean:
	rm *.o 