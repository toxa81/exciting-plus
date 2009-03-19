include ../make.inc

.SUFFIXES: .o .f90
.f90.o:
	$(F77) $(F77_OPTS) -c -o $(*F).o $<

all: writeinput

#SRC_writeinput = ../src/modmain.f90 
#writeinput.f90
#OBJ_writeinput = $(SRC_writeinput:.f90=.o) 
modmain.o:
	$(F77) $(F77_OPTS) -c ../src/modmain.f90
readinput.o:
	$(F77) $(F77_OPTS) -c ../src/readinput.f90
readspecies.o:
	$(F77) $(F77_OPTS) -c ../src/readspecies.f90
r3minv.o:
	$(F77) $(F77_OPTS) -c ../src/r3minv.f90
r3mv.o:
	$(F77) $(F77_OPTS) -c ../src/r3mv.f90
	
OBJ_writeinput =  modmain.o readinput.o readspecies.o r3minv.o r3mv.o writeinput.o 

writeinput: $(OBJ_writeinput)
	$(F77) -o $@ $(OBJ_writeinput)
		
clean:
	rm *.o *.mod