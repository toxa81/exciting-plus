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
OBJ_writeinput = modmain.o

writeinput: $(OBJ_writeinput)
	$(F77) -o $@ $(OBJ_writeinput)
		
