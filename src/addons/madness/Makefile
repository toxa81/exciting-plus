include ../../../make.inc

OBJ = madness_intrf.o

.cpp.o:
	$(CXX) $(CXX_OPTS) $(MADNESS_INC) -std=c++0x -c $<

all: $(OBJ)
	ar -r libmad.a $(OBJ)

clean:
	rm *.o *.a
