include ../../../make.inc

.SUFFIXES: .cpp

OBJ = compact_apwalm.o lapw.o lapw_init.o lapw_load_global.o lapw_load_species.o \
      lapw_seceqn.o lapw_set_h.o lapw_set_o.o
      
.cpp.o:
	$(CXX) $(CPP_OPTS) -Wall -O0 -g -c -o $(*D)/$(*F).o $<
all: $(OBJ)
	ar -r libapw.a $(OBJ)

clean:
	rm *.o
