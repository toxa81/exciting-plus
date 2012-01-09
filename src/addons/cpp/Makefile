include ../../../make.inc

.SUFFIXES: .cpp

OBJ = compact_apwalm.o lapw.o lapw_init.o lapw_load_global.o lapw_load_species.o \
      lapw_seceqn.o lapw_set_h.o lapw_set_o.o lapw_fvmt.o lapw_test_fvmt.o lapw_set_sv.o \
      lapw_seceqn_fv.o lapw_seceqn_fv_old.o lapw_cmpwf.o lapw_cmpbf.o lapw_load_kpoint.o \
      lapw_seceqn_init.o        
.cpp.o:
	$(CXX) $(CPP_OPTS) -Wall -O0 -g -c -o $(*D)/$(*F).o $<
all: $(OBJ)
	ar -r libapw.a $(OBJ)

clean:
	rm *.o
