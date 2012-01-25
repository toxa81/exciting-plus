include ../../../../make.inc

.SUFFIXES: .cpp

OBJ = lapw.o lapw_init.o lapw_load_global.o lapw_load_species.o lapw_set_sv.o \
      lapw_load_kpoint.o lapw_seceqn_init.o lapw_fft.o lapw_eigen_states.o \
      timer.o lapw_timers.o lapw_execute.o
      
.cpp.o:
	$(CXX) $(CPP_OPTS) -Wall -O3 -c -o $(*D)/$(*F).o $<
all: $(OBJ)
	ar -r libapw.a $(OBJ)

clean:
	rm *.o
