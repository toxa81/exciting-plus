
#-------------------------------------------------------------------------------
# defaults
#-------------------------------------------------------------------------------
MAKE = make
AR = ar

include ../make.inc

#-------------------------------------------------------------------------------
# Source files
#-------------------------------------------------------------------------------

# modules
SRC_modules = ./addons/mod_mpi_grid.f90 \
              ./addons/mod_timer.f90 \
              ./addons/mod_addons.f90 \
              ./addons/mod_addons_q.f90 \
              ./addons/mod_papi.f90 \
              ./addons/mod_hdf5.f90 \
              ./modmain.f90 \
              ./addons/wann/mod_wannier.f90 \
              ./addons/expigqr/mod_expigqr.f90 \
              ./addons/linresp/mod_linresp.f90 \
              ./addons/mod_nrkp.f90 \
              ./addons/mod_madness.f90 \
              ./modqpt.f90 \
              ./libxc_funcs.f90 \
              ./libxc.f90 \
              ./libxcifc.f90 \
              ./modxcifc.f90 \
              ./modldapu.f90 \
              ./modrdm.f90 \
              ./modphonon.f90 \
              ./modtest.f90 \
              ./addons/sic/mod_sic.f90
              
OBJ_modules = $(SRC_modules:.f90=.o)
SRC = $(shell python ../gensrclog.py; mv srclog.f90 ./addons; find . -name \*.f90)
OBJ = $(SRC:.f90=.o)

EXE = elk

#-------------------------------------------------------------------------------
# Suffix rules
#-------------------------------------------------------------------------------
.SUFFIXES: .f90

.f90.o:
	$(F90) $(F90_OPTS) -c -o $(*D)/$(*F).o $<

#-------------------------------------------------------------------------------
# Targets
#-------------------------------------------------------------------------------
all: gensrc elk

serial: F90_OPTS:=$(subst -D_MPI_,,$(F90_OPTS)) 
serial: gensrc elk lib

madelk: F90_OPTS:=$(F90_OPTS) -D_MAD_
madelk: gensrc madapi elk1

madapi:
	cd addons/madness; make

gensrc:
	python ../gensrclog.py
	mv srclog.f90 ./addons
	python ../ftemplate.py <./addons/mod_mpi_grid.tpl >./addons/mod_mpi_grid.f90
	
elk:	$(OBJ_modules) $(OBJ)
	$(F90) $(F90_OPTS) -o $(EXE) $(OBJ) $(LIB_XC) $(LIB_SYS) $(LIB_LAPACK)

elk1:	$(OBJ_modules) $(OBJ)
	$(F90) $(F90_OPTS) -o $(EXE) $(OBJ) addons/madness/madness_intrf.o $(MADNESS_LIB) $(LIB_XC) \
	$(LIB_SYS) $(LIB_LAPACK)

lib:
	ar r libelk.a $(OBJ)
	
clean:
	find . -name \*.o -exec rm {} \;
	rm -f *.mod *~ fort.* ifc* *.gcno gmon.out *.aux *.dvi *.log \
        elk.pdf *.tex *.toc elk.lst elk.out $(EXE)

doc:
	rm -f elk.aux elk.bbl elk.blg elk.log elk.pdf elk.tex elk.toc elk.lst
	ls $(SRC) > elk.lst
	./protex -F -s $(SRC_main) $$(cat elk.lst) > elk.tex
	pdflatex elk;pdflatex elk;pdflatex elk

backup:
	tar -czf elk.tgz $(SRC) fftlib eos spacegroup junk BLAS LAPACK \
         Makefile notes.txt docs ../examples ../tests protex COPYING \
         README elk_silhouette.pdf depend ../release gentest genalltests \
         checklist.txt rmspaces

VERSION = $$(awk -F"/" '/data version/ {print $$2}' modmain.f90 | sed 's/ //g;s/,/./g')$
RELEASE = ../release/elk

release:
	rm -rf $(RELEASE)
	mkdir $(RELEASE)
	cd fftlib; $(MAKE) clean
	cd eos; $(MAKE) clean
	cd spacegroup; $(MAKE) clean
	mkdir $(RELEASE)/src
	cp -r fftlib $(RELEASE)/src
	cp -r eos $(RELEASE)/src
	cp -r spacegroup $(RELEASE)/src
	cp $(SRC) Makefile protex elk_silhouette.pdf $(RELEASE)/src
	cp libxc.inc libxc_funcs.f90 libxc.f90 libxcifc.f90 $(RELEASE)/src
	cd ../examples; find . -type f -name *~ -exec rm -f {} \;
	cp -r ../examples $(RELEASE)
	cd ../tests; find . -type f -name *~ -exec rm -f {} \;
	cp -r ../tests $(RELEASE)
	mkdir $(RELEASE)/species
	cp ../species/*.in $(RELEASE)/species
	cp -r ../utilities $(RELEASE)
	cp ../COPYING $(RELEASE)
	cp ../README $(RELEASE)
	cp release_notes.txt $(RELEASE)
	cp ../Makefile $(RELEASE)
	cp ../setup $(RELEASE)
	$(MAKE) doc
	cp elk.pdf ../docs
	cd spacegroup;$(MAKE) doc;cp spacegroup.pdf ../../docs;$(MAKE) clean
	cp -r ../docs $(RELEASE)
	tar -C ../release -czf ../release/elk-$(VERSION).tgz elk
	cp ../docs/elk.pdf ../release
	cp ../docs/spacegroup.pdf ../release

lines:
	cat $(SRC) | wc -l

spaces:
	./rmspaces $(SRC)

# DO NOT DELETE
