import os
import shutil
import subprocess

first_revision = 1625
last_revision = 1630

make_inc="""
MAKE = make
F90 = mpif90
CXX = mpicxx
CC = mpicc
CPP_OPTS = -D_MPI_ -D_LIBAPW_
  
F90_OPTS = -O3 -Wall -cpp $(CPP_OPTS) -fopenmp -I/Users/anton/local/include 
F90_LINK_OPTS=-fopenmp -lstdc++

#LAPACK_LIB = -llapack -lblas
LAPACK_LIB = $(HOME)/local/lib/liblapack.a $(HOME)/local/lib/libblas.a

LIBAPW = ./addons/cpp/libapw.a

# === collect all libraries under one name ===
LIBS = $(LAPACK_LIB) $(HDF5_LIB) $(XC_LIB) $(NFFT_LIB) $(MADNESS_LIB) $(LIBAPW)
"""

elk_in="""
tasks
  0

nempty
  20
  
maxscl
  100

avec
  1.0  1.0 -1.0
  1.0 -1.0  1.0
 -1.0  1.0  1.0

scale
  2.708

atoms
  1                                   : nspecies
  'Fe.in'                             : spfname
  1                                   : natoms
  0.0  0.0  0.0    0.0  0.0  0.1      : atposl, bfcmt

ngridk
  4  4  4
"""

Fe_in="""'Fe'                                       : spsymb
 'iron'                                     : spname
  -26.0000                                  : spzn
   101799.2074                              : spmass
  0.392232E-06    2.0000   32.8043   750    : sprmin, rmt, sprmax, nrmt
  10                                        : spnst
   1   0   1   2.00000    T                 : spn, spl, spk, spocc, spcore
   2   0   1   2.00000    T
   2   1   1   2.00000    T
   2   1   2   4.00000    T
   3   0   1   2.00000    F
   3   1   1   2.00000    F
   3   1   2   4.00000    F
   3   2   2   4.00000    F
   3   2   3   2.00000    F
   4   0   1   2.00000    F
   1                                        : apword
  0.1500   0  F                             : apwe0, apwdm, apwve
   4                                        : nlx
   0   1                                    : l, apword
  0.1500   0  T                             : apwe0, apwdm, apwve
   1   1                                    : l, apword
  0.1500   0  T                             : apwe0, apwdm, apwve
   2   1                                    : l, apword
  0.1500   0  T                             : apwe0, apwdm, apwve
   3   1                                    : l, apword
  0.1500   0  T                             : apwe0, apwdm, apwve
   6                                        : nlorb
   0   2                                    : lorbl, lorbord
  0.1500   0  T                             : lorbe0, lorbdm, lorbve
  0.1500   1  T
   1   2                                    : lorbl, lorbord
  0.1500   0  T                             : lorbe0, lorbdm, lorbve
  0.1500   1  T
   2   2                                    : lorbl, lorbord
  0.1500   0  T                             : lorbe0, lorbdm, lorbve
  0.1500   1  T
   3   2                                    : lorbl, lorbord
  0.1500   0  T                             : lorbe0, lorbdm, lorbve
  0.1500   1  T
   0   3                                    : lorbl, lorbord
  0.1500   0  F                             : lorbe0, lorbdm, lorbve
  0.1500   1  F
 -3.4344   0  T
   1   3                                    : lorbl, lorbord
  0.1500   0  F                             : lorbe0, lorbdm, lorbve
  0.1500   1  F
 -2.1817   0  T"""


#def initial_checkout():
#    shutil.rmtree("trunk-tmp", 1)   
#    os.system("svn checkout --revision " + str(first_revision) + " http://exciting-plus.googlecode.com/svn/trunk/ trunk-tmp ") 

def checkout(rev):
    subprocess.call(["svn","checkout","--revision",str(rev),"http://exciting-plus.googlecode.com/svn/trunk/","trunk-tmp"]) 

def add_make_inc():
    fout=open ("trunk-tmp/make.inc", "w")
    fout.write(make_inc)
    fout.close()
    
def make():
    add_make_inc()
    subprocess.call(["make","-C","./trunk-tmp","clean"])
    subprocess.call(["make","-C","./trunk-tmp"])
    return os.path.isfile("./trunk-tmp/src/elk")

def prepare_input(spinpol, exactrho):
    shutil.rmtree("run-tmp", 1)   
    os.mkdir("run-tmp")
    fout=open ("run-tmp/elk.in","w")
    fout.write(elk_in)
    if spinpol:
        fout.write("spinpol\n.true.\n")
    if exactrho:
        fout.write("exactrho\n.true.\n")
    fout.close()
    fout=open ("run-tmp/Fe.in","w")
    fout.write(Fe_in)
    fout.close()
    
def execute():
    wd = os.getcwd() + "/run-tmp"
    ex = os.getcwd() + "/trunk-tmp/src/elk"
    p = subprocess.Popen(ex,cwd=wd)
    p.wait()
    
def get_results(results, testid):
    fin = open("./run-tmp/TOTENERGY.OUT","r")
    lines = fin.readlines()
    fin.close()
    results[testid] = lines[-1]

def run_elk_default_nm(results):
    prepare_input(False, False)
    execute()
    get_results(results, "default_nm")

def run_elk_exactrho_nm(results):
    prepare_input(False, True)
    execute()
    get_results(results, "exactrho_nm")

def run_elk_default_mag(results):
    prepare_input(True, False)
    execute()
    get_results(results, "default_mag")

def run_elk_exactrho_mag(results):
    prepare_input(True, True)
    execute()
    get_results(results, "exactrho_mag")

def run_tests(rout):
    results = {}
    run_elk_default_nm(results)
    #run_elk_exactrho_nm(results)
    run_elk_default_mag(results)
    #run_elk_exactrho_mag(results)
    rout.write("    default_nm : " + results["default_nm"] + "\n") 
    #rout.write("    exactrho_nm : " + results["exactrho_nm"] + "\n") 
    rout.write("    default_mag : " + results["default_mag"] + "\n") 
    #rout.write("    exactrho_mag : " + results["exactrho_mag"] + "\n") 

def all_clean():
    #shutil.rmtree("trunk-tmp", 1) 
    #shutil.rmtree("run-tmp", 1)       
    try:
        os.remove("elk_results.txt")
    except:
        pass
        

all_clean()
for r in range(first_revision, last_revision + 1):
    rout = open("elk_results.txt","a+")
    rout.write("revision : " + str(r) + "\n")
    checkout(r)
    if (make()):
        run_tests(rout)
    else:
        rout.write("    compilation error\n")
    rout.write("\n")
    rout.close()
    
    


