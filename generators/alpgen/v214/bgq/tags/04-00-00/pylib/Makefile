# DEFINE FORTRAN COMPILATION DEFAULTS
include ../compile.mk
#
# DEFINE DIRECTORY TREE STRUCTURE:
py= .
alp=../alplib

LIBPY = libpythia.a
# DEFINE FILE GROUPS:
#

py_version=6428
#py_version=6325

# Files used for the shower evolution
SHOWER= $(py)/*.f $(alp)/alpsho.f

# Include files
INC= $(py)/hepevt.inc $(py)/hnt.inc $(py)/pythia6300.inc

# object files
#OBJ=$(SHOWER:.f=.o) 


# compilation
%.o: %.f
	$(FFF) -c -o $*.o $*.f \

# LINKING GENERIC EXECUTABLE FOR PYTHIA EVENT-EVOLUTION
pyuser: $(py)/pyuser.o $(py)/atopyt.o \
        $(py)/pythia$(py_version).o  alpsho.o
	$(FFF) -o $(py)/pyuser $(py)/pyuser.o \
        $(py)/atopyt.o $(py)/pythia$(py_version).o $(py)/alpsho.o -L/users/hpcusers/svn/tools/lhapdf-5.9.1/local/lib -lLHAPDF
pydbg: $(py)/pydbg.o $(py)/atopyt.o \
        $(py)/pythia$(py_version).o  $(alp)/alpsho-dbg.o $(alp)/alputi.o 
	$(FFF) -o $(py)/pydbg $(py)/pydbg.o \
        $(py)/atopyt.o $(py)/pythia$(py_version).o $(alp)/alpsho-dbg.o $(alp)/alputi.o 
pyben: $(py)/pyben.o $(py)/atopyt.o \
        $(py)/pythia$(py_version).o  $(alp)/alpsho.o $(alp)/alputi.o 
	$(FFF) -o $(py)/pyben $(py)/pyben.o \
        $(py)/atopyt.o $(py)/pythia$(py_version).o $(alp)/alpsho.o $(alp)/alputi.o 
# $(py)/ntuple.o  
#	`cernlib mathlib`  

# DIRECTORY CLEANUP UTILITIES:
#
# remove object files only
cleanobj:
	-rm $(OBJ) $(py)/*.o  $(alp)/alpsho.o

# remove object files, etc
cleanall:
	-rm $(OBJ) $(py)/fort.* $(py)/*.top $(py)/*.par \
	$(py)/*.wgt $(py)/*.unw $(py)/*.mon \
	$(py)/*.stat $(alp)/alpsho.o
