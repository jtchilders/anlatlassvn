# DO NOT EDIT FROM HERE ON:
#
# DEFINE DIRECTORY AND FILE ALIASES
alp= ../alplib
her= ../herlib
prclib= ../$(prc)lib
prcusr=.
prcfile=$(prclib)/$(prc)
ifeq ($(strip $(usrfile)),)
  execfile=$(prc)gen
  usrfile=$(prcusr)/$(prc)usr
else
  execfile=$(usrfile)
endif

# DEFINE FILE GROUPS:
# Files used for the parton-level event genertaion:
ALPGEN= $(alp)/Acp.f $(alp)/Aint.f $(alp)/alputi.f \
	$(alp)/alppdf.f $(alp)/Asu3.f $(alp)/Aint90.f 
PARTON= $(prcfile).f $(usrfile).f $(prclib)/Aproc.f $(ALPGEN)

# Include files
INC=  $(prclib)/$(prc).inc $(alp)/alpgen.inc

# include files' dependencies
$(PARTON): $(INC)

# object files
OBJ=$(PARTON:.f=.o) $(PARTON90:.f90=.o)


# compilation
%.o: %.f $(PARTON) $(INC)
	$(FFF) -c -o $*.o $*.f 
$(prclib)/XXX.o90 : $(alp)/A90.f90 $(prclib)/ini_$(prc).f90
	cd $(prclib); cp $(alp)/A90.f90 XXX.f90; \
	cat $(prclib)/ini_$(prc).f90 >> XXX.f90;\
	$(FF90)  -c XXX.f90; cp XXX.o XXX.o90

$(prclib)/XXX.o90V : $(alp)/A90.f90 $(prclib)/ini_$(prc).f90
	cd $(VF90); cp $(alp)/A90.f90 XXX.f90;\
        cat $(prclib)/ini_$(prc).f90 >> XXX.f90;\
	$(FF90V) -c XXX.f90; cp XXX.o $(prclib)/XXX.o90V; mv XXX.o $(prclib); \
	rm -f *.vo; rm -f V*.inc;

# fortran77 version
gen: $(OBJ)
	$(FFF) -c $(BGQ_SUFFIX) -o $(alp)/alpgen.o $(alp)/alpgen.f
	$(FFFL) -o $(execfile) $(usrfile).o $(prcfile).o \
	$(alp)/alpgen.o $(alp)/alputi.o $(alp)/alppdf.o \
	$(alp)/Acp.o $(alp)/Asu3.o $(alp)/Aint.o
# fortran90 version
gen90: $(usrfile).o $(prcfile).o $(prclib)/$(prc).inc \
	$(alp)/alputi.o $(alp)/alppdf.o \
	$(alp)/Aint90.o $(prclib)/XXX.o90 $(alp)/alpgen.inc
	$(FFF) -c $(BGQ_SUFFIX) -o $(alp)/alpgen.o $(alp)/alpgen.f
	$(FFL90) -o $(execfile)90 $(usrfile).o $(prcfile).o \
	$(alp)/alpgen.o $(alp)/alputi.o $(alp)/alppdf.o \
	$(alp)/Aint90.o $(prclib)/XXX.o
# fortran90 version, Vast/Veridian compyler
gen90V: $(usrfile).o $(prcfile).o $(prclib)/$(prc).inc\
	$(alp)/alputi.o $(alp)/alppdf.o \
	$(alp)/Aint90.o $(prclib)/XXX.o90V $(alp)/alpgen.inc
	$(FFF) -c -o $(alp)/alpgen.o $(alp)/alpgen.f
	$(FFFL) -o $(execfile)90V $(usrfile).o $(prcfile).o \
	$(alp)/alpgen.o $(alp)/alputi.o $(alp)/alppdf.o \
	$(prclib)/XXX.o $(alp)/Aint90.o $(VF90)/libvast90.a

# MPI VERSIONS
#
# fortran77 version
gen_mpi: $(OBJ)
	$(FFF) -c $(BGQ_SUFFIX) $(BGQ_CPP_DEFINE)-DUSE_MPI=1 -o $(alp)/alppdf.o $(alp)/alppdf.f
	$(FFF) -c $(BGQ_SUFFIX) $(BGQ_CPP_DEFINE)-DUSE_MPI=1 -o $(alp)/alpgen.o $(alp)/alpgen.f
	$(FFFL) -o $(execfile)_mpi $(usrfile).o $(prcfile).o \
	          $(alp)/alpgen.o $(alp)/alputi.o $(alp)/alppdf.o \
	          $(alp)/Acp.o $(alp)/Asu3.o $(alp)/Aint.o
# fortran90 version
gen90_mpi: $(usrfile).o $(prcfile).o $(prclib)/$(prc).inc\
	       $(alp)/alputi.o  \
	       $(alp)/Aint90.o $(prclib)/XXX.o90 $(alp)/alpgen.inc
	$(FFF) -c $(BGQ_SUFFIX) $(BGQ_CPP_DEFINE)-DUSE_MPI=1 -o $(alp)/alppdf.o $(alp)/alppdf.f
	$(FFF) -c $(BGQ_SUFFIX) $(BGQ_CPP_DEFINE)-DUSE_MPI=1 -o $(alp)/alpgen.o $(alp)/alpgen.f
	          $(FFL90) -o $(execfile)90_mpi $(usrfile).o $(prcfile).o \
	          $(alp)/alpgen.o $(alp)/alputi.o $(alp)/alppdf.o \
	          $(alp)/Aint90.o $(prclib)/XXX.o

# fortran 77 with MPI and no MRST PDFs
gen_mpi_nomrstpdfs: $(OBJ)
	$(FFF) -c $(BGQ_SUFFIX) $(BGQ_CPP_DEFINE)-DUSE_MPI=1 $(BGQ_CPP_DEFINE)-DDISABLE_MRST_PDFS -o $(alp)/alppdf.o $(alp)/alppdf.f
	$(FFF) -c $(BGQ_SUFFIX) $(BGQ_CPP_DEFINE)-DUSE_MPI=1  -o $(alp)/alpgen.o $(alp)/alpgen.f
	$(FFFL) -o $(execfile)_mpi_nomrstpdfs $(usrfile).o $(prcfile).o \
	          $(alp)/alpgen.o $(alp)/alputi.o $(alp)/alppdf.o \
	          $(alp)/Aint90.o $(prclib)/XXX.o

# fortran 90 with MPI and no MRST PDFs
gen90_mpi_nomrstpdfs: $(usrfile).o $(prcfile).o $(prclib)/$(prc).inc\
	$(alp)/alputi.o \
	$(alp)/Aint90.o $(prclib)/XXX.o90 $(alp)/alpgen.inc
	$(FFF) -c $(BGQ_SUFFIX) $(BGQ_CPP_DEFINE)-DUSE_MPI=1 $(BGQ_CPP_DEFINE)-DDISABLE_MRST_PDFS -o $(alp)/alppdf.o $(alp)/alppdf.f
	$(FFF) -c $(BGQ_SUFFIX) $(BGQ_CPP_DEFINE)-DUSE_MPI=1  -o $(alp)/alpgen.o $(alp)/alpgen.f
	$(FFL90) -o $(execfile)90_mpi_nomrstpdfs $(usrfile).o $(prcfile).o \
	$(alp)/alpgen.o $(alp)/alputi.o $(alp)/alppdf.o \
	$(alp)/Aint90.o $(prclib)/XXX.o

# DIRECTORY CLEANUP UTILITIES:
#
# remove object files only
cleanobj:
	-rm -f $(PARTON:.f=.o) $(PARTON90:.f90=.o) $(prcusr)/../*/*.o90* $(alp)/alpgen.o

# remove object files, etc
cleanall:
	-rm -f $(OBJ) $(prcusr)/fort.* $(prcusr)/*.top $(prcusr)/*.par \
	$(prcusr)/*.wgt $(prcusr)/*.unw $(prcusr)/*.mon \
	$(prcusr)/*.stat $(prcusr)/../*/*.o90* $(alp)/alpgen.o

