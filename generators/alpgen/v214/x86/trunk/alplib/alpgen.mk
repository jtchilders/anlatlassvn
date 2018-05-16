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
ALPGEN= $(alp)/alpgen.f $(alp)/Acp.f $(alp)/Aint.f $(alp)/alputi.f \
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
	$(FFF) -o $(execfile) $(usrfile).o $(prcfile).o \
	$(alp)/alpgen.o $(alp)/alputi.o $(alp)/alppdf.o \
	$(alp)/Acp.o $(alp)/Asu3.o $(alp)/Aint.o
# fortran90 version
gen90: $(usrfile).o $(prcfile).o $(prclib)/$(prc).inc\
	$(alp)/alpgen.o $(alp)/alputi.o $(alp)/alppdf.o \
	$(alp)/Aint90.o $(prclib)/XXX.o90 $(alp)/alpgen.inc 
	$(FF90) -o $(execfile)90 $(usrfile).o $(prcfile).o \
	$(alp)/alpgen.o $(alp)/alputi.o $(alp)/alppdf.o \
	$(alp)/Aint90.o $(prclib)/XXX.o
# fortran90 version, Vast/Veridian compyler
gen90V: $(usrfile).o $(prcfile).o $(prclib)/$(prc).inc\
	$(alp)/alpgen.o $(alp)/alputi.o $(alp)/alppdf.o \
	$(alp)/Aint90.o $(prclib)/XXX.o90V $(alp)/alpgen.inc 
	$(FFF) -o $(execfile)90V $(usrfile).o $(prcfile).o \
	$(alp)/alpgen.o $(alp)/alputi.o $(alp)/alppdf.o \
	$(prclib)/XXX.o $(alp)/Aint90.o $(VF90)/libvast90.a

# DIRECTORY CLEANUP UTILITIES:
#
# remove object files only
cleanobj:
	-rm $(PARTON:.f=.o) $(PARTON90:.f90=.o) $(prcusr)/../*/*.o90*

# remove object files, etc
cleanall:
	-rm $(OBJ) $(prcusr)/fort.* $(prcusr)/*.top $(prcusr)/*.par \
	$(prcusr)/*.wgt $(prcusr)/*.unw $(prcusr)/*.mon \
	$(prcusr)/*.stat $(prcusr)/../*/*.o90*

