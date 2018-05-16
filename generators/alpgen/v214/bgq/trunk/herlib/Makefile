# DEFINE FORTRAN COMPILATION DEFAULTS
include ../compile.mk
#
# DEFINE DIRECTORY TREE STRUCTURE:
her= .
alp= ../alplib

# DEFINE FILE GROUPS:
#

# Files used for the shower evolution
SHOWER= $(her)/*.f $(alp)/alpsho.f

# Include files
INC= $(her)/alpsho.inc \
	$(her)/HERWIG65.INC $(her)/herwig6510.inc 
#$(her)/HERWIG64.INC $(her)/herwig6400.inc 

# object files
OBJ=$(SHOWER:.f=.o) 

# include files' dependencies
$(SHOWER): $(her)/alpsho.inc \
	$(her)/HERWIG65.INC $(her)/herwig6510.inc 
#$(her)/HERWIG64.INC $(her)/herwig6400.inc 

# compilation
%.o: %.f
	$(FFF) -c -o $*.o $*.f 


# LINKING GENERIC EXECUTABLE FOR HERWIG EVENT-EVOLUTION
hwuser: $(her)/hwuser.o $(her)/atoher.o \
        $(her)/herwig6510.o $(her)/pdfdummy.o $(alp)/alpsho.o
	$(FFF) -o $(her)/hwuser \
        $(her)/hwuser.o $(her)/atoher.o \
        $(her)/herwig6510.o $(her)/pdfdummy.o $(alp)/alpsho.o 

# DIRECTORY CLEANUP UTILITIES:
#
# remove object files only
cleanobj:
	-rm $(OBJ) $(alp)/*.o $(alp)/alpsho.o

# remove object files, etc
cleanall:
	-rm $(OBJ) $(her)/fort.* $(her)/*.top $(her)/*.par \
	$(her)/*.wgt $(her)/*.unw $(her)/*.mon \
	$(her)/*.stat $(alp)/alpsho.o





