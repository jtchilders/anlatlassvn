prc=hjet
include ../compile.mk
include ../alplib/alpgen.mk
$(prclib)/VX77.o : $(prclib)/77/VX77.f
	$(FFF) -c  $(prclib)/77/VX77.f -o $(prclib)/VX77.o
gen: $(usrfile).o $(prcfile).o $(prclib)/$(prc).inc\
	$(alp)/alpgen.o $(alp)/alputi.o $(alp)/alppdf.o \
	$(alp)/Aint90.o $(prclib)/VX77.o $(alp)/alpgen.inc 
	$(FFF) -o $(execfile) $(usrfile).o $(prcfile).o \
	$(alp)/alpgen.o $(alp)/alputi.o $(alp)/alppdf.o \
	$(prclib)/VX77.o $(alp)/Aint90.o 
