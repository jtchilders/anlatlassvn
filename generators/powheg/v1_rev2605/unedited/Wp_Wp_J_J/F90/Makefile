
ifeq ($(FC),ifort)
MAKEFILE=Makefile.ifort
endif
ifeq ($(FC),gfortran)
MAKEFILE=Makefile.gfortran
endif

all:
	make -f $(MAKEFILE)



clean:
	make -f $(MAKEFILE) realclean
