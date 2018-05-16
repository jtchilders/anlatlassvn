#-*- Makefile -*-
## Choose compiler: gfortran,g77,ifort
COMPILER=gfortran
FC=$(COMPILER)

## Choose PDF: native,lhapdf
## LHAPDF package has to be installed separately
PDF=lhapdf

# Choose Analysis: none, default
# default analysis may require FASTJET package, that has to be installed 
# separately (see below)
ANALYSIS=default
## For static linking uncomment the following
#STATIC= -static
#
ifeq ("$(COMPILER)","gfortran")	
#FFLAGS= -Wall 
## For floating point exception trapping  uncomment the following 
# FPEOBJ=trapfpe.o
F77= gfortran -fno-automatic -ffixed-line-length-none 
# 	-Wall -fbounds-check 
#	-Wsurprising -Wno-unused -finit-real=NAN -finit-integer=-1000000 -Wtabs 
#	-ffpe-trap=invalid,zero,overflow,underflow -fbacktrace 
OPT=-O2 
## For debugging uncomment the following
#DEBUG= -ggdb -pg
endif

ifeq ("$(FC)","g77")
F77= g77 -fno-automatic 
#FFLAGS= -Wall -ffortran-bounds-check
## For floating point exception trapping  uncomment the following 
#FPEOBJ=trapfpe.o
OPT=-O3
## For debugging uncomment the following
#DEBUG= -ggdb -pg
endif


ifeq ("$(FC)","ifort")
F77 = ifort -save -132 
#  -extend_source
#CXX = icpc
LIBS = -limf
#FFLAGS =  -check
## For floating point exception trapping  uncomment the following 
#FPE = -fpe0
OPT = -O3 #-fast
## For debugging uncomment the following
#DEBUG= -debug -g
endif



PWD=$(shell pwd)
WDNAME=$(shell basename $(PWD))
VPATH= ./:../:obj/

INCLUDE0=$(PWD)
INCLUDE1=$(PWD)/include
INCLUDE2=$(shell dirname $(PWD))/include
INCLUDE3=$(PWD)/vbfnlo-files
FF=$(F77) $(FFLAGS) $(FPE) $(OPT) $(DEBUG) -I$(INCLUDE0) -I$(INCLUDE1) -I$(INCLUDE2) -I$(INCLUDE3)


INCLUDE =$(wildcard ../include/*.h *.h include/*.h vbfnlo-files/*.inc vbfnlo-files/*.h)

ifeq ("$(PDF)","lhapdf")
# LHAPDF may be compiled with ifort or gfortran. The user might need to 
# set up the correct path.
LHAPDF_CONFIG=lhapdf-config
PDFPACK=lhapdfif.o
LIBSLHAPDF= -Wl,-rpath,$(shell $(LHAPDF_CONFIG) --libdir)  -L$(shell $(LHAPDF_CONFIG) --libdir) -lLHAPDF
ifeq  ("$(STATIC)","-static") 
## If LHAPDF has been compiled with gfortran and you want to link it statically, you have to include
## libgfortran as well. The same holds for libstdc++. 
## One possible solution is to use fastjet, since $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
## does perform this inclusion. The path has to be set by the user. 
 LIBGFORTRANPATH=/usr/lib/gcc/x86_64-redhat-linux/4.1.2
 LIBSTDCPP=/lib64
 LIBSLHAPDF+=  -L$(LIBGFORTRANPATH)  -lgfortranbegin -lgfortran -L$(LIBSTDCPP) -lstdc++
endif
LIBS+=$(LIBSLHAPDF)
else
PDFPACK=mlmpdfif.o hvqpdfpho.o
endif

ifeq ("$(ANALYSIS)","default")
##To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config) 
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_analysis.o fastjetfortran.o
# Also add required Fastjet drivers to PWHGANAL (examples are reported)
PWHGANAL+= fastjetsisconewrap.o fastjetktwrap.o #fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o
else
PWHGANAL=pwhg_analysis-dummy.o
endif

LIBS+=$(LHAPDFLIB) -Lvbfnlo-files -lvbfnlo 

%.o: %.f $(INCLUDE)
	$(FF) -c -o obj/$@ $<

%.o: %.c
	$(CC) $(DEBUG) -c -o obj/$@ $^ 

%.o: %.cc
	$(CC) $(DEBUG) -c -o obj/$@ $^ $(FJCXXFLAGS)


USER=init_couplings.o init_processes.o  Born.o virtual.o  calc_invariants.o\
	 breit.o Born_phsp2.o Born_phspfile.o boost.o \
	convert_coup.o real.o $(PWHGANAL)

PWHG=pwhg_main.o pwhg_init.o bbinit.o btilde.o lhefwrite.o            \
	LesHouches.o LesHouchesreg.o gen_Born_phsp.o find_regions.o   \
        test_Sudakov.o pt2maxreg.o sigborn.o gen_real_phsp.o maxrat.o \
	gen_index.o gen_radiation.o Bornzerodamp.o sigremnants.o      \
	random.o boostrot.o bra_ket_subroutines.o cernroutines.o      \
	init_phys.o powheginput.o pdfcalls.o sigreal.o sigcollremn.o  \
        pwhg_bookhist.o pwhg_analysis_driver.o checkmomzero.o         \
	reshufflemoms.o setstrongcoupl.o integrator.o newunit.o       \
	mwarn.o sigsoftvirt.o sigcollsoft.o sigvirtual.o              \
	$(PDFPACK) $(USER) $(SYSOBJ)

# target to generate LHEF output
pwhg_main:$(PWHG) vbfnlo-files/libvbfnlo.a
	(cd vbfnlo-files; make;)
	$(FF) $(patsubst %,obj/%, $(PWHG)) $(LIBS) $(LIBSFASTJET) -o $@  $(STATIC)


vbfnlo-files/libvbfnlo.a: vbfnlo-files/*f  vbfnlo-files/*/*f
	cd vbfnlo-files ; make

LHEFANAL=lhef_analysis.o boostrot.o random.o cernroutines.o         \
	 opencount.o powheginput.o pwhg_bookhist.o  \
	 lhefread.o newunit.o pwhg_analysis_driver.o \
	 $(PWHGANAL) $(SYSOBJ)


# target to analyze LHEF output
lhef_analysis:$(LHEFANAL) vbfnlo-files/libvbfnlo.a
	(cd vbfnlo-files; make;)
	$(FF) $(patsubst %,obj/%,$(LHEFANAL)) -o $@ $(LIBS) $(LIBSFASTJET) $(STATIC) 


# target to read event file, shower events with HERWIG + analysis

HERWIG=main-HERWIG.o setup-HERWIG-lhef.o herwig6510.o \
	boostrot.o powheginput.o \
        pwhg_bookhist.o lhefread.o pdfdummies.o newunit.o \
	pwhg_analysis_driver.o \
	random.o cernroutines.o opencount.o  $(PWHGANAL) $(SYSOBJ)


main-HERWIG-lhef: $(HERWIG)vbfnlo-files/libvbfnlo.a
	(cd vbfnlo-files; make;)
	$(FF) $(patsubst %,obj/%,$(HERWIG)) -o $@ $(LIBSFASTJET) $(LIBS)


PYTHIA=main-PYTHIA.o setup-PYTHIA-lhef.o pythia-6.4.25.o \
	boostrot.o powheginput.o \
        pwhg_bookhist.o lhefread.o newunit.o \
	pwhg_analysis_driver.o pdfdummies.o \
	random.o cernroutines.o opencount.o  $(PWHGANAL) $(SYSOBJ)

main-PYTHIA-lhef: $(PYTHIA) vbfnlo-files/libvbfnlo.a
	(cd vbfnlo-files; make;)
	$(FF) $(patsubst %,obj/%,$(PYTHIA)) $(LIBS) $(LIBSFASTJET)  -o $@ $(STATIC)

clean:
	rm -f *.o obj_ew/*.o obj/*.o \
        vbfnlo-files/*.o  vbfnlo-files/libvbfnlo.a  vbfnlo-files/*/*.o \
        pwhg_main lhef_analysis main-HERWIG-lhef \
	main-PYTHIA-lhef
