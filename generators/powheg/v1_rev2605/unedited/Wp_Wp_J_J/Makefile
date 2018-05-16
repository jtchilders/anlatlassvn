#-*- Makefile -*-
## Choose compiler: gfortran,g77,ifort
COMPILER=gfortran
FC=$(COMPILER)
## Choose PDF: native,lhapdf
## LHAPDF package has to be installed separately
PDF=native
#Choose Analysis: none, default
## default analysis may require FASTJET package, that has to be installed separately (see below)
ANALYSIS=none
## For static linking uncomment the following
#STATIC= -static
#
ifeq ("$(FC)","gfortran")	
F77= gfortran -fno-automatic 	
## -fbounds-check sometimes causes a weird error due to non-lazy evaluation
## of boolean in gfortran.
#FFLAGS= -Wall -Wimplicit-interface -fbounds-check
## For floating point exception trapping  uncomment the following 
#FPE=-ffpe-trap=invalid,zero,overflow,underflow 
## gfortran 4.4.1 optimized with -O3 yields erroneous results
## Use -O2 to be on the safe side
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
F77 = ifort -save
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

ifdef DEBUG
#FPE=-ffpe-trap=invalid,zero,overflow
#,underflow
OPT=-O0
endif



PWD=$(shell pwd)
WDNAME=$(shell basename $(PWD))
VPATH= ./:../:obj-$(FC)/

INCLUDE0=$(PWD)
INCLUDE1=$(PWD)/include
INCLUDE2=$(shell dirname $(PWD))/include
FF=$(F77) $(FFLAGS) $(FPE) $(OPT) $(DEBUG) -I$(INCLUDE0) -I$(INCLUDE1) -I$(INCLUDE2) -I F90 -I F90/obj-$(FC)


INCLUDE =$(wildcard ../include/*.h *.h include/*.h)

ifeq ("$(PDF)","lhapdf")
# LHAPDF may be compiled with ifort or gfortran. The user should
# set up the correct path.
LHAPDF_CONFIG=~/Pheno/PDFpacks/lhapdf-5.8.4-$(FC)/bin/lhapdf-config
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
## Also add required Fastjet drivers to PWHGANAL (examples are reported)
#PWHGANAL+= fastjetsisconewrap.o fastjetktwrap.o fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o
else
PWHGANAL=pwhg_analysis-dummy.o
endif


QCDLOOP=QCDLoop-1.9-$(FC)

LIBS+= -LF90/obj-$(FC) -lWpWp2j -L$(QCDLOOP)/ql -L$(QCDLOOP)/ff -lqcdloop -lff

%.o: %.f $(INCLUDE)
	$(FF) -c $(DEBUG) -o obj-$(FC)/$@ $<

%.o: %.c
	$(CC) $(DEBUG) -c -o obj-$(FC)/$@ $^ 

%.o: %.cc
	$(CXX) $(DEBUG) -c -o obj-$(FC)/$@ $^ $(FJCXXFLAGS)

USER=init_couplings.o init_processes.o Born_phsp.o Born.o virtual.o	\
     real.o $(PWHGANAL)	boost.o phi1_2.o phi3m0.o breit.o

PWHG=pwhg_main.o pwhg_init.o bbinit.o btilde.o lhefwrite.o		\
	LesHouches.o LesHouchesreg.o gen_Born_phsp.o find_regions.o	\
	test_Sudakov.o pt2maxreg.o sigborn.o gen_real_phsp.o maxrat.o	\
	gen_index.o gen_radiation.o Bornzerodamp.o sigremnants.o	\
	random.o boostrot.o bra_ket_subroutines.o cernroutines.o	\
	init_phys.o powheginput.o pdfcalls.o sigreal.o sigcollremn.o	\
	pwhg_bookhist.o pwhg_analysis_driver.o checkmomzero.o		\
	setstrongcoupl.o integrator.o newunit.o mwarn.o sigsoftvirt.o	\
	sigcollsoft.o sigvirtual.o  mint_upb.o  $(PDFPACK) $(USER) $(FPEOBJ)


# target to generate LHEF output
pwhg_main: F90/obj-$(FC)/libWpWp2j.a $(QCDLOOP)/ql/libqcdloop.a $(QCDLOOP)/ff/libff.a  $(PWHG)
	$(FF) $(patsubst %,obj-$(FC)/%,$(PWHG)) $(LIBS) $(LIBSFASTJET) -o $@ $(IFORTLIBS) $(STATIC)


F90/obj-$(FC)/libWpWp2j.a:F90/*.f* F90/$(FC)/*.f*
	(cd F90 ; make FC=$(FC))

$(QCDLOOP)/ql/libqcdloop.a:
	(cd $(QCDLOOP)/ql; make FC=$(FC))

$(QCDLOOP)/ff/libff.a:
	(cd $(QCDLOOP)/ff; make FC=$(FC))


LHEF=lhef_analysis.o boostrot.o random.o cernroutines.o		\
     opencount.o powheginput.o pwhg_bookhist.o $(PWHGANAL)	\
     lhefread.o newunit.o pwhg_analysis_driver.o $(FPEOBJ)

# target to analyze LHEF output
lhef_analysis:$(LHEF)
	$(FF) $(patsubst %,obj/%,$(LHEF)) $(LIBS) $(LIBSFASTJET) $(STATIC)  -o $@ 



# target to read event file, shower events with HERWIG + analysis
HERWIG=main-HERWIG.o setup-HERWIG-lhef.o herwig.o boostrot.o	\
	powheginput.o pwhg_bookhist.o $(PWHGANAL) lhefread.o	\
	pdfdummies.o opencount.o $(FPEOBJ) 

main-HERWIG-lhef: $(HERWIG)
	$(FF) $(patsubst %,obj/%,$(HERWIG))  $(LIBSFASTJET)  $(STATIC) -o $@

# target to read event file, shower events with PYTHIA + analysis
PYTHIA=main-PYTHIA.o setup-PYTHIA-lhef.o pythia.o boostrot.o powheginput.o		\
	pwhg_bookhist.o $(PWHGANAL) lhefread.o newunit.o pdfdummies.o	\
	pwhg_analysis_driver.o random.o cernroutines.o opencount.o	\
	$(FPEOBJ)

main-PYTHIA-lhef: $(PYTHIA)
	$(FF) $(patsubst %,obj/%,$(PYTHIA)) $(LIBSFASTJET)  $(STATIC) -o $@

# target to cleanup
.PHONY: clean

clean:
	cd QCDLoop-1.9-$(FC); make clean ; cd ql ; make clean ; cd ../ff/ ; make clean
	cd F90 ; make FC=$(FC) clean
	rm -f obj-$(FC)/*.o pwhg_main lhef_analysis main-HERWIG-lhef main-PYTHIA-lhef
