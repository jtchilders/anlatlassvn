#-*- Makefile -*-
## Default intel compiler suite make COMPILER=gnu for gfortran/g++
ifndef COMPILER
COMPILER=gnu
endif
## Choose PDF: native,lhapdf
## LHAPDF package has to be installed separately
PDF=lhapdf
#Choose Analysis: none,  antikt
## default analysis may require FASTJET package, that has to be installed separately (see below)
ANALYSIS=antikt
#
#
ifeq ("$(COMPILER)","gnu")
FC = gfortran
F77 = $(FC) -fno-automatic -ffixed-line-length-132
CXX = g++
CC = gcc	
## -fbounds-check sometimes causes a weird error due to non-lazy evaluation
## of boolean in gfortran.
FFLAGS= -fbounds-check -Wall #-Wimplicit-procedure 
## For floating point exception trapping  uncomment the following 
FPE=-ffpe-trap=invalid,zero,overflow
## gfortran 4.4.1 optimized with -O3 yields erroneous results
## Use -O2 to be on the safe side
OPT=-O2
## For debugging uncomment the following
#DEBUG= -ggdb #-pg
##
endif

ifeq ("$(COMPILER)","intel")
FC = ifort
F77 = $(FC) -save -assume bscc -extend_source
CXX = icpc
CC = icc
LIBS = -limf
FFLAGS =  -check #-warn all -module obj-$(COMPILER)/mod 
## For floating point exception trapping  uncomment the following 
#FPE = -fpe0
OPT = -O3 #-fast
## For debugging uncomment the following
#DEBUG= -debug -g
#PP=-fpp #-DHELAS_CHECK
endif




PWD=$(shell pwd)
WDNAME=$(shell basename $(PWD))
DIRVIRTUALS=$(PWD)/Virtuals
VPATH= ./:../:obj-$(COMPILER)/:madgraph/:madgraph/born/:madgraph/real/:madgraph/born_decay/:madgraph/real_decay/


INCLUDE0=$(PWD)
INCLUDE1=$(PWD)/include
INCLUDE2=$(shell dirname $(PWD))/include
INCLUDE3=$(shell dirname $(PWD))

FF=$(F77) $(FPE) $(OPT) $(PP) $(DEBUG) -I$(INCLUDE0) -I$(INCLUDE1) -I$(INCLUDE2)  -I$(INCLUDE3)


INCLUDE =$(wildcard ../include/*.h *.h include/*.h madgraph/*.inc)


ifeq ("$(PDF)","lhapdf")
LHAPDF_CONFIG=lhapdf-config
PDFPACK=lhapdfif.o
LIBSLHAPDF= -Wl,-rpath,$(shell $(LHAPDF_CONFIG) --libdir)  -L$(shell $(LHAPDF_CONFIG) --libdir)  -lLHAPDF 
## If LHAPDF has been compiled with gfortran and you want to link it statically, you have to include
## libgfortran as well. The same holds for libstdc++. 
## One possible solution is to use fastjet, since $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
## does perform this inclusion. The path has to be set by the user. 
# LIBGFORTRANPATH=/usr/lib/gcc/x86_64-redhat-linux/4.1.2
# LIBSTDCPP=/usr/lib/gcc/x86_64-redhat-linux/4.1.1
# LIBSLHAPDF+=  -L$(LIBGFORTRANPATH)  -lgfortranbegin -lgfortran -L$(LIBSTDCPP) -lstdc++
LIBS+=$(LIBSLHAPDF)
else
PDFPACK=mlmpdfif.o hvqpdfpho.o
endif


ifeq ("$(ANALYSIS)","antikt")
##To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins )  -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_analysis-antikt.o 
PWHGANAL+= fastjetfortran.o
else
PWHGANAL=pwhg_analysis-dummy.o
endif

LIBS+=$(LIBSFASTJET)



LIBDHELAS=madgraph/DHELAS/libdhelas3-$(FC).a
LIBS += $(LIBDHELAS)

# Libraries added  before here can be linked statically .
# From now on only dinamically 
LIBVIRTUALS=libvirtual.so
LIBVIRTUALS_PATH=$(DIRVIRTUALS)/obj-$(COMPILER)
LIBS += -Wl,-eh-frame-hdr -shared-libgcc -Wl,-Bdynamic -Wl,-rpath,$(LIBVIRTUALS_PATH) -L$(LIBVIRTUALS_PATH) -lvirtual


%.o: %.f $(INCLUDE)
	$(FF) $(FFLAGS) -c -o obj-$(COMPILER)/$@ $<

%.o: %.c
	$(CC) $(DEBUG) -c -o obj-$(COMPILER)/$@ $^ 

%.o: %.cc
	$(CXX) $(DEBUG) -c -o obj-$(COMPILER)/$@ $^ $(FJCXXFLAGS)

%.o: %.cpp FiniteVirtuals.h 
	$(CXX) $(DEBUG) -c -o obj-$(COMPILER)/$@ $< $(FJCXXFLAGS) $(IFLAGS)


MADGRAPH=functions.o setparam.o 

# Matrix elements needed for NLO calculation
MADGRAPH+=$(shell echo madgraph/born/*.f | sed -e 's/.f/.o/g;s/madgraph\/born\///g') 
MADGRAPH+=$(shell echo madgraph/real/*.f | sed -e 's/.f/.o/g;s/madgraph\/real\///g') 

# Matrix elements needed for decay
MADGRAPH+=$(shell echo madgraph/born_decay/*.f | sed -e 's/.f/.o/g;s/madgraph\/born\_decay\///g') 
MADGRAPH+=$(shell echo madgraph/real_decay/*.f | sed -e 's/.f/.o/g;s/madgraph\/real\_decay\///g') 


USER=init_couplings.o init_processes.o Born_phsp.o Born.o virtual.o \
        real.o $(PWHGANAL) $(MADGRAPH)

PWHG= pwhg_main.o pwhg_init.o bbinit.o btilde.o lhefwrite.o		\
	ttjdecay.o LesHouches.o LesHouchesreg.o gen_Born_phsp.o		\
	find_regions.o test_Sudakov.o pt2maxreg.o sigborn.o		\
	gen_real_phsp.o maxrat.o gen_index.o gen_radiation.o		\
	Bornzerodamp.o sigremnants.o random.o boostrot.o		\
	lhefroutines.o bra_ket_subroutines.o cernroutines.o		\
	init_phys.o powheginput.o pdfcalls.o sigreal.o sigcollremn.o	\
	pwhg_bookhist.o pwhg_bookhist_add.o pwhg_analysis_driver.o	\
	reshufflemoms.o checkmomzero.o setstrongcoupl.o integrator.o	\
	newunit.o mwarn.o sigsoftvirt.o sigcollsoft.o sigvirtual.o	\
	$(PDFPACK) $(USER) $(FPEOBJ)

# main target, to generate the hardest emission and write events on a LHEF file
pwhg_main: $(PWHG) $(LIBDHELAS) $(LIBVIRTUALS) 
	$(FF) $(patsubst %.o,obj-$(COMPILER)/%.o,$(PWHG)) $(LIBS) -o $@-$(COMPILER)

LHEF=lhef_analysis.o boostrot.o random.o cernroutines.o opencount.o	\
	powheginput.o pwhg_bookhist.o pwhg_bookhist_add.o lhefread.o	\
	newunit.o pwhg_analysis_driver.o printleshouches.o		\
	$(PWHGANAL) $(FPEOBJ)

# target to analyze LHEF output
lhef_analysis:$(LHEF)
	$(FF) $(patsubst %.o,obj-$(COMPILER)/%.o,$(LHEF))  $(LIBSFASTJET)  -o $@-$(COMPILER)


# target to read event file, shower events with HERWIG + analysis
HERWIG=main-HERWIG.o herwig.o boostrot.o powheginput.o newunit.o			\
	pwhg_bookhist.o pwhg_bookhist_add.o $(PWHGANAL) lhefread.o		\
	pdfdummies.o cernroutines.o opencount.o printleshouches.o

HERWIG_DRIVER_NOUE= setup-HERWIG-lhef.o

main-HERWIG-lhef: $(HERWIG) $(HERWIG_DRIVER_NOUE) 
	$(FF) $(patsubst %,obj-$(COMPILER)/%, $(HERWIG) $(HERWIG_DRIVER_NOUE))  $(LIBSFASTJET)  -o $@-$(COMPILER)

HERWIG_DRIVER_UE= setup-HERWIG-UE-lhef.o

main-HERWIG-UE-lhef: $(HERWIG) $(HERWIG_DRIVER_UE) 
	$(FF) $(patsubst %,obj-$(COMPILER)/%, $(HERWIG) $(HERWIG_DRIVER_UE))  $(LIBSFASTJET)  -o $@-$(COMPILER)


# target to read event file, shower events with PYTHIA + analysis
PYTHIA=main-PYTHIA.o  pythia.o boostrot.o powheginput.o	\
	pwhg_bookhist.o pwhg_bookhist_add.o $(PWHGANAL) lhefread.o pdfdummies.o			\
        cernroutines.o opencount.o printleshouches.o newunit.o

PYTHIA_DRIVER_NOMPI=setup-PYTHIA-lhef.o

main-PYTHIA-lhef: $(PYTHIA) $(PYTHIA_DRIVER_NOMPI)
	$(FF) $(patsubst %,obj-$(COMPILER)/%, $(PYTHIA) $(PYTHIA_DRIVER_NOMPI))  $(LIBSFASTJET)  -o $@-$(COMPILER)

PYTHIA_DRIVER_MPI=setup-PYTHIA-MPI-lhef.o

main-PYTHIA-MPI-lhef: $(PYTHIA) $(PYTHIA_DRIVER_MPI)
	$(FF) $(patsubst %,obj-$(COMPILER)/%, $(PYTHIA) $(PYTHIA_DRIVER_MPI))  $(LIBSFASTJET)  -o $@-$(COMPILER)


# 
all: pwhg_main lhef_analysis main-HERWIG-lhef main-PYTHIA-lhef main-HERWIG-UE-lhef main-PYTHIA-MPI-lhef

# target to cleanup
.PHONY: dist-clean clean
clean:
	rm -f obj-$(COMPILER)/*.o madgraph/DHELAS/obj-$(FC)/*.o madgraph/DHELAS/*$(FC).a pwhg_main-$(COMPILER) lhef_analysis-$(COMPILER) main-HERWIG-lhef-$(COMPILER)	\
	main-PYTHIA-lhef-$(COMPILER)

dist-clean: clean
	rm -f obj*/*.o obj*/*/* madgraph/DHELAS/*/*.o madgraph/DHELAS/*.a \
	pwhg_main-* lhef_analysis-* main-HERWIG-lhef-* main-PYTHIA-lhef-*
	@echo looking into subdir $(DIRVIRTUALS)
	$(MAKE) dist-clean --directory=$(DIRVIRTUALS) COMPILER=$(COMPILER)

###### to compile dependencies
$(LIBDHELAS): force_look
	@echo looking into subdir madgraph/DHELAS
	$(MAKE) FC=$(FC) FF='$(FF)' --directory=madgraph/DHELAS 

$(LIBVIRTUALS): force_look
	@echo looking into subdir $(DIRVIRTUALS)
	$(MAKE) --directory=$(DIRVIRTUALS) COMPILER=$(COMPILER) DEBUG='$(DEBUG)'
force_look :
	true
######

