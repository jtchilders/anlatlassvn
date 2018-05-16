#-*- Makefile -*-

# various files included here are based on POWHEG-BOX/Z

# choose compiler: gfortran,ifort
COMPILER=gfortran
#
# choose PDF: native,lhapdf
# LHAPDF package has to be installed separately
PDF = lhapdf
#
# choose analysis: none, default, all
# default+all: requires FASTJET package, has to be installed separately (see below)
# all: requires MT2 code by Cheng and Han, has to be downloaded to MT2DIR (see below)
ANALYSIS=default
#
# for static linking uncomment the following
#STATIC=-static
#
# paths to libraries and header directories
PWD=$(shell pwd)
# please modify if you prefer to use your own version of SLHALib:
LIBSLHA = $(PWD)/SLHALib-2.2/build/libSLHA.a
INCSLHA = $(PWD)/SLHALib-2.2/src
# comment out to avoid autobuild:
BUILDSLHA = $(PWD)/SLHALib-2.2
# please modify if you prefer to use your own version of QCDLoop:
LIBQCDLOOP = $(PWD)/QCDLoop-1.9/ql/libqcdloop.a
LIBFF =      $(PWD)/QCDLoop-1.9/ff/libff.a
# comment out to avoid autobuild:
BUILDQCDLOOP = $(PWD)/QCDLoop-1.9
# directory for mt2_bisect.[h|cpp] files, only needed for ANALYSIS=all
# (see file mt2_bisect_wrap.cpp for details on how to obtain these files)
MT2DIR=$(PWD)/mt2_bisect-1.01a


ifeq ("$(COMPILER)","gfortran")	
F77 = gfortran -fno-automatic -ffixed-line-length-none	
CXX = g++
# -fbounds-check sometimes causes a weird error due to non-lazy evaluation
# of boolean in gfortran.
#FFLAGS= -Wall -Wno-implicit-interface -Wtabs -fbounds-check
# for floating point exception trapping  uncomment the following 
#FPE=-ffpe-trap=invalid,zero,overflow,underflow 
# gfortran 4.4.1 with -O3 yields erroneous results, use -O2 to be on safe side
OPT=-O2
# for debugging uncomment the following
#DEBUG= -ggdb -pg
endif
ifeq ("$(COMPILER)","ifort")
F77 = ifort -save -extend_source
#CXX = icpc
CXX = g++
LIBS = -limf
#FFLAGS =  -check
# for floating point exception trapping  uncomment the following 
#FPE = -fpe0
OPT = -O3 #-fast
# for debugging uncomment the following
#DEBUG= -debug -g
endif


ifeq ("$(PDF)","lhapdf")
LHAPDF_CONFIG=lhapdf-config
PDFPACK=lhapdfif.o
LIBSLHAPDF= -Wl,-rpath,$(shell $(LHAPDF_CONFIG) --libdir)  -L$(shell $(LHAPDF_CONFIG) --libdir) -lLHAPDF
ifeq  ("$(STATIC)","-static") 
# if LHAPDF has been compiled with gfortran and you want to link it statically, you have to include
# libgfortran as well. The same holds for libstdc++. 
# One possible solution is to use fastjet, since $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
# does perform this inclusion. The path has to be set by the user. 
 LIBGFORTRANPATH=/usr/lib/gcc/x86_64-redhat-linux/4.1.2
 LIBSTDCPP=/lib64
 LIBSLHAPDF+= -L$(LIBGFORTRANPATH) -lgfortranbegin -lgfortran -L$(LIBSTDCPP) -lstdc++
endif
LIBS+=$(LIBSLHAPDF)
else
PDFPACK=mlmpdfif.o hvqpdfpho.o
endif


ifeq ("$(ANALYSIS)","default")
USE_FASTJET=yes
PWHGANAL=pwhg_analysis.o
endif
ifeq ("$(ANALYSIS)","all")
USE_FASTJET=yes
USE_MT2=yes
PWHGANAL=pwhg_analysis.o
endif
ifeq ("$(ANALYSIS)","none")
PWHGANAL=pwhg_analysis-dummy.o
endif

ifdef USE_FASTJET
FASTJET_CONFIG=$(shell which fastjet-config)
ifeq ($(FASTJET_CONFIG),)
$(error cannot find fastjet-config (did you install fastjet ?))
endif
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL+= fastjetfortran.o fastjetsisconewrap.o fastjetktwrap.o
endif

ifdef USE_MT2
PWHGANAL+= mt2_bisect_wrap.o mt2_bisect.o
FFLAGS+= -DHAVE_MT2CODE
INCLUDES+= -I$(MT2DIR)
CXXFLAGS+= -I$(MT2DIR)
endif


VPATH= ./:obj/:
HEADERS =$(wildcard ../include/*.h *.h include/*.h)
INCLUDES+=-I$(PWD) -I$(shell dirname $(PWD))/include -I$(INCSLHA)
FF=$(F77) $(FFLAGS) $(FPE) $(OPT) $(DEBUG) $(INCLUDES)



# build these targets by default
all: pwhg_main lhef_analysis main-PYTHIA-lhef


# target to generate LHEF output
PWHG=pwhg_main.o pwhg_init.o bbinit.o btilde.o lhefwrite.o		\
	LesHouches.o LesHouchesreg.o gen_Born_phsp.o find_regions.o	\
	test_Sudakov.o pt2maxreg.o sigborn.o gen_real_phsp.o maxrat.o	\
	gen_index.o gen_radiation.o Bornzerodamp.o sigremnants.o	\
	random.o boostrot.o bra_ket_subroutines.o cernroutines.o	\
	init_phys.o powheginput.o pdfcalls.o sigreal.o sigcollremn.o	\
	pwhg_bookhist.o pwhg_analysis_driver.o checkmomzero.o		\
	setstrongcoupl.o integrator.o newunit.o mwarn.o sigsoftvirt.o	\
	sigcollsoft.o sigvirtual.o 	\
    init_couplings.o init_processes.o \
    Born_phsp.o Born.o Born_lh.o virtual.o real.o \
    funcprocess.o funcbasic.o funcloop.o \
    $(PWHGANAL) $(PDFPACK) $(FPEOBJ)
PWHGEXTRA=$(LIBQCDLOOP) $(LIBFF) $(LIBSLHA)
pwhg_main:$(PWHG) $(PWHGEXTRA) Makefile
	$(FF) $(patsubst %,obj/%,$(PWHG)) $(PWHGEXTRA) $(LIBS) $(LIBSFASTJET) $(STATIC) -o $@

# target to analyze LHEF output
LHEF=lhef_analysis.o boostrot.o random.o cernroutines.o	\
     opencount.o powheginput.o pwhg_bookhist.o \
     lhefread.o newunit.o pwhg_analysis_driver.o \
     funcprocess.o init_couplings.o \
     $(PWHGANAL) $(FPEOBJ)
LHEFEXTRA=$(LIBSLHA)
lhef_analysis:$(LHEF) $(LHEFEXTRA)
	$(FF) $(patsubst %,obj/%,$(LHEF)) $(LHEFEXTRA) $(LIBS) $(LIBSFASTJET) $(STATIC)  -o $@ 

# target to read event file, shower events with PYTHIA + analysis
# NOTE: pythia-6.4.25-decay.f is based on POWHEG-BOX/pythia-6.4.25.f,
# to which we applied a simple patch to enable slepton decays
# (for clarity a separate patch file is included which is not needed otherwise)
PYTHIA=main-PYTHIA.o setup-PYTHIA-lhef.o pythia-6.4.25-decay.o \
	boostrot.o powheginput.o init_couplings.o funcprocess.o \
	pwhg_bookhist.o lhefread.o newunit.o pdfdummies.o \
	pwhg_analysis_driver.o random.o cernroutines.o opencount.o	\
	$(PWHGANAL) $(FPEOBJ)
PYTHIAEXTRA=$(LIBSLHA)
main-PYTHIA-lhef: $(PYTHIA) $(PYTHIAEXTRA)
	$(FF) $(patsubst %,obj/%,$(PYTHIA)) $(PYTHIAEXTRA) $(LIBSFASTJET)  $(STATIC) -o $@


# autobuild targets for scalar loop integral libraries FF and QCDLoop
ifdef BUILDQCDLOOP
$(LIBFF):
	cd $(BUILDQCDLOOP)/ff && $(MAKE) FC=$(COMPILER)
$(LIBQCDLOOP): $(LIBFF)
	cd $(BUILDQCDLOOP)/ql && $(MAKE) FC=$(COMPILER)
endif

# autobuild targets for SLHALib needed for SUSY parameter I/O
ifdef BUILDSLHA
$(LIBSLHA): 
	cd $(BUILDSLHA) && ./configure FC=$(COMPILER) && $(MAKE)
endif

# MT2 code
mt2_bisect.o: $(MT2DIR)/mt2_bisect.cpp
	$(CXX) $(CXXFLAGS) -c -o obj/$@ $^


# project specific Fortran files
%.o: %.f $(HEADERS)
	$(FF) -c -o obj/$@ $<

# project specific Fortran files which need the preprocessor (.F instead of .f)
%.o: %.F $(HEADERS)
	$(FF) -c -o obj/$@ $<

# project specific C++ files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o obj/$@ $^

# shared Fortran files
%.o: ../%.f $(HEADERS)
	$(FF) -c -o obj/$@ $<

# shared C files
%.o: ../%.c
	$(CC) $(DEBUG) -c -o obj/$@ $^ 

# shared C++ files
%.o: ../%.cc
	$(CXX) $(DEBUG) -c -o obj/$@ $^ $(FJCXXFLAGS)


# target to cleanup everything except for external libraries
.PHONY: mostlyclean
mostlyclean:
	rm -f obj/*.o \
	      pwhg_main lhef_analysis main-PYTHIA-lhef

# target to cleanup everything
.PHONY: clean
clean: mostlyclean
ifdef BUILDQCDLOOP
	cd $(BUILDQCDLOOP) && make clean && cd ql && make clean && cd ../ff/ && make clean
endif
ifdef BUILDSLHA
	cd $(BUILDSLHA) && (test -f makefile && make clean && rm makefile || true)
endif
	