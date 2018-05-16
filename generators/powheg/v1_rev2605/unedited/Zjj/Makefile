#-*- Makefile -*-

VIRTUAL=none
# none,BH,MCFM

GRID=false
BH=false

BHPATH=/usersVol2/re/C_Cpp/Blackhat
QDPATH=/usersVol2/re/.Packages_ER/qd_v_2_3_13

## Choose compiler: gfortran,g77,ifort
COMPILER=gfortran
## Choose PDF: native,lhapdf
## LHAPDF package has to be installed separately
PDF=native
#native,lhapdf
#Choose Analysis: none
ANALYSIS=none

## For static linking uncomment the following
#STATIC= -static
#
ifeq ("$(COMPILER)","gfortran")	
F77= gfortran -fno-automatic 	
## -fbounds-check sometimes causes a weird error due to non-lazy evaluation
## of boolean in gfortran.
FFLAGS= -Wall \
	-ffixed-line-length-132
#	-Wimplicit-interface
## For floating point exception trapping  uncomment the following 
#FPE=-ffpe-trap=invalid,zero 
#,overflow,underflow DONT USE THEM UNLESS THERE IS A SERIOUS REASON
## gfortran 4.4.1 optimized with -O3 yields erroneous results
## Use -O2 to be on the safe side
OPT= -O2 #-O3
## For debugging uncomment the following
DEBUG=  # -pg #-ggdb #-pg -g
endif

ifeq ("$(COMPILER)","g77")
F77= g77 -fno-automatic 
FFLAGS= -Wall -ffortran-bounds-check \
	-ffixed-line-length-132
## For floating point exception trapping  uncomment the following 
#FPEOBJ=trapfpe.o
OPT=-O3
## For debugging uncomment the following
#DEBUG= -ggdb -pg
endif


ifeq ("$(COMPILER)","ifort")
F77 = ifort -save  -extend_source
CXX = icpc
LIBS = -limf
FFLAGS =  -check
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
FF=$(F77) $(FFLAGS) $(FPE) $(OPT) $(DEBUG) -I$(INCLUDE0) -I$(INCLUDE1) -I$(INCLUDE2)


INCLUDE =$(wildcard ../include/*.h *.h include/*.h)

ifeq ("$(PDF)","lhapdf")
LHAPDF_CONFIG=lhapdf-config

ifeq ("$(GRID)","true") 
LHAPDF_CONFIG=/mt/data-grid/***/LHAPDF/install_dir/bin/lhapdf-config
endif

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


ifeq ("$(ANALYSIS)","none")
PWHGANAL=pwhg_analysis-dummy.o 
endif

########## 
# This section is needed to compile and link routines relevant for madgraph:
# add libdhelas3.2a to the list of linked libraries:
#LIBSMADGRAPH=madgraph/DHELAS32/libdhelas3.2.a
#LIBS+=$(LIBSMADGRAPH)
# clear and reassign VPATH (needs also directories relevant to link madgraph):
vpath
VPATH=./:../:obj/:madgraph/:madgraph/real/
# list of madgraph routines:
MADGRAPH=$(shell echo madgraph/real/*.f | sed -e 's/.f/.o/g;s/madgraph\/real\///g') 
MADGRAPH+=functions.o setparam.o $(LIBSMADGRAPH) 
MADGRAPH+=helas.o	
##########

##########
# This is needed to compile and link some source files taken from MCFM
ifeq ("$(VIRTUAL)","MCFM")
VPATH+=MCFM_src
MCFM=$(shell echo MCFM_src/*.f | sed -e 's/\.f/\.o/g;s/MCFM_src\///g') 
endif
##########


ifeq ("$(BH)","true")

BH_CONFIG=$(BHPATH)/bin/blackhat-config
LIBSBLACKHAT = -Wl,-rpath,$(BHPATH)/lib  $(shell $(BH_CONFIG) --libs) -lstdc++
QD_CONFIG=$(QDPATH)/bin/qd-config
LIBSQD = -Wl,-rpath,$(QDPATH)/lib  $(shell $(QD_CONFIG) --libs)

LIBS+=$(LIBSBLACKHAT) $(LIBSQD)

## remember to set LD_LIBRARY_PATH as enviromental variable,
## if the -Wl command is not used when linking with BH
##LD_LIBRARY_PATH=/mt/data-grid/daniel/BH/install_dir/lib:$LD_LIBRARY_PATH
##export LD_LIBRARY_PATH

endif

ifeq ("$(VIRTUAL)","none")
VIRTUALF=virtual.o
endif

ifeq ("$(VIRTUAL)","BH")
VIRTUALF=virtual.o
endif

ifeq ("$(VIRTUAL)","MCFM")
VIRTUALF=virtual_mcfm.o
endif


%.o: %.f $(INCLUDE)
	$(FF) -c -o obj/$@ $<

%.o: %.c
	$(CC) $(DEBUG) -c -o obj/$@ $^ 

%.o: %.cc
	$(CXX) $(DEBUG) -c -o obj/$@ $^ $(FJCXXFLAGS)

#MCFM=$MCFM real_mcfm.o

USER=init_couplings.o init_processes.o Bornquad_planar.o $(VIRTUALF)	\
     real.o $(PWHGANAL) $(MADGRAPH) $(MCFM) Born_phsp.o

PWHG=pwhg_main.o pwhg_init.o bbinit.o btilde.o lhefwrite.o		\
	LesHouches.o LesHouchesreg.o gen_Born_phsp.o find_regions.o	\
	test_Sudakov.o pt2maxreg.o sigborn.o gen_real_phsp.o maxrat.o	\
	gen_index.o gen_radiation.o Bornzerodamp.o sigremnants.o	\
	random.o boostrot.o bra_ket_subroutines.o cernroutines.o	\
	init_phys.o powheginput.o pdfcalls.o sigreal.o sigcollremn.o	\
	pwhg_bookhist.o pwhg_analysis_driver.o checkmomzero.o		\
	setstrongcoupl.o integrator.o newunit.o mwarn.o sigsoftvirt.o	\
	sigcollsoft.o sigvirtual.o reshufflemoms.o $(PDFPACK) $(USER) $(FPEOBJ)

# target to generate LHEF output
pwhg_main: $(PWHG)
	$(FF) $(patsubst %.o,obj/%.o,$(PWHG)) $(LIBS) $(LIBSFASTJET) $(STATIC) -o $@

# target to analyze LHEF output
LHEF=lhef_analysis.o boostrot.o random.o cernroutines.o	      \
	opencount.o powheginput.o pwhg_bookhist.o lhefread.o  \
	newunit.o pwhg_analysis_driver.o $(FPEOBJ)            \
	$(PWHGANAL)

lhef_analysis: $(LHEF)
	$(FF) $(patsubst %,obj/%,$(LHEF)) $(LIBSFASTJET) $(STATIC)  -o $@ 

# target to read event file, shower events with HERWIG + analysis
HERWIG=main-HERWIG.o setup-HERWIG-lhef.o herwig.o boostrot.o powheginput.o \
	pwhg_bookhist.o lhefread.o opencount.o pdfdummies.o $(FPEOBJ)      \
	$(PWHGANAL)

main-HERWIG-lhef: $(HERWIG)
	$(FF) $(patsubst %,obj/%,$(HERWIG)) $(LIBSFASTJET) $(STATIC) -o $@

# target to read event file, shower events with PYTHIA + analysis
PYTHIA=main-PYTHIA.o setup-PYTHIA-lhef.o pythia.o boostrot.o powheginput.o \
	pwhg_bookhist.o lhefread.o opencount.o pdfdummies.o $(FPEOBJ)      \
	$(PWHGANAL)

main-PYTHIA-lhef: $(PYTHIA)
	$(FF) $(patsubst %,obj/%,$(PYTHIA)) $(LIBSFASTJET) $(STATIC) -o $@

# target to cleanup
.PHONY: clean
clean:
	rm -f obj/*.o pwhg_main lhef_analysis main-HERWIG-lhef	\
	main-PYTHIA-lhef \
	madgraph/DHELAS32/*.o madgraph/DHELAS32/*.a

###################################################################
# From here on there are targets specific to the Zjj implementation
##################################################################

##########
# command to compile dhelas3.2 library:
$(LIBSMADGRAPH): force_look
	@echo looking into subdir madgraph/DHELAS32
	$(MAKE) FC='$(F77) $(FFLAGS) ${DEBUG} $(FPE) $(OPT)' --directory=madgraph/DHELAS32
force_look :
	true
##########
