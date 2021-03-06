#-*- Makefile -*-
############################################### 
#                  Reminder:                  #
# The code in this folder makes use of a      #
# modified sigreal.f file. See comments in    #
# real.f and in the local sigreal.f file.     #
# Keep this file checked when major revisions #
# are made to the repository.                 #
###############################################

## Choose compiler: gfortran,g77,ifort
COMPILER=gfortran
## Choose PDF: native,lhapdf
## LHAPDF package has to be installed separately
PDF=native
#Choose Analysis: none, default
## default analysis may require FASTJET package, that has to be installed separately (see below)
ANALYSIS=none
## For static linking uncomment the following
#STATIC= -static
#
ifeq ("$(COMPILER)","gfortran")	
F77= gfortran -fno-automatic 	
## -fbounds-check sometimes causes a weird error due to non-lazy evaluation
## of boolean in gfortran.
FFLAGS= -Wall -Wimplicit-interface \
	-ffixed-line-length-132
## For floating point exception trapping  uncomment the following 
#FPE=-ffpe-trap=invalid,zero,overflow,underflow 
## gfortran 4.4.1 optimized with -O3 yields erroneous results
## Use -O2 to be on the safe side
OPT=-O2
## For debugging uncomment the following
#DEBUG= -ggdb -pg
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
PWHGANAL=pwhg_analysis.o 
## Also add required Fastjet drivers to PWHGANAL (examples are reported)
PWHGANAL+= fastjetktwrap.o # fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o
else
PWHGANAL=pwhg_analysis-dummy.o
endif

########## 
# This section is needed to compile and link routines
# relevant for madgraph and QCDloop (that needs to be downloaded).
# add libdhelas3.2a to the list of linked libraries:
LIBSMADGRAPH=madgraph/dhelas3.2/libdhelas3.2.a
LIBS+=$(LIBSMADGRAPH)
# add libqcdloop.a and libff.a to the list of linked libraries:
LIBSQCDLOOP= QCDLoop-1.9/ql/libqcdloop.a QCDLoop-1.9/ff/libff.a
LIBS += $(LIBSQCDLOOP)
# clear and reassign VPATH (needs also directories relevant to link madgraph):
vpath
VPATH= ./:../:obj/:madgraph/:madgraph/ME_decayed:QCDLoop-1.9/ql:QCDLoop-1.9/ff
# list of madgraph routines:
MADGRAPH= $(LIBSMADGRAPH) my_setpara.o switchmom.o \
	MAD_born_routines.o \
	MAD_real_routines_DS_tt.o MAD_real_routines_DS_full.o \
	MAD_real_routines_DR.o \
	MAD_born_undec_routines.o MAD_real_undec_routines_DS_full.o \
	MAD_born_dec_routines.o MAD_real_dec_routines.o #printout.o
##########


%.o: %.f $(INCLUDE)
	$(FF) -c -o obj/$@ $<

%.o: %.c
	$(CC) $(DEBUG) -c -o obj/$@ $^ 

%.o: %.cc
	$(CXX) $(DEBUG) -c -o obj/$@ $^ $(FJCXXFLAGS)

USER=init_couplings.o init_processes.o Born_phsp.o Born.o virtual.o	\
     real.o $(PWHGANAL) bmunu2.o boost.o $(MADGRAPH)

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
	$(PWHGANAL) boost.o

lhef_analysis: $(LHEF)
	$(FF) $(patsubst %,obj/%,$(LHEF)) $(LIBSFASTJET) $(STATIC)  -o $@ 

# target to read event file, shower events with HERWIG + analysis
HERWIG=main-HERWIG.o setup-HERWIG-lhef.o herwig.o boostrot.o powheginput.o \
	pwhg_bookhist.o lhefread.o opencount.o pdfdummies.o $(FPEOBJ)      \
	$(PWHGANAL) boost.o

main-HERWIG-lhef: $(HERWIG)
	$(FF) $(patsubst %,obj/%,$(HERWIG)) $(LIBSFASTJET) $(STATIC) -o $@

# target to read event file, shower events with PYTHIA + analysis
PYTHIA=main-PYTHIA.o setup-PYTHIA-lhef.o pythia.o boostrot.o powheginput.o \
	pwhg_bookhist.o lhefread.o opencount.o pdfdummies.o $(FPEOBJ)      \
	$(PWHGANAL) boost.o

main-PYTHIA-lhef: $(PYTHIA)
	$(FF) $(patsubst %,obj/%,$(PYTHIA)) $(LIBSFASTJET) $(STATIC) -o $@

# target to cleanup
.PHONY: clean
clean:
	rm -f obj/*.o pwhg_main lhef_analysis main-HERWIG-lhef	\
	main-PYTHIA-lhef \
	madgraph/dhelas3.2/*.o madgraph/dhelas3.2/*.a \
	merge_t_tb

###################################################################
# From here on there are targets specific to the ST implementation
##################################################################

##########
# to generate the sigreal.f modified version in process subdirectory
sigreal.o: force_look
	@echo "WARNING: For this process the common sigreal.f file cannot be used"
	@echo "WARNING: Creating now a proper local copy of sigreal.f"
	sh make_sigreal.sh
	@echo "WARNING: local copy of sigreal.f created in current subdirectory"
	$(FF) -c -o obj/sigreal.o sigreal.f
##########

##########
# command to compile dhelas3.2 library:
$(LIBSMADGRAPH): force_look
	@echo looking into subdir madgraph/dhelas3.2
	$(MAKE) FC='$(F77) -ffixed-line-length-132' --directory=madgraph/dhelas3.2 
force_look :
	true
##########

##########
# to build executable that merge t and tbar events
# according to relative cross sections:

MERGE_t_tb=merge_t_tb.o newunit.o random.o cernroutines.o

#  target to merge t and tbar samples
merge_t_tb: $(MERGE_t_tb)
	$(FF) $(patsubst %.o,obj/%.o,$(MERGE_t_tb)) -o $@
##########