#-*- Makefile -*-
## Choose compiler: gfortran,g77,ifort
COMPILER=gfortran
## Choose PDF: native,lhapdf
## LHAPDF package has to be installed separately
PDF=lhapdf
#Choose Analysis: none, default
## default analysis may require FASTJET package, that has to be installed separately (see below)
ANALYSIS=default5f
## For static linking uncomment the following
#STATIC=-static
#
ifeq ("$(COMPILER)","gfortran")	
F77= gfortran -fno-automatic 	
## -fbounds-check sometimes causes a weird error due to non-lazy evaluation
## of boolean in gfortran.
FFLAGS= -Wall \
	-ffixed-line-length-132 #\
		-fbounds-check
## For floating point exception trapping  uncomment the following 
FPE= #-ffpe-trap=invalid,zero #,overflow,underflow 
## gfortran 4.4.1 optimized with -O3 yields erroneous results
## Use -O2 to be on the safe side
OPT= -O2
## For debugging uncomment the following
DEBUG= #-ggdb #-pg
endif

ifeq ("$(COMPILER)","g77")
F77= g77 -fno-automatic 
FFLAGS= -Wall -ffortran-bounds-check \
	-ffixed-line-length-132
## For floating point exception trapping  uncomment the following 
#FPEOBJ=trapfpe.o
OPT=-O2
## For debugging uncomment the following
#DEBUG= -ggdb -pg
endif


ifeq ("$(COMPILER)","ifort")
F77 = /usr/local/intel/fc/11.0.74/bin/ia32/ifort -save  -extend_source
#CXX = icpc
LIBS = #-limf
FFLAGS =  #-check
## For floating point exception trapping  uncomment the following 
#FPE = -fpe0
OPT = #-O3 #-fast
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
LHAPDF_CONFIG=$(shell which lhapdf-config)
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


ifeq ("$(ANALYSIS)","dummy")
PWHGANAL=pwhg_analysis-dummy.o
else
ifeq ("$(ANALYSIS)","default")
##To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_analysis.o 
## Also add required Fastjet drivers to PWHGANAL (examples are reported)
PWHGANAL+= fastjetktwrap.o # fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o
endif
ifeq ("$(ANALYSIS)","default5f")
##To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_analysis5f.o 
## Also add required Fastjet drivers to PWHGANAL (examples are reported)
PWHGANAL+= fastjetktwrap.o # fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o
endif
ifeq ("$(ANALYSIS)","NLOvalidation")
#To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_analysis_NLOvalidation.o 
## Also add required Fastjet drivers to PWHGANAL (examples are reported)
PWHGANAL+= fastjetktwrap.o # fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o
endif
ifeq ("$(ANALYSIS)","Bjet")
#To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_analysis_Bjet.o 
## Also add required Fastjet drivers to PWHGANAL (examples are reported)
PWHGANAL+= fastjetktwrap.o # fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o
endif
ifeq ("$(ANALYSIS)","Bjetbugfix")
#To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_analysis_Bjet_bugfix.o 
## Also add required Fastjet drivers to PWHGANAL (examples are reported)
PWHGANAL+= fastjetktwrap.o # fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o
endif
ifeq ("$(ANALYSIS)","BjetbugfixalsoNLO")
#To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_analysis_Bjet_bugfix_alsoNLO.o 
## Also add required Fastjet drivers to PWHGANAL (examples are reported)
PWHGANAL+= fastjetktwrap.o # fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o
endif
ifeq ("$(ANALYSIS)","Bjet-new")
#To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_analysis_Bjet-new.o 
## Also add required Fastjet drivers to PWHGANAL (examples are reported)
PWHGANAL+= fastjetktwrap.o # fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o
endif

ifeq ("$(ANALYSIS)","Bjet-realistic")
#To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_analysis_Bjet-realistic.o 
## Also add required Fastjet drivers to PWHGANAL (examples are reported)
PWHGANAL+= fastjetktwrap.o # fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o
endif

ifeq ("$(ANALYSIS)","hwanst4f_v3")
#To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=powheg_hwanst4f_v3.o 
## Also add required Fastjet drivers to PWHGANAL (examples are reported)
PWHGANAL+= fastjetktwrap.o # fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o
endif

ifeq ("$(ANALYSIS)","hwanst4f_v7")
#To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=powheg_hwanst4f_v7.o 
## Also add required Fastjet drivers to PWHGANAL (examples are reported)
PWHGANAL+= fastjetktwrap.o # fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o
endif

ifeq ("$(ANALYSIS)","hwanst4f_v8")
#To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=powheg_hwanst4f_v8.o 
## Also add required Fastjet drivers to PWHGANAL (examples are reported)
PWHGANAL+= fastjetktwrap.o # fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o
endif


endif



########## 
# This section is needed to compile and link routines relevant for madgraph:
# add libdhelas3.2a to the list of linked libraries:
#LIBSMADGRAPH=madgraph/DHELAS32/libdhelas3.2.a
#LIBS+=$(LIBSMADGRAPH)
# clear and reassign VPATH (needs also directories relevant to link madgraph):
vpath
VPATH=./:../:obj/:madgraph/:madgraph/real/:madgraph/born/:madgraph/real_tb/:madgraph/born_tb/
# list of madgraph routines:
MADGRAPH=$(shell echo madgraph/born/*.f | sed -e 's/.f/.o/g;s/madgraph\/born\///g') 
MADGRAPH+=$(shell echo madgraph/real/*.f | sed -e 's/.f/.o/g;s/madgraph\/real\///g') 
MADGRAPH+=$(shell echo madgraph/born_tb/*.f | sed -e 's/.f/.o/g;s/madgraph\/born_tb\///g') 
MADGRAPH+=$(shell echo madgraph/real_tb/*.f | sed -e 's/.f/.o/g;s/madgraph\/real_tb\///g') 
MADGRAPH+=functions.o setparam.o $(LIBSMADGRAPH) 
MADGRAPH+=helas.o	
##########
# This is needed to compile and link some source files taken from MCFM
VPATH+=MCFM_src
MCFM=$(shell echo MCFM_src/*.f | sed -e 's/\.f/\.o/g;s/MCFM_src\///g') 
##########




%.o: %.f $(INCLUDE)
	$(FF) -c -o obj/$@ $<

%.o: %.c
	$(CC) $(DEBUG) -c -o obj/$@ $^ 

%.o: %.cc
	$(CXX) $(DEBUG) -c -o obj/$@ $^ $(FJCXXFLAGS)

USER=init_couplings.o init_processes.o Born_phsp.o Born.o virtual.o	\
     real.o $(PWHGANAL) $(MADGRAPH) bmunu2.o $(MCFM)

PWHG=pwhg_main.o pwhg_init.o bbinit.o btilde.o lhefwrite.o		\
	LesHouches.o LesHouchesreg.o gen_Born_phsp.o find_regions.o	\
	test_Sudakov.o pt2maxreg.o sigborn.o gen_real_phsp.o maxrat.o	\
	gen_index.o gen_radiation.o Bornzerodamp.o sigremnants.o	\
	random.o boostrot.o bra_ket_subroutines.o cernroutines.o	\
	init_phys.o powheginput.o pdfcalls.o sigreal.o sigcollremn.o	\
	pwhg_bookhist.o pwhg_analysis_driver.o checkmomzero.o		\
	setstrongcoupl.o integrator.o newunit.o mwarn.o sigsoftvirt.o	\
	sigcollsoft.o sigvirtual.o  $(PDFPACK) $(USER) $(FPEOBJ)

# target to generate LHEF output
pwhg_main: $(PWHG)
	$(FF) $(patsubst %.o,obj/%.o,$(PWHG)) $(LIBS) $(LIBSFASTJET) $(STATIC) -o $@

# target to analyze LHEF output
LHEF=lhef_analysis.o boostrot.o random.o cernroutines.o	      \
	opencount.o powheginput.o pwhg_bookhist.o lhefread.o  \
	newunit.o pwhg_analysis_driver.o $(FPEOBJ)            \
	$(PWHGANAL) #boost.o

lhef_analysis: $(LHEF)
	$(FF) $(patsubst %,obj/%,$(LHEF)) $(LIBSFASTJET) $(STATIC)  -o $@ 

# target to read event file, shower events with HERWIG + analysis
HERWIG=main-HERWIG.o setup-HERWIG-lhef.o herwig.o boostrot.o powheginput.o \
	pwhg_bookhist.o lhefread.o opencount.o $(FPEOBJ)      \
	$(PWHGANAL) #boost.o #pdfdummies.o #running herwig with external pdf!

main-HERWIG-lhef: $(HERWIG)
	$(FF) $(patsubst %,obj/%,$(HERWIG)) $(LIBS) $(LIBSFASTJET) $(STATIC) -o $@

# target to read event file, shower events with PYTHIA + analysis
PYTHIA=main-PYTHIA.o setup-PYTHIA-lhef.o pythia.o boostrot.o powheginput.o \
	pwhg_bookhist.o lhefread.o opencount.o pdfdummies.o $(FPEOBJ)      \
	$(PWHGANAL) #boost.o

main-PYTHIA-lhef: $(PYTHIA)
	$(FF) $(patsubst %,obj/%,$(PYTHIA)) $(LIBSFASTJET) $(STATIC) -o $@

#############################################
# target to read event file, shower events with HERWIG6520 + analysis
HERWIG6520=main-HERWIG.o setup-HERWIG-lhef.o herwig6520.o boostrot.o powheginput.o \
	pwhg_bookhist.o lhefread.o opencount.o $(FPEOBJ)      \
	$(PWHGANAL) pdfdummies.o #boost.o #running herwig with external pdf!

main-HERWIG6520-lhef: $(HERWIG6520)
	$(FF) $(patsubst %,obj/%,$(HERWIG6520)) $(LIBS) $(LIBSFASTJET) $(STATIC) -o $@

# target to read event file, shower events with PYTHIA-TUNED+MPI + analysis
PYTHIA_tune=main-PYTHIA.o setup-PYTHIA_tune-lhef.o pythia.o boostrot.o powheginput.o \
	pwhg_bookhist.o lhefread.o opencount.o pdfdummies.o $(FPEOBJ)      \
	$(PWHGANAL) #boost.o

main-PYTHIA_tune-lhef: $(PYTHIA_tune)
	$(FF) $(patsubst %,obj/%,$(PYTHIA_tune)) $(LIBSFASTJET) $(STATIC) -o $@

#############################################

# target to cleanup
.PHONY: clean
clean:
	rm -f obj/*.o pwhg_main lhef_analysis main-HERWIG-lhef	\
	main-PYTHIA-lhef main-HERWIG6520-lhef \
	madgraph/DHELAS32/*.o madgraph/DHELAS32/*.a \
	merge_t_tb

###################################################################
# From here on there are targets specific to this implementation
##################################################################

##########
# # command to compile dhelas3.2 library:
# $(LIBSMADGRAPH): force_look
# 	@echo looking into subdir madgraph/DHELAS32
# 	$(MAKE) FC='$(F77) $(FFLAGS) $(FPE) $(OPT)' --directory=madgraph/DHELAS32
# force_look :
# 	true
##########



##########
# to build executable that merge t and tbar events
# according to relative cross sections:

MERGE_t_tb=merge_t_tb.o newunit.o random.o cernroutines.o

#  target to merge t and tbar samples
merge_t_tb: $(MERGE_t_tb)
	$(FF) $(patsubst %.o,obj/%.o,$(MERGE_t_tb)) -o $@
##########