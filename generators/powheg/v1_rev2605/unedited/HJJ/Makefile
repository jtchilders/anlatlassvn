#-*- Makefile -*-
## Choose compiler: gfortran,ifort (g77 not supported, F90 constructs in use!)
COMPILER=gfortran
FC=$(COMPILER)
## Choose PDF: native,lhapdf
## LHAPDF package has to be installed separately
PDF=lhapdf
#Choose Analysis: none, default
## default analysis may require FASTJET package, that has to be installed separately (see below)
#ANALYSIS=dummy
ANALYSIS=alt
#ANALYSIS=HnJ
#ANALYSIS=none 
## For static linking uncomment the following
#STATIC= -static
#
ifeq ("$(COMPILER)","gfortran")	
F77=gfortran -fno-automatic 	
## -fbounds-check sometimes causes a weird error due to non-lazy evaluation
## of boolean in gfortran.
#FFLAGS= -Wall -Wimplicit-interface -fbounds-check
## For floating point exception trapping  uncomment the following 
#FPE=-ffpe-trap=invalid,zero,overflow
#,underflow 
## gfortran 4.4.1 optimized with -O3 yields erroneous results
## Use -O2 to be on the safe side
OPT=-O2
## For debugging uncomment the following
#DEBUG= -ggdb
ifdef DEBUG
FPE=-ffpe-trap=invalid,zero,overflow
#,underflow
OPT=-O0
endif
endif

PWD=$(shell pwd)
OBJ=$(PWD)/obj-$(COMPILER)
WDNAME=$(shell basename $(PWD))

VPATH= ./:../Version-pre2-1/progress/:../Version-pre2-1/:../:./Madlib/:./MODEL/:./DHELAS/:./HJJvirt/:$(OBJ)/


INCLUDE0=./
INCLUDE0.1=../Version-pre2-1/progress/include
INCLUDE0.2=../Version-pre2-1/include
INCLUDE2=$(shell dirname $(PWD))/include
INCLUDE3=./MCFM_Include 
FF=$(F77) $(FFLAGS) $(FPE) $(OPT) $(DEBUG) -I$(INCLUDE0)  -I$(INCLUDE0.1)  -I$(INCLUDE0.2) -I$(INCLUDE2) -I$(INCLUDE3)

LIBFILES=$(shell  for dir in ./Madlib ./MODEL ./DHELAS ; do cd $$dir ; echo *.[fF] ' ' | sed 's/[fF] /o /g' ; cd .. ; done  )

#INCLUDE =$(wildcard ../include/*.h *.h include/*.h)

ifeq ("$(PDF)","lhapdf")
#LHAPDF_CONFIG=~/Pheno/PDFpacks/lhapdf-5.8.4-$(FC)/bin/lhapdf-config
LHAPDF_CONFIG=lhapdf-config
PDFPACK=lhapdfif.o
LIBSLHAPDF= -Wl,-rpath,$(shell $(LHAPDF_CONFIG) --libdir)  -L$(shell $(LHAPDF_CONFIG) --libdir) -lLHAPDF
ifeq  ("$(STATIC)","-static") 
## If LHAPDF has been compiled with gfortran and you want to link it statically, you have to include
## libgfortran as well. The same holds for libstdc++. 
## One possible solution is to use fastjet, since $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
## does perform this inclusion. The path has to be set by the user. 
# LIBGFORTRANPATH= #/usr/lib/gcc/x86_64-redhat-linux/4.1.2
# LIBSTDCPP=/lib64
# LIBSLHAPDF+=  -L$(LIBGFORTRANPATH)  -lgfortranbegin -lgfortran -L$(LIBSTDCPP) -lstdc++
endif
LIBS+=$(LIBSLHAPDF)
else
PDFPACK=mlmpdfif.o hvqpdfpho.o
endif

# dummy analysis
PWHGANAL=pwhg_bookhist-multi.o pwhg_analysis-dummy.o

ifeq ("$(ANALYSIS)","default")
##To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
#FASTJET_CONFIG=~/lib/fastjet242/bin/fastjet-config
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_bookhist-multi.o pwhg_analysis.o fastjetfortran.o
## Also add required Fastjet drivers to PWHGANAL (examples are reported)
#PWHGANAL+= fastjetsisconewrap.o fastjetktwrap.o fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o
endif

ifeq ("$(ANALYSIS)","HnJ")
##To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
#FASTJET_CONFIG=~/lib/fastjet242/bin/fastjet-config
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_bookhist-multi.o pwhg_analysis-HnJ.o fastjetfortran.o
## Also add required Fastjet drivers to PWHGANAL (examples are reported)#PWHGANAL+= fastjetsisconewrap.o fastjetktwrap.o fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o
endif


ifeq ("$(ANALYSIS)","alt")
##To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
#FASTJET_CONFIG=~/lib/fastjet242/bin/fastjet-config
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_bookhist-multi.o pwhg_analysis-alt-H2J.o fastjetfortran.o
## Also add required Fastjet drivers to PWHGANAL (examples are reported)#PWHGANAL+= fastjetsisconewrap.o fastjetktwrap.o fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o
endif




%.o: %.f $(INCLUDE)
	$(FF) -c -o $(OBJ)/$@ $<

%.o: %.F $(INCLUDE)
	$(FF) -c -o $(OBJ)/$@ $<

%.o: %.c
	$(CC) $(DEBUG) -c -o $(OBJ)/$@ $^ 

%.o: %.cc
	$(CXX) $(DEBUG) -c -o $(OBJ)/$@ $^ $(FJCXXFLAGS)

USER=init_couplings.o init_processes.o Born_phsp.o Born.o virtual.o	\
     real.o matchcolour.o $(PWHGANAL)

H2JVIRTFILES=Ampvirt_AQgg.o GZHggggvsqPoles.o HAQggvsqanal.o Hqaggvsq.o Hqaqavsqanal.o\
	Hqarbvsqanal.o	checkegzres.o	gg_hgg_vi.o Ampvirt_gggg.o GZHqaggvsqPoles.o\
	Hggggvsq.o Hqaqasq.o Hqarbsq.o checkscheme.o qqgghn.o\
	Ampvirtsq_AQaq.o H4prenorm.o Hggggvsqanal.o Hqaqavsq.o Hqarbvsq.o amplonumer.o\
	gg_hgg_eval.o H4plo.o lfunctions.o spinoru.o h4g.o h4gnew.o h4q.o dclaus.o hqqggdfm.o\
	hqqgg.o Amplo_AQgg.o i3m.o dot.o\

H2JVFILES= A0HAQgg.o A0phiAQgg.o A0phiggggpppp.o Lsm1DS.o A1phiAQggmpmm.o A1phiAgQgmppp.o\
	A41phiAQgg.o A43phiAQggmpmp.o Aslcphiqarbmpmp.o F41mF.o A0Hggggmmpp.o A0phiAgQg.o\
	A0phiqarb.o A1phiAQggmpmp.o A1phiggggmmmm.o A41phiqarb.o A43phiAQggmppm.o Aslcphiqarbmppm.o\
	F42me.o A0Hggggmpmp.o A0phiggggmmmm.o A1Hggggmmpp.o A1phiAQggmppm.o A1phiggggmmpp.o A42Hqarb.o Afphiqarbmpmp.o\
	BGRL.o F42meF.o A0Hggggpmmm.o A0phiggggmmpp.o A1Hggggmpmp.o A1phiAQggmppp.o A1phiggggmpmp.o A42phiqarb.o\
	Afphiqarbmppm.o F31m.o F42mhF.o A0Hggggpppp.o A0phiggggmpmp.o A1Hggggpmmm.o A1phiAgQgmmpp.o A41HAQgg.o\
	A43HAQgg.o Alcphiqarbmpmp.o F33m.o A0Hqarb.o A0phiggggpmmm.o A1Hggggpppp.o A1phiAgQgmppm.o A41Hqarb.o A43phiAQggmpmm.o\
	Alcphiqarbmppm.o F41m.o lnrat.o \


MCFMFILES=i2mcfm_2_POWHEG.o set_interface_MCFM.o i2mcfm_fill.o \
	dotem.o set_epinv.o coupling.o


PWHG=pwhg_main.o pwhg_init.o bbinit.o btilde.o lhefwrite.o		\
	LesHouches.o LesHouchesreg.o gen_Born_phsp.o find_regions.o	\
	test_Sudakov.o pt2maxreg.o sigborn.o gen_real_phsp.o maxrat.o	\
	gen_index.o gen_radiation.o Bornzerodamp.o sigremnants.o	\
	random.o boostrot.o bra_ket_subroutines.o cernroutines.o	\
	init_phys.o powheginput.o pdfcalls.o sigreal.o sigcollremn.o	\
	pwhg_analysis_driver.o checkmomzero.o pwhgreweight.o            \
	setstrongcoupl.o integrator.o newunit.o mwarn.o sigsoftvirt.o	\
	reshufflemoms.o setlocalscales.o validflav.o mint_upb.o         \
	sigcollsoft.o sigvirtual.o opencount.o ubprojections.o          \
        $(PDFPACK) $(USER)  $(FPEOBJ) $(MCFMFILES) $(H2JVIRTFILES) $(H2JVFILES)

LIBDIRMG=.
LINKMGLIBS =  -L$(LIBDIRMG)  -lmadgraph -lmodel -ldhelas3 

MADLIBS=libdhelas3.a libmadgraph.a libmodel.a

# target to generate LHEF output
pwhg_main:$(PWHG) libfiles.a
	$(FF) $(patsubst %,$(OBJ)/%,$(PWHG)) $(OBJ)/libfiles.a $(LIBS) $(LIBSFASTJET) $(STATIC) -o $@


libfiles.a: $(LIBFILES)
	cd $(OBJ) ; \rm libfiles.a ; ar cru libfiles.a $(LIBFILES)

LHEF=lhef_analysis.o boostrot.o random.o cernroutines.o		\
     opencount.o powheginput.o $(PWHGANAL)	\
     lhefread.o newunit.o pwhg_analysis_driver.o bra_ket_subroutines.o $(FPEOBJ)

# target to analyze LHEF output
lhef_analysis:$(LHEF)
	$(FF) $(patsubst %,$(OBJ)/%,$(LHEF)) $(LIBS) $(LIBSFASTJET) $(STATIC)  -o $@ 



# target to read event file, shower events with HERWIG + analysis
HERWIG=main-HERWIG.o setup-HERWIG-lhef.o herwig.o boostrot.o	\
	powheginput.o $(PWHGANAL) lhefread.o	\
	pdfdummies.o opencount.o $(FPEOBJ) 

main-HERWIG-lhef: $(HERWIG)
	$(FF) $(patsubst %,$(OBJ)/%,$(HERWIG))  $(LIBSFASTJET)  $(STATIC) -o $@

# target to read event file, shower events with PYTHIA + analysis
PYTHIA=main-PYTHIA.o setup-PYTHIA-lhef.o pythia.o boostrot.o powheginput.o		\
	$(PWHGANAL) lhefread.o newunit.o pdfdummies.o	\
	pwhg_analysis_driver.o random.o cernroutines.o opencount.o	\
	$(FPEOBJ)

main-PYTHIA-lhef: $(PYTHIA)
	$(FF) $(patsubst %,$(OBJ)/%,$(PYTHIA)) $(LIBSFASTJET)  $(STATIC) -o $@

# target to cleanup
.PHONY: clean veryclean $(OBJ)

$(OBJ):
	if ! [ -d $(OBJ) ] ; then mkdir $(OBJ) ; fi

clean:
	rm -f pwhg_main lhef_analysis main-HERWIG-lhef main-PYTHIA-lhef ; \
        cd $(OBJ) ; rm -f $(PWHG) $(HERWIG) $(PYTHIA) $(LHEF)


veryclean:
	cd $(OBJ) ; \rm *


ifeq ("$(COMPILER)","gfortran")
XFFLAGS +=-ffixed-line-length-132
else
XFFLAGS +=-extend-source
endif
