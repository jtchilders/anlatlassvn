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
ANALYSIS=ATLAS
## For static linking uncomment the following
#STATIC= -static
#
ifeq ("$(COMPILER)","gfortran")	
F77= gfortran
USEDLIBFILES=-Wl,--print-map | grep '^obj-gfortran/libfiles.a' | sed 's/obj-gfortran\/libfiles\.a(// ; s/)//' > libfileslist.txt
OBJ=obj-gfortran
## -fbounds-check sometimes causes a weird error due to non-lazy evaluation
## of boolean in gfortran.
#FFLAGS= -Wall -Wimplicit-interface -fbounds-check
## For floating point exception trapping  uncomment the following 
#FPE=-ffpe-trap=invalid,zero,overflow,underflow 
## gfortran 4.4.1 optimized with -O3 yields erroneous results
## Use -O2 to be on the safe side
OPT=-O3
## For debugging uncomment the following
#DEBUG= -ggdb
ifdef DEBUG
OPT=-O0 -pg
FPE=-ffpe-trap=invalid,zero,overflow,underflow
endif
endif

ifeq ("$(COMPILER)","ifort")
F77 = ifort   -extend_source
OBJ=obj-ifort
## For floating point exception trapping  uncomment the following 
#FPE = -fpe0
OPT = -O3 
#-O3 #-fast
## For debugging uncomment the following
#DEBUG= -debug -g
ifdef DEBUG
OPT=-O0
endif
endif

LINKMGLIBS = 

PWD=$(shell pwd)
WDNAME=$(shell basename $(PWD))

DEPDIRS=./MCFMinterface/ ./MCFM/src/pwgplots/ ./MCFM/src/Need/ \
./MCFM/src/Lib/ ./MCFM/src/Phase/ ./MCFM/src/Z2jet/ \
./MCFM/src/Zbb/  ./MCFM/src/W2jet/ ./MCFM/src/W2jetvirt/ ./MCFM/src/Wbb/ 

VPATH= ./:../Version-pre2-1/progress:../Version-pre2-1:../:$(OBJ)/:./MCFMinterface:./MCFM/src/pwgplots/:./MCFM/src/Need/:./MCFM/src/Lib/:./MCFM/src/Phase/:./MCFM/src/Z2jet/:./MCFM/src/Zbb/:./MCFM/src/W2jet/:./MCFM/src/W2jetvirt/:./MCFM/src/Wbb/


INCLUDEA=../Version-pre2-1/progress/include
INCLUDEC=../Version-pre2-1/include
INCLUDE0=$(PWD)
INCLUDE1=$(PWD)/include
INCLUDE2=$(shell dirname $(PWD))/include
INCLUDE3=$(PWD)/MCFM/src/Inc/
INCLUDE4=./MCFMinterface/Inc
FF=$(F77) $(FFLAGS) $(FPE) $(OPT) $(DEBUG) -I$(INCLUDEA) -I$(INCLUDEC) -I$(INCLUDE0) -I$(INCLUDE1) -I$(INCLUDE2) -I$(INCLUDE3) -I$(INCLUDE4)



#LIBFILES=$(shell  for dir in $(DEPDIRS) ; do if cd $$dir ; then echo *.[fF] ' ' | sed 's/coupling.f//g ; s/reader.f//g ; s/mdata_madgraph.f//g ; s/dipoles_mass_old.f//g ; s/initialize3.f//g ; s/alfacut.f//g ; s/redconstants.f//g ; s/ eq.f//g ; s/[fF] /o /g' ; fi ; cd $(PWD) ; done  )


ifeq ("$(PDF)","lhapdf")
LHAPDF_CONFIG=lhapdf-config
PDFPACK=lhapdfif.o
LIBSLHAPDF= -Wl,-rpath,$(shell $(LHAPDF_CONFIG) --libdir)  -L$(shell $(LHAPDF_CONFIG) --libdir) -lLHAPDF
ifeq  ("$(STATIC)","-static") 
## If LHAPDF has been compiled with gfortran and you want to link it statically, you have to include
## libgfortran as well. The same holds for libstdc++. 
## One possible solution is to use fastjet, since $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
## does perform this inclusion. The path has to be set by the user. 
# LIBSTDCPP=/lib64
# LIBSLHAPDF+= -lgfortranbegin -lgfortran -L$(LIBSTDCPP) 
# -lstdc++
endif
LIBS+=$(LIBSLHAPDF)
else
PDFPACK=mlmpdfif.o hvqpdfpho.o
endif

PWHGANAL=pwhg_bookhist-multi.o pwhg_analysis-dummy.o

ifeq ("$(ANALYSIS)","ATLAS")
##To include Fastjet configuration uncomment the following lines. 
#FASTJET_CONFIG=/users/melia/fastj/bin/fastjet-config
FASTJET_CONFIG=$(shell which fastjet-config)
#FASTJET_CONFIG=~/lib/fastjet242/bin/fastjet-config
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_bookhist-multi.o pwhg_analysis-ATLAS.o fastjetfortran.o
## Also add required Fastjet drivers to PWHGANAL (examples are reported)
#PWHGANAL+= fastjetsisconewrap.o fastjetktwrap.o fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o
endif

#PWHGANAL=./pwhg_analysis.o

%.o: %.f
	$(FF) -c -o $(OBJ)/$@ $<

%.o: %.c
	$(CC) $(DEBUG) -c -o $(OBJ)/$@ $^ 

%.o: %.cc
	$(CXX) $(DEBUG) -c -o $(OBJ)/$@ $^ $(FJCXXFLAGS)

USERAUX=conjcolor.o clinkqa.o clinkqga.o isodd.o iseven.o setupBmunu.o 
USER=$(USERAUX) $(PWHGANAL) init_processes.o init_processes_born_MCFM.o init_processes_real_MCFM.o \
     setborn.o setvirtual.o real.o init_couplings.o Born_phsp.o borncolour_lh.o borncolourqagg.o \
     reorderZ0lord.o reorderZ0real.o qqb_z2jet_pwhg.o qqb_z2jet_v_pwhg.o qqb_z2jet_g_pwhg.o \
     qqb_z2jet_cs_pwhg.o compare.o comparejk.o \
     fillMCFMCommon.o Z2jetampmu.o fillZBmunu.o 


MCFMfiles=spinormu.o couplz.o spinoru.o a61z.o a62z.o a63z.o atreez.o msq_ZqqQQg.o storecsz.o z2jetsq.o \
          ampqqb_qqb.o aqqb_zbb.o xzqqgg_v.o xzqqggg.o subqcd.o a61g.o a6g.o a6treeg.o fax.o faxsl.o fcc.o \
          fcc_qpgmgpqm.o fcc_qpgpgmqm.o fcc_qpgpgpqm.o fcc_qpgpqmgm.o fcc_qpgpqmgp.o fcc_qpqmgmgp.o fcc_qpqmgpgm.o \
          fcc_qpqmgpgp.o fsc.o fsc1.o fsc2.o fsc3.o fsc4.o fsc5.o fsc6.o fsc7.o fsc8.o fvf.o fvs.o vvg.o a61.o a62.o \
          a6routine.o atree.o i3m.o lfunctions.o t.o lnrat.o dclaus.o a63.o a6ax.o fmt.o fzip.o makemb.o amp_qqggg.o \
          xzqqgg.o a6.o fpm.o fpp.o fsl.o vv.o makem.o nagyqqqqg.o

PWHG=pwhg_main.o pwhg_init.o bbinit.o btilde.o lhefwrite.o		\
	LesHouches.o LesHouchesreg.o gen_Born_phsp.o find_regions.o	\
	test_Sudakov.o pt2maxreg.o sigborn.o gen_real_phsp.o maxrat.o	\
	gen_index.o gen_radiation.o Bornzerodamp.o sigremnants.o	\
	random.o boostrot.o bra_ket_subroutines.o cernroutines.o 	\
	init_phys.o powheginput.o pdfcalls.o sigreal.o sigcollremn.o	\
	pwhg_analysis_driver.o checkmomzero.o 				\
	setstrongcoupl.o integrator.o newunit.o mwarn.o sigsoftvirt.o	\
	sigcollsoft.o sigvirtual.o mint_upb.o opencount.o pwhgreweight.o \
        setlocalscales.o validflav.o ubprojections.o \
        $(PDFPACK) $(USER) $(MCFMfiles) $(FPEOBJ)



# target to generate LHEF output
#pwhg_main:  libfiles.a $(PWHG)
#	$(FF) $(patsubst %,$(OBJ)/%,$(PWHG)) $(LIBS)  $(OBJ)/libfiles.a $(LIBSFASTJET) $(STATIC) -o $@ $(USEDLIBFILES)

pwhg_main:  $(PWHG) 
	$(FF) $(patsubst %,$(OBJ)/%,$(PWHG)) $(LIBS) $(LIBSFASTJET) $(STATIC) -o $@

#libfiles.a: $(LIBFILES)
#	cd $(OBJ) ; \rm libfiles.a ; ar cru libfiles.a $(LIBFILES)


LHEF=lhef_analysis.o boostrot.o random.o cernroutines.o		\
     opencount.o powheginput.o $(PWHGANAL)	\
     lhefread.o newunit.o pwhg_analysis_driver.o $(FPEOBJ)

# target to analyze LHEF output
lhef_analysis:$(LHEF)
	$(FF) $(patsubst %,$(OBJ)/%,$(LHEF)) $(LIBS) $(LIBSFASTJET) $(STATIC)  -o $@ 



# target to read event file, shower events with HERWIG + analysis
HERWIG=main-HERWIG.o setup-HERWIG-lhef.o herwig.o powheginput.o\
	$(PWHGANAL) lhefread.o opencount.o pdfdummies.o $(FPEOBJ)

main-HERWIG-lhef: $(HERWIG)
	$(FF) $(patsubst %,$(OBJ)/%,$(HERWIG))  $(LIBSFASTJET)  $(STATIC) -o $@

# target to read event file, shower events with PYTHIA + analysis
PYTHIA=main-PYTHIA.o setup-PYTHIA-lhef.o pythia.o boostrot.o powheginput.o \
	$(PWHGANAL) lhefread.o newunit.o bra_ket_subroutines.o \
	pwhg_analysis_driver.o random.o cernroutines.o opencount.o	\
	$(FPEOBJ)

main-PYTHIA-lhef: $(PYTHIA)
	$(FF) $(patsubst %,$(OBJ)/%,$(PYTHIA)) $(LIBSFASTJET) $(LIBS) $(STATIC) -o $@

# target to read event file, shower events with PYTHIA + analysis

# target to cleanup
.PHONY: clean
clean:
	rm -f $(OBJ)/*.o pwhg_main lhef_analysis main-HERWIG-lhef	\
	main-PYTHIA-lhef

