#-*- Makefile -*-
## Choose compiler: gfortran,g77,ifort
COMPILER=gfortran
## Choose PDF: native,lhapdf
## LHAPDF package has to be installed separately
PDF=lhapdf
#Choose Analysis
# default analysis may require FASTJET package, that has to be installed separately (see below)
ANALYSIS=ATLAS
#ANALYSIS=dummy
## For static linking uncomment the following
#STATIC= -static
#
ifeq ("$(COMPILER)","gfortran")	
F77= gfortran
USEDLIBFILES=-Wl,--print-map | grep '^obj-gfortran/libfiles.a' | sed 's/obj-gfortran\/libfiles\.a(// ; s/)//' > libfileslist.txt
## -fbounds-check sometimes causes a weird error due to non-lazy evaluation
## of boolean in gfortran.
#FFLAGS= -Wall -Wimplicit-interface -fbounds-check
## For floating point exception trapping  uncomment the following 
#FPE=-ffpe-trap=invalid,zero,overflow,underflow 
## gfortran 4.4.1 optimized with -O3 yields erroneous results
## Use -O2 to be on the safe side
OPT=-O3
OBJ=obj-gfortran
## For debugging uncomment the following
#DEBUG= -ggdb
ifdef DEBUG
OPT=-O0
#FPE=-ffpe-trap=invalid,zero,overflow
endif
endif

ifeq ("$(COMPILER)","ifort")
F77 = ifort -extend_source
LIBS = -limf
FFLAGS =  #-check
## For floating point exception trapping  uncomment the following 
#FPE = -fpe0
OPT = -O3 #-fast
OBJ=obj-ifort
## For debugging uncomment the following
#DEBUG= -debug -g
ifdef DEBUG
OPT=-O0
endif
endif


PWD=$(shell pwd)
WDNAME=$(shell basename $(PWD))

VPATH= :./:../Version-pre2-1/progress:../Version-pre2-1:../:$(OBJ)/:./MCFMinterface:./MCFM/src/pwgplots/:./MCFM/src/Need/:./MCFM/src/Lib/:./MCFM/src/Phase/:./MCFM/src/W2jet/:./MCFM/src/W2jetvirt/:./MCFM/src/Wbb/:./MCFM/src/Zbb/

INCLUDEC=../Version-pre2-1/include
INCLUDE0=$(PWD)
INCLUDE1=$(PWD)/include
INCLUDE2=$(shell dirname $(PWD))/include
INCLUDE3=$(PWD)/MCFM/src/Inc/
INCLUDE4=./MCFMinterface/Inc
FF=$(F77) $(FFLAGS) $(FPE) $(OPT) $(DEBUG) -I$(INCLUDEC) -I$(INCLUDE0) -I$(INCLUDE1) -I$(INCLUDE2) -I$(INCLUDE3) -I$(INCLUDE4)

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

PWHGANAL=pwhg_analysis-dummy.o pwhg_bookhist-multi.o

ifeq ("$(ANALYSIS)","ATLAS")
#To include Fastjet configuration uncomment the following lines. 
FASTJET_CONFIG=$(shell which fastjet-config)
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_bookhist-multi.o pwhg_analysis-ATLAS.o fastjetfortran.o
endif

PYTHIA8LOCATION=/home/nason/Pheno/pythia8165
FJCXXFLAGS+=-I$(PYTHIA8LOCATION)/include
LIBPYTHIA8=-L$(PYTHIA8LOCATION)/lib/archive -lpythia8 -llhapdfdummy -lstdc++


%.o: %.f
	$(FF) -c -o $(OBJ)/$@ $<

%.o: %.c
	$(CC) $(DEBUG) -c -o $(OBJ)/$@ $^ 

%.o: %.cc
	$(CXX) $(DEBUG) -c -o $(OBJ)/$@ $^ $(FJCXXFLAGS)

USERAUX=conjcolor.o clinkqa.o clinkqga.o isodd.o iseven.o setupBmunu.o 
MCFMfiles=setupckmallowed.o ckmfill.o couplz.o spinoru.o qqb_w2jet.o qqb_wp2jet_g.o \
          qqbw2j_loop.o subqcd.o w2jetsq.o xwqqgg_v.o xwqqggg.o a61LLL.o a61LRL.o \
          a62.o a6routine.o atrLLL.o atrLRL.o atree.o t.o amp_qqggg.o aqqb_zbb.o \
          xzqqgg.o xzqqgg_v.o lnrat.o ampwqqqqg.o a61g.o a6g.o a6treeg.o fax.o \
          faxsl.o fcc.o fcc_qpgmgpqm.o fcc_qpgpgmqm.o fcc_qpgpgpqm.o fcc_qpgpqmgm.o \
          fcc_qpgpqmgp.o fcc_qpqmgmgp.o fcc_qpqmgpgm.o fcc_qpqmgpgp.o fsc.o fsc1.o \
          fsc2.o fsc3.o fsc4.o fsc5.o fsc6.o fsc7.o fsc8.o fvf.o fvs.o vvg.o a6.o \
          a61.o fpm.o fpp.o fsl.o i3m.o lfunctions.o vv.o dclaus.o
USER=$(USERAUX) $(PWHGANAL) init_processes.o init_processes_born_MCFM.o \
     init_processes_real_MCFM.o reorderlord.o reorderreal.o qqb_w2jet_pwhg.o \
     qqb_w2jet_v_pwhg.o qqb_w2jet_g_pwhg.o setborn.o setvirtual.o real.o \
     init_couplings.o Born_phsp.o borncolour_lh.o qqb_w2jet_cs_pwhg.o \
     borncolourqagg.o w2jetsqlc.o fillMCFMCommon.o \
     fillPhysPars.o fillBmunu.o W2jetampmu.o spinormu.o 


PWHG=pwhg_main.o pwhg_init.o bbinit.o btilde.o lhefwrite.o		\
	LesHouches.o LesHouchesreg.o gen_Born_phsp.o find_regions.o	\
	test_Sudakov.o pt2maxreg.o sigborn.o gen_real_phsp.o maxrat.o	\
	gen_index.o gen_radiation.o Bornzerodamp.o sigremnants.o	\
	random.o boostrot.o bra_ket_subroutines.o cernroutines.o 	\
	init_phys.o powheginput.o pdfcalls.o sigreal.o sigcollremn.o	\
	pwhg_analysis_driver.o checkmomzero.o		\
	setstrongcoupl.o integrator.o newunit.o mwarn.o sigsoftvirt.o	\
	sigcollsoft.o sigvirtual.o mint_upb.o opencount.o pwhgreweight.o \
        setlocalscales.o validflav.o  ubprojections.o  \
        $(PDFPACK) $(USER) $(MCFMfiles) $(FPEOBJ)



# target to generate LHEF output
pwhg_main: $(PWHG)
	$(FF) $(patsubst %,$(OBJ)/%,$(PWHG)) $(LIBS) $(LINKMGLIBS) $(LIBSFASTJET) $(STATIC) -o $@

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
PYTHIA=main-PYTHIA.o setup-PYTHIA-lhef.o pythia.o boostrot.o \
        powheginput.o \
	$(PWHGANAL) lhefread.o newunit.o bra_ket_subroutines.o \
	pwhg_analysis_driver.o random.o cernroutines.o opencount.o	\
	$(FPEOBJ)

main-PYTHIA-lhef: $(PYTHIA)
	$(FF) $(patsubst %,$(OBJ)/%,$(PYTHIA)) $(LIBS) $(LIBSFASTJET)  $(STATIC) -o $@

PYTHIA8=main-PYTHIA8.o pythia8F77.o boostrot.o powheginput.o \
	$(PWHGANAL) lhefread.o newunit.o pdfdummies.o		\
	pwhg_analysis_driver.o random.o cernroutines.o opencount.o	\
	$(FPEOBJ)

main-PYTHIA8-lhef: $(PYTHIA8)
	$(FF) $(patsubst %,$(OBJ)/%,$(PYTHIA8)) $(LIBSFASTJET) $(LIBPYTHIA8) $(STATIC) -o $@




# target to cleanup
.PHONY: clean
clean:
	rm -f $(OBJ)/* pwhg_main lhef_analysis main-HERWIG-lhef	\
	main-PYTHIA-lhef

