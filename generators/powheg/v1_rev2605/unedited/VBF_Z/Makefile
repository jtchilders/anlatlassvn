#-- Makefile to be used for VBF pp->Zjj

## Choose PDF: native,lhapdf
## LHAPDF package has to be installed separately
PDF=lhapdf
#PDF=native

# Choose Analysis: none, default
# default analysis requires FASTJET package, that has to be installed 
# separately (see below)
#ANALYSIS=none
ANALYSIS=default
## For static linking uncomment the following
#STATIC= -static

F77= gfortran -fno-automatic -ffixed-line-length-none
OPT=-O2

PWD=$(shell pwd)
WDNAME=$(shell basename $(PWD))
VPATH= ./:../:obj/:mg_zfiles:

INCLUDE0=$(PWD)
INCLUDE1=$(PWD)/include
INCLUDE2=$(shell dirname $(PWD))/include
FF=$(F77) $(FFLAGS) $(FPE) $(OPT) $(DEBUG) -I$(INCLUDE0) -I$(INCLUDE1) -I$(INCLUDE2) -I nlo_files


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
# LIBGFORTRANPATH=/usr/lib/gcc/x86_64-redhat-linux/4.1.2
# LIBSTDCPP=/lib64
 LIBSLHAPDF+=  -L$(LIBGFORTRANPATH)  -lgfortranbegin -lgfortran -L$(LIBSTDCPP) -lstdc++
endif
LIBS+=$(LIBSLHAPDF)
else
PDFPACK=mlmpdfif.o hvqpdfpho.o
endif

ifeq ("$(ANALYSIS)","default")
FASTJET_CONFIG=$(shell which fastjet-config)
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_analysis.o fastjetfortran.o
# Also add required Fastjet drivers to PWHGANAL (examples are reported)
PWHGANAL+= fastjetsisconewrap.o fastjetktwrap.o #fastjetCDFMidPointwrap.o fastjetD0RunIIConewrap.o fastjetfortran.o
else
PWHGANAL=pwhg_analysis-dummy.o
endif

LIBS+=$(LHAPDFLIB) -Lnlo_files -lnlo 

%.o: %.f $(INCLUDE)
	$(FF) -c -o obj/$@ $<

%.o: %.c
	$(CC) $(DEBUG) -c -o obj/$@ $^ 

%.o: %.cc
	$(CC) $(DEBUG) -c -o obj/$@ $^ $(FJCXXFLAGS)


USER=init_couplings.o init_processes.o Born_phsp.o Born.o virtual.o  \
	calc_invariants.o compute_tensors_z.o compute_tensors_z_real.o \
	convert_coup.o real.o boost.o breit.o  $(PWHGANAL)

PWHG=pwhg_main.o pwhg_init.o bbinit.o btilde.o lhefwrite.o            \
	LesHouches.o LesHouchesreg.o gen_Born_phsp.o find_regions.o   \
        test_Sudakov.o pt2maxreg.o sigborn.o gen_real_phsp_fix.o maxrat.o \
	gen_index.o gen_radiation.o Bornzerodamp.o sigremnants.o      \
	random.o boostrot.o bra_ket_subroutines.o cernroutines.o      \
	init_phys.o powheginput.o pdfcalls.o sigreal.o sigcollremn.o  \
        pwhg_bookhist.o pwhg_analysis_driver.o checkmomzero.o         \
	reshufflemoms.o setstrongcoupl.o integrator.o newunit.o       \
	mwarn.o sigsoftvirt.o sigcollsoft.o sigvirtual.o  \
	$(PDFPACK) $(USER) $(SYSOBJ)

# target to generate LHEF output
pwhg_main:$(PWHG) nlo_files/libnlo.a 
	(cd nlo_files; make;)
	$(FF) $(patsubst %,obj/%, $(PWHG)) $(LIBS) $(LIBSFASTJET) -o $@


nlo_files/libnlo.a: nlo_files/*f 
	cd nlo_files ; make

LHEFANAL=lhef_analysis.o boostrot.o random.o cernroutines.o         \
	 opencount.o powheginput.o pwhg_bookhist.o  \
	 lhefread.o newunit.o pwhg_analysis_driver.o \
	 $(PWHGANAL) $(SYSOBJ)


# target to analyze LHEF output
lhef_analysis:$(LHEFANAL)
	$(FF) $(patsubst %,obj/%,$(LHEFANAL)) -o $@ $(LIBS) $(LIBSFASTJET)


# target to read event file, shower events with HERWIG + analysis
HERWIG=main-HERWIG.o setup-HERWIG-lhef.o herwig6510.o \
	boostrot.o powheginput.o \
        pwhg_bookhist.o lhefread.o pdfdummies.o newunit.o \
	pwhg_analysis_driver.o \
	random.o cernroutines.o opencount.o  $(PWHGANAL) $(SYSOBJ)


main-HERWIG-lhef: $(HERWIG)
	$(FF) $(patsubst %,obj/%,$(HERWIG)) -o $@ $(LIBSFASTJET)

# target to read event file, shower events with PYTHIA + analysis
PYTHIA=main-PYTHIA.o setup-PYTHIA-lhef.o pythia-6.4.25.o \
	boostrot.o powheginput.o \
        pwhg_bookhist.o lhefread.o newunit.o \
	pwhg_analysis_driver.o pdfdummies.o \
	random.o cernroutines.o opencount.o  $(PWHGANAL) $(SYSOBJ)


main-PYTHIA-lhef: $(PYTHIA)
	$(FF) $(patsubst %,obj/%,$(PYTHIA)) -o $@ $(LIBS) $(LIBSFASTJET) $(LIBS)

clean:
	rm -f *.o *~ obj_ew/*.o obj/*.o \
        nlo_files/*.o nlo_files/libnlo.a \
        nlo_files/*/*.o \
        pwhg_main lhef_analysis main-HERWIG-lhef \
	main-PYTHIA-lhef
