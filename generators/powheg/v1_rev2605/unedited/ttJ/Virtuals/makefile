ifndef COMPILER
COMPILER=gnu
endif

OBJSDIR = obj-$(COMPILER)
IFLAGS = -Iauto/ -IReduction-1.1/ -Iff/ -Iinclude/ 
PWD=$(shell pwd)

ifeq ("$(COMPILER)", "intel")
  CXX = icpc
  CC = icc
  FC = ifort
  SETCXXOUT = -c -o
  SETLDOUT = -o
  LFLAGS = -lifcore -lirc -limf  -lpthread  -lifport -lifcoremt
  FFLAGS := $(FFLAGS) -132 -fPIC
  CFLAGS := $(CFLAGS) -ftz  	
  LIBDIR := -L/opt/products/ifort/11.0/lib/intel64/
endif

ifeq ("$(COMPILER)", "gnu")
  CXX = g++
  CC = gcc
  FC = gfortran
  SETCXXOUT = -c -o
  SETLDOUT = -o
  LFLAGS = -lgfortranbegin -lgfortran -lm
  FFLAGS := $(FFLAGS) -fno-second-underscore -ffixed-line-length-132 -fPIC
endif

CFLAGS := $(CFLAGS) $(DEBUG) -fPIC 
FFLAGS := $(FFLAGS) $(DEBUG)
IFLAGS := -Iauto -Imadgraph $(IFLAGS)
SETCXXOUT = -c -o

LFLAGS := -lreduction -lz -lff -lreal  $(LFLAGS)

LIBDIR := -L$(PWD)/oldf77/$(OBJSDIR) -L$(PWD)/Reduction-1.1/$(OBJSDIR) -L$(PWD)/ff/$(OBJSDIR) $(LIBDIR) 

VPATH = $(OBJSDIR)

OBJECTS =gggtt2.o qqttg2.o qgttq2.o gqbttqb2.o gggtt2virtfin.o		\
	qqttg2virtfin.o qgttq2virtfin.o gqbttqb2virtfin.o EvalLOf.o	\
	MasslessFermionLoops_gg.o MasslessFermionLoops_qq.o		\
	TopLoops_gg.o TopLoops_qq.o pentagon_gg.o pentagon_qq.o		\
	trapfpe.o RemainingDiags_gg.o RemainingDiags_qq.o		\
	masscounterterm_gg.o masscounterterm_qq.o Ghosts_gg.o		\
	Ghosts_qq.o EvalSMExx.o qqgtthel.o gggtthel.o Color.o		\
	Subtractions.o SysInfo.o wrap_virtual.o OLP_interface.o         \
        Histogram.o


ALLLIBS= $(LIBOLDF77)  $(LIBFF) $(LIBREDUCTION) 
LIBOLDF77=$(PWD)/oldf77/$(OBJSDIR)/libreal.a
LIBREDUCTION=$(PWD)/Reduction-1.1/$(OBJSDIR)/libreduction.a
LIBFF=$(PWD)/ff/$(OBJSDIR)/libff.a



library: $(ALLLIBS)
	$(MAKE) libvirtual.so.1.0.0


force_look :
	true

$(LIBOLDF77): force_look
	@echo looking into $(PWD)/oldf77/
	$(MAKE)  --directory=$(PWD)/oldf77/ COMPILER=$(COMPILER) DEBUG='$(DEBUG)'

$(LIBFF): force_look
	@echo  looking into $(PWD)/ff/
	$(MAKE)  --directory=$(PWD)/ff/ COMPILER=$(COMPILER) DEBUG='$(DEBUG)'

$(LIBREDUCTION):  force_look
	@echo looking into $(PWD)/Reduction-1.1/
	$(MAKE)  --directory=$(PWD)/Reduction-1.1/ COMPILER=$(COMPILER) DEBUG='$(DEBUG)'


#
# Write down explicit rules to switch off optimization (may cause problems)
#


EvalSMExx.o: EvalSMExx.cpp
	$(CXX) -O0 $(CFLAGS) $(IFLAGS) $(SETCXXOUT) $(OBJSDIR)/$@ $<

masscounterterm_gg.o: masscounterterm_gg.cpp
	$(CXX) -O0 $(CFLAGS) $(IFLAGS) $(SETCXXOUT) $(OBJSDIR)/$@ $<

RemainingDiags_gg.o: RemainingDiags_gg.cpp
	$(CXX) $(CFLAGS) -O0 $(IFLAGS) $(SETCXXOUT) $(OBJSDIR)/$@ $<

RemainingDiags_qq.o: RemainingDiags_qq.cpp
	$(CXX) $(CFLAGS) -O0 $(IFLAGS) $(SETCXXOUT) $(OBJSDIR)/$@ $<

Ghosts_gg.o: Ghosts_gg.cpp
	$(CXX) -O0 $(CFLAGS) $(IFLAGS) $(SETCXXOUT) $(OBJSDIR)/$@ $<

MasslessFermionLoops_gg.o: MasslessFermionLoops_gg.cpp
	$(CXX) -O0 $(CFLAGS) $(IFLAGS) $(SETCXXOUT) $(OBJSDIR)/$@ $<

TopLoops_gg.o: TopLoops_gg.cpp
	$(CXX) -O0 $(CFLAGS) $(IFLAGS) $(SETCXXOUT) $(OBJSDIR)/$@ $<

diag231.o: diag231.cpp
	$(CXX) -O0 $(CFLAGS) $(IFLAGS) $(SETCXXOUT) $(OBJSDIR)/$@ $<

diag343.o: diag343.cpp
	$(CXX) -O0 $(CFLAGS) $(IFLAGS) $(SETCXXOUT) $(OBJSDIR)/$@ $<

pentagon_gg.o: pentagon_gg.cpp
	$(CXX) $(CFLAGS) -O0 $(IFLAGS) $(SETCXXOUT) $(OBJSDIR)/$@ $<

A345div.o: A345div.cpp
	$(CXX) -O0 $(CFLAGS) $(IFLAGS) $(SETCXXOUT) $(OBJSDIR)/$@ $<


qqttg.o: qqttg.F
	$(FC) $(FFLAGS) qqttg.F $(SETCXXOUT) $(OBJSDIR)/$@ $<

%.o: %.cpp %.h
	$(CXX) $(CFLAGS) $(IFLAGS) $(SETCXXOUT) $(OBJSDIR)/$@ $<

%.o: %.cpp
	$(CXX) $(CFLAGS) $(IFLAGS) $(SETCXXOUT) $(OBJSDIR)/$@ $<

%.o: %.c
	$(CC) $(CFLAGS) $(IFLAGS) $(SETCXXOUT) $(OBJSDIR)/$@ $<

%.o: %.F
	$(FC) $(FFLAGS) $(IFLAGS) $(SETCXXOUT) $(OBJSDIR)/$@ $<

clean:
	cd $(OBJSDIR); $(RM) $(RMOPTS) $(OBJECTS) $(FASTJET_WRAP) *.o

dist-clean: clean libraryclean
	$(RM) $(RMOPTS) obj-*/*.o *.exe *~
	for subdir in $(dir $(ALLLIBS))  ; do \
	  cd $$subdir/../; $(MAKE) COMPILER=$(COMPILER) $@ ; \
	done


LVFLAGS := -Wl,-soname,libvirtual.so.1.0.0 -Wl,--version-script,version_script $(LFLAGS) 


libvirtual.so.1.0.0: $(OBJECTS) $(ALLLIBS)
	$(CXX) $(CFLAGS) -shared -o  $(OBJSDIR)/$@ $(IFLAGS) \
	$(patsubst %.o,$(OBJSDIR)/%.o,$(OBJECTS)) $(LIBDIR) $(LVFLAGS) 

libraryclean: clean
	cd $(OBJSDIR); \
	$(RM) $(RMOPTS) $(ALLLIBS) libvirtual.so.1.0.0


