OBJSDIR = obj-$(COMPILER)

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

VPATH = $(OBJSDIR)/

OBJECTS = Reduction.o ScalarInt.o GGReduction.o TopologyData.o \
	ScalarInt.o FourMomentum.o svdcmp.o Matrix.o \
	StandardModelParameters.o li.o

LIB=libreduction.a

all: $(LIB)

$(LIB): $(OBJECTS)
	cd $(OBJSDIR); 	\
	ar crvs libreduction.a $(OBJECTS); \
	ranlib libreduction.a

# explicit target added to avoid ifort optimization, which 
# generate an error
Reduction.o: Reduction.cpp Reduction.h
	$(CXX) -O0 $(CFLAGS) $(IFLAGS) $(SETCXXOUT) $(OBJSDIR)/$@ $<

%.o: %.cpp %.h
	$(CXX) $(CFLAGS) $(IFLAGS) $(SETCXXOUT) $(OBJSDIR)/$@ $< 

%.o: %.cpp
	$(CXX) $(CFLAGS) $(IFLAGS) $(SETCXXOUT) $(OBJSDIR)/$@ $<

%.o: %.c
	$(CXX) $(CFLAGS) $(IFLAGS) $(SETCXXOUT) $(OBJSDIR)/$@ $<

%.o: %.F
	$(FC)  $(FFLAGS) $(SETCXXOUT) $(OBJSDIR)/$@ $<

%.o: %.f
	$(FC)  $(FFLAGS) $(SETCXXOUT) $(OBJSDIR)/$@ $<

clean:
	cd $(OBJSDIR); 	\
	$(RM) $(RMOPTS) $(OBJECTS) test.exe test.o a.out  \
	GGRedTest.o GGRedTest.exe GGTest.exe

dist-clean: clean	
	cd $(OBJSDIR); 	\
	$(RM) $(RMOPTS) libreduction.a




