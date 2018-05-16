
OBJSDIR = obj-$(COMPILER)
ifeq ("$(COMPILER)", "intel")
  CXX = icpc
  CC = icc
  FC = ifort
  SETCXXOUT = -c -o
  SETLDOUT = -o
  LFLAGS = -lifcore -lirc -limf  -lpthread  -lifport -lifcoremt
  FFLAGS := $(FFLAGS) -132 -fPIC  	
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
IFLAGS := -I./auto -I../auto $(IFLAGS)

OBJECTS =  evalqr.o spinor.o ggggtt.o gggtt.o color.o Dggg.o Dqgq.o Soft.o \
   subtractions.o SoftConf.o CollConf.o CheckKinematics.o Mgggtt2.o \
   Mggggtt2.o umgggtt2.o observables.o durham.o ampsquare.o  \
   gggtthel.o ggggtthel.o ColorCorr.o phasespace.o Convolution.o \
   event4.o eventn.o dggggtt2.o Tools.o

VPATH = $(OBJSDIR)/

LIBRARY = libreal.a

all: $(LIBRARY)

$(LIBRARY): $(OBJECTS) 
	cd $(OBJSDIR); 	ar crvs libreal.a $(OBJECTS); ranlib libreal.a 

color.o: SetupColor.F color.f; 
	$(FC) $(CFLAGS) $(DEBUG) -c -o $(OBJSDIR)/color.o $<

EventGen.a: spinor.o eventn.o;
	cd $(OBJSDIR); 	\
	ar crvs EventGen.a spinor.o eventn.o; \
	ranlib EventGen.a           
clean:
	cd $(OBJSDIR); \
	$(RM) $(RMOPTS) $(OBJECTS) 

dist-clean: clean
	cd $(OBJSDIR); $(RM) $(RMOPTS) libreal.a CalcObs++ CalcObs++.o
#Dggg.dec Dggg.F Soft.dec Soft.F Dqgq.dec Dqgq.F 




%.o: %.cpp %.h
	$(CXX) $(CFLAGS) $(IFLAGS) $(SETCXXOUT) $(OBJSDIR)/$@ $<

%.o: %.cpp
	$(CXX) $(CFLAGS) $(IFLAGS) $(SETCXXOUT) $(OBJSDIR)/$@ $<

%.o: %.c
	$(CC) $(CFLAGS) $(IFLAGS) $(SETCXXOUT) $(OBJSDIR)/$@ $<

%.o: %.F %.dec
	$(FC)  $(FFLAGS) $(IFLAGS)  $(SETCXXOUT) $(OBJSDIR)/$@ $<
%.o: %.F 
	$(FC)  $(FFLAGS) $(IFLAGS)  $(SETCXXOUT) $(OBJSDIR)/$@ $<


