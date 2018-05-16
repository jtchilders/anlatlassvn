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

VPATH = $(OBJSDIR)/

LIBRARY	      = libff.a

AAOBJS	      = aaxbx.o \
		aaxcx.o \
		aaxdx.o \
		aaxex.o \
		aaxinv.o \
		aacbc.o \
		aaccc.o \
		aacinv.o \
		spence.o \
		npoin.o

OBJS	      = ff2dl2.o \
		ffabcd.o \
		ffca0.o \
		ffcb0.o \
		ffcb1.o \
		ffcb2.o \
		ffcb2p.o \
		ffcdb0.o \
		ffcc0.o \
		ffcc0p.o \
		ffcc1.o \
		ffcel2.o \
		ffcel3.o \
		ffcel4.o \
		ffcel5.o \
		ffceta.o \
		ffcli2.o \
		ffcrr.o \
		ffcxr.o \
		ffcxs3.o \
		ffcxs4.o \
		ffcxyz.o \
		ffdcc0.o \
		ffdcxs.o \
		ffdel2.o \
		ffdel3.o \
		ffdel4.o \
		ffdel5.o \
		ffdel6.o \
		ffdl2i.o \
		ffdl5p.o \
		ffdxc0.o \
		ffinit.o \
		ffrcvr.o \
		fftran.o \
		ffxb0.o \
		ffxb1.o \
		ffxb2p.o \
		ffxc0.o \
		ffxc0i.o \
		ffxc0p.o \
		ffxc1.o \
		ffxd0.o \
		ffxd0h.o \
		ffxd0i.o \
		ffxd0p.o \
		ffxd1.o \
		ffxdb0.o \
		ffxdbd.o \
		ffxdi.o \
		ffxdpv.o \
		ffxe0.o \
		ffxe1.o \
		ffxf0.o \
		ffxf0h.o \
		ffxli2.o \
		ffxxyz.o

# older software can be made to run with this dummy routine
#		ffini.o

AASRCS	      = aaxbx.f \
		aaxcx.f \
		aaxdx.f \
		aaxex.f \
		aaxinv.f \
		aacbc.f \
		aaccc.f \
		aacinv.f \
		spence.f \
		npoin.f

SRCS	      = ff2dl2.f \
		ffabcd.f \
		ffca0.f \
		ffcb0.f \
		ffcb1.f \
		ffcb2.f \
		ffcb2p.f \
		ffcdb0.f \
		ffcc0.f \
		ffcc0p.f \
		ffcc1.f \
		ffcdb0.f \
		ffcel2.f \
		ffcel3.f \
		ffcel4.f \
		ffcel5.f \
		ffceta.f \
		ffcli2.f \
		ffcrr.f \
		ffcxr.f \
		ffcxs3.f \
		ffcxs4.f \
		ffcxyz.f \
		ffdcc0.f \
		ffdcxs.f \
		ffdel2.f \
		ffdel3.f \
		ffdel4.f \
		ffdel5.f \
		ffdel6.f \
		ffdl2i.f \
		ffdl5p.f \
		ffdxc0.f \
		ffinit.f \
		ffrcvr.f \
		fftran.f \
		ffxb0.f \
		ffxb1.f \
		ffxb2p.f \
		ffxc0.f \
		ffxc0i.f \
		ffxc0p.f \
		ffxc1.f \
		ffxd0.f \
		ffxd0h.f \
		ffxd0i.f \
		ffxd0p.f \
		ffxd1.f \
		ffxdb0.f \
		ffxdbd.f \
		ffxdi.f \
		ffxdpv.f \
		ffxe0.f \
		ffxe1.f \
		ffxf0.f \
		ffxf0h.f \
		ffxli2.f \
		ffxxyz.f

all: $(LIBRARY)

$(LIBRARY): $(OBJS) $(AAOBJS)
	cd $(OBJSDIR); ar cru $(LIBRARY) $(OBJS) $(AAOBJS); ranlib $(LIBRARY)

ffinit.o: ffinit.f
	@sh set_ff_path.sh 
	$(FC) $(FFLAGS) $(IFLAGS) $(SETCXXOUT) $(OBJSDIR)/$@ $<

%.o: %.f
	$(FC) $(FFLAGS) $(IFLAGS) $(SETCXXOUT) $(OBJSDIR)/$@ $<

clean:	
	cd $(OBJSDIR); \
	rm -f $(OBJS) $(AAOBJS) npointes.o npointes

dist-clean: clean
	cd $(OBJSDIR); \
	rm -f $(LIBRARY) 
