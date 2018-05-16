# DEFINE FORTRAN COMPILATION DEFAULTS

ifeq ($(shell uname),AIX)
 FFF = xlf -O4
 FF90 = xlf90 -O4 -qsuffix=f=f90:o=o90
endif


ifeq ($(shell uname),Linux)
  ifneq (,$(findstring edison,$(shell hostname)))
    FFF = ftn -fno-automatic -cpp -O2 -g
    FFFL = ${FFF}
    FF90 = ftn -fno-automatic -cpp -O2 -g #-Wall -Wshadow -Wcast-align
    FFL90 = ${FF90}
    BGQ_CPP_DEFINE=
    BGQ_SUFFIX=
  else ifeq ($(shell uname -m),ppc64)
    COMMON_OPTS=-O5 -qmaxmem=100000 -qarch=qp -qtune=qp -qsave -qstrict -qsimd=auto -qhot=level=1 -qprefetch -qunroll=yes -qnoipa -qfloat=rsqrt:hssngl:fltint
    COPTS=-fno-automatic $(COMMON_OPTS)
    LOPTS=$(COMMON_OPTS)
    FFF = mpixlf77 $(COPTS)
    FFFL = mpixlf77 $(LOPTS)
    FF90 = mpixlf90 $(COPTS)
    FFL90 = mpixlf90 $(LOPTS)
    BGQ_CPP_DEFINE=-WF,# must prefix -D definitions with this for BGQ to pass it to the preprocessor
    BGQ_SUFFIX=-qsuffix=f=f:cpp=f# treat suffix f as F so preporcess is run on BGQ
    #FFF = mpif77 -cpp -O2
    #FF90 = mpif90 -cpp -O2
  else ifeq ($(shell uname -m),x86_64)
    FFF = mpif77 -fno-automatic -cpp -O2 -g -fbacktrace #-Wall -Wshadow -Wcast-align
    FFFL = ${FFF}
    FF90 = mpif90 -fno-automatic -cpp -O2 -g -fbacktrace #-Wall -Wshadow -Wcast-align
    FFL90 = ${FF90}
    BGQ_CPP_DEFINE=
    BGQ_SUFFIX=
  endif
endif

ifeq ($(shell uname),Darwin)
  ifeq  ($(shell uname -m),Power Macintosh) 
     FFF = g77 -O1 -Wall -fno-automatic -Wno-globals -fno-backslash \
	 -ffast-math
  endif
  ifeq  ($(shell uname -m),i386) 
     FFF = ifort -V -noautomatic -warn
     FF90 = ifort -V -noautomatic -warn
     FFF = gfortran -fno-automatic
	FF90 = gfortran -fno-automatic
  endif
endif
# Unix-ALPHA fortran
# FFF = f90 -fast
# FF90 = f90 -fast  


