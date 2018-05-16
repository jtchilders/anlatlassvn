# DEFINE FORTRAN COMPILATION DEFAULTS

ifeq ($(shell uname),AIX)
 FFF = xlf -O4
 FF90 = xlf90 -O4 -qsuffix=f=f90:o=o90
endif

ifeq ($(shell uname),Linux)
  ifeq ($(shell uname -m),ppc64)
    FFF = mpixlf77 -fno-automatic -O5 -qarch=qp -qtune=qp -qsave -qsimd=auto -qhot=level=1 -qprefetch -qunroll=yes -qnoipa -qfloat=rsqrt:hssngl:fltint
    FFFL = mpixlf77 -O5 -qarch=qp -qtune=qp -qsave -qsimd=auto -qhot=level=1 -qprefetch -qunroll=yes -qnoipa -qfloat=rsqrt:hssngl:fltint
    FF90 = mpixlf90 -fno-automatic -O5 -qarch=qp -qtune=qp -qsave -qsimd=auto -qhot=level=1 -qprefetch -qunroll=yes -qnoipa -qfloat=rsqrt:hssngl:fltint
    FFL90 = mpixlf90 -O5 -qarch=qp -qtune=qp -qsave -qsimd=auto -qhot=level=1 -qprefetch -qunroll=yes -qnoipa -qfloat=rsqrt:hssngl:fltint
	 BGQ_CPP_DEFINE=-WF,# must prefix -D definitions with this for BGQ to pass it to the preprocessor
	 BGQ_SUFFIX=-qsuffix=f=f:cpp=f# treat suffix f as F so preporcess is run on BGQ
	 #FFF = mpif77 -cpp -O2
    #FF90 = mpif90 -cpp -O2
  endif
  ifeq ($(shell uname -m),x86_64)
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


