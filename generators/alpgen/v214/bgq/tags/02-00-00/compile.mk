# DEFINE FORTRAN COMPILATION DEFAULTS

ifeq ($(shell uname),AIX)
 FFF = xlf -O4
 FF90 = xlf90 -O4 -qsuffix=f=f90:o=o90
endif

ifeq ($(shell uname),Linux)
  ifeq ($(shell uname -m),ppc64)
    FFF = mpixlf77 -O5 -qarch=qp -qtune=qp -qsave -qsimd=auto -qhot=level=1 -qprefetch -qunroll=yes -qnoipa -qfloat=rsqrt:hssngl:fltint
    FF90 = mpixlf90 -O5 -qarch=qp -qtune=qp -qsave -qsimd=auto -qhot=level=1 -qprefetch -qunroll=yes -qnoipa -qfloat=rsqrt:hssngl:fltint
  endif
  ifeq ($(shell uname -m),x86_64)
    FFF = mpif77 -fno-automatic 
    FF90 = mpif90 -fno-automatic
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


