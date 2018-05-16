# DEFINE FORTRAN COMPILATION DEFAULTS

ifeq ($(shell uname),AIX)
 FFF = xlf -O4
 FF90 = xlf90 -O4 -qsuffix=f=f90:o=o90
endif

ifeq ($(shell uname),Linux)
  FFF = mpif77 -fno-automatic
  FFL = mpixlf77
  FF90 = mpif90 -fno-automatic
#FFF = g77 -O1 -Wall -fno-automatic -Wno-globals -fno-backslash \
#	 -ffast-math
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


