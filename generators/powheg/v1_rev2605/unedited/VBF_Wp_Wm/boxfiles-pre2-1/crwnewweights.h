c -*- fortran -*-
      integer nscmax,nscales
      parameter (nscmax=10)
      real * 8 rwnewweight,rwrenfact,rwfacfact,rwnpdf1,rwnpdf2
      common/crwnewweights/rwnewweight(nscmax),rwrenfact(nscmax),
     1         rwfacfact(nscmax),rwnpdf1(nscmax),rwnpdf2(nscmax),nscales
