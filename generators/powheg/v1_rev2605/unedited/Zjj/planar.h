c -*- Fortran -*-
      integer ncstructmax
      parameter (ncstructmax=2)
      integer ncstruct,nclines,clineofpart(ncstructmax,nlegborn)
      double precision sigmaofncstruct(ncstructmax)
      common/planar/ncstruct,nclines,clineofpart,sigmaofncstruct
      
      
