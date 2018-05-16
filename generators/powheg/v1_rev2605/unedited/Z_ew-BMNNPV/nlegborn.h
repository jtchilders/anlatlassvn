c -*- Fortran -*-

c The user must set nlegborn to the appropriate value for his process.
      integer nlegborn,nlegreal
      
      parameter (nlegborn=4)
      parameter (nlegreal=nlegborn+1)

c     ndiminteg is the dimensionality of the full real integral
c     ndiminteg=(nlegreal-2)*3-4+2-1
c     if there are undecayed resonances, we need extra variables to pilot
c     the resonance's masses

      integer ndiminteg
      parameter (ndiminteg=(nlegreal-2)*3-4+2-1+1)

c     the last +1 is for double channel Z/gamma
 

      integer maxprocborn,maxprocreal
      parameter (maxprocborn=200,maxprocreal=500)


	  integer nparton
	  parameter (nparton=22)


	  logical ifdis
	  common/ifdiscommon/ifdis
