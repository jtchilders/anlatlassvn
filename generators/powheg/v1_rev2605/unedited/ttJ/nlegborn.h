c -*- Fortran -*-

c The user must set nlegborn to the appropriate value for his process.
      integer nlegborn,nlegreal
      
      parameter (nlegborn=5)
      parameter (nlegreal=nlegborn+1)

c     ndiminteg is the dimensionality of the full real integral
c     ndiminteg=(nlegreal-2)*3-4+2-1
c     if there are undecayed resonances, we need extra variables to pilot
c     the resonance's masses

      integer ndiminteg
      parameter (ndiminteg=(nlegreal-2)*3-4+2-1
     1     +0  )          ! top and antitop masses kept fixed
c(FOLLOWING CHOICES ARE NOT YET SUPPORTED BY ROUTINES FOR VIRTUAL AMPLITUDES)   c     2    +1  )          ! top and antitop masses along a BW, but equal  
c     3    + 2 )          ! top and antitop are two different resonances 

 
      integer maxprocborn,maxprocreal
      parameter (maxprocborn=200,maxprocreal=500)
