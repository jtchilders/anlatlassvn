c -*- Fortran -*-
      
      integer nlegborn,nlegreal
      parameter (nlegborn=           5 )
      parameter (nlegreal=nlegborn+1)
      integer ndiminteg
      parameter (ndiminteg=(nlegreal-2)*3-4+2-1+1
c     one extra variable to separate ISR and FSR for second jet
     1 +1 )
      
      integer maxprocborn,maxprocreal
      parameter (maxprocborn=999,maxprocreal=999)
 
