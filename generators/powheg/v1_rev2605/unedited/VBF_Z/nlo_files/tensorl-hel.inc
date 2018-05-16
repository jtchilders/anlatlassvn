c  precalculated lepton tensors for Z --> l lbar 
c
c  lepton spinors and Z polarization vectors

      double complex lep(6,-1:1),lem(6,-1:1)
      common /lspinor/ lep,lem
      
      integer ie
      
      double complex ze,ae
      double precision qe

      common /zemspinor/ ze(6,-1:1), ae(6,-1:1),  qe(0:4)

c  leptonic tensors
      double complex aaee,azee,zaee,zzee,CCee, CCee6
      common /tensz/  aaee(0:3,0:3,3,-1:1),azee(0:3,0:3,3,-1:1),
     3		      zaee(0:3,0:3,3,-1:1),zzee(0:3,0:3,3,-1:1),
     3  	      CCee(0:3,0:3,3,-1:1), 
     3  	      CCee6(0:3,0:3,3,-1:1)