c event kinematics, ordered according to generation ordering
      real*8 pt(maxpar),eta(maxpar),dr(maxpar,maxpar),
     +  dphi(maxpar,maxpar)
      common/evtparr/pt,eta,dr,dphi
c decay information
      character*1 idecay
      common/tdecflag/idecay
      real*8 idec(4,4,maxpar)
      common/decay/idec     ! idec(1:4,i,j) four momentum of the I-th 
C            decay product of the J-th particle; 
C            the four momentum is (p_x,p_y,p_z,E)
C            for the W and Z I=1 is fermion I=2 is fermionbar 
C            I=3,4 only relevant for higgs decay
C
c      integer zdmode,zfl(maxpar-2)
c      common/zfldec/zdmode,zfl
