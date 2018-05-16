      function born_bmunu(mu,nu,k1,k2,k3,k4,k5)
      implicit none
      include "PhysPars.h"
      integer mu,nu
      double precision k1(0:3),k2(0:3),k3(0:3),k4(0:3),k5(0:3)
      double precision k1k2,k1k3,k1k4,k1k5,k2k3,k2k4,k2k5,k3k4,k3k5,k4k5
      double precision PropW2,PropB2,PropT2,PropBT


      double precision born_bmunu
      double complex bmunutmp

      double precision Pair,Eps,EpsComp,dotp
      external Pair,Eps,EpsComp,dotp

      double precision metric(0:3,0:3)
      data metric/1d0, 0d0, 0d0, 0d0,
     #            0d0,-1d0, 0d0, 0d0,
     #            0d0, 0d0,-1d0, 0d0,
     #            0d0, 0d0, 0d0,-1d0/
      save metric

      
      k1k2=dotp(k1,k2)
      k1k3=dotp(k1,k3)
      k1k4=dotp(k1,k4)
      k1k5=dotp(k1,k5)
      k2k3=dotp(k2,k3)
      k2k4=dotp(k2,k4)
      k2k5=dotp(k2,k5)
      k3k4=dotp(k3,k4)
      k3k5=dotp(k3,k5)
      k4k5=dotp(k4,k5)


      PropW2=1./((-2.*k1k5-ph_Wmass**2)**2 + ph_Wmass**2 * ph_Wwidth**2)
      PropT2=1./((-2.*k2k3)**2 + topmass_pow**2 * topwidth_pow**2)
      PropB2=1./(-2.*k2k4)**2
      PropBT=1./(-2.*k2k3)/(-2.*k2k4)



      bmunutmp=
     $ PropW2*(512*k3k5*PropB2*
     -     (-(k1k2*k2k4*metric(mu,nu)) + 
     -       Pair(mu,k4)*(k2k4*Pair(nu,k1) + 
     -          (k1k2 - k1k4)*(Pair(nu,k2) - 2*Pair(nu,k4))) + 
     -       (k2k4*Pair(mu,k1) + (k1k2 - k1k4)*Pair(mu,k2))*Pair(nu,k4))
     -      + 128*PropBT*(4*k1k5*k2k3*k2k4*metric(mu,nu) - 
     -       4*k1k3*k2k4*k2k5*metric(mu,nu) + 
     -       4*k1k2*k2k5*k3k4*metric(mu,nu) - 
     -       4*k1k2*k2k3*k4k5*metric(mu,nu) + 
     -       (0,2)*k2k5*Eps(k1,k2,k3,k4)*metric(mu,nu) - 
     -       (0,1)*k2k4*EpsComp(nu,k2,k3,k5)*Pair(mu,k1) + 
     -       (0,1)*k2k3*EpsComp(nu,k2,k4,k5)*Pair(mu,k1) + 
     -       (0,1)*k2k5*EpsComp(nu,k1,k3,k4)*Pair(mu,k2) + 
     -       (0,1)*k1k2*EpsComp(nu,k3,k4,k5)*Pair(mu,k2) - 
     -       (0,1)*k2k5*EpsComp(nu,k1,k2,k4)*Pair(mu,k3) - 
     -       (0,1)*k1k2*EpsComp(nu,k2,k4,k5)*Pair(mu,k3) + 
     -       (0,1)*k2k5*EpsComp(nu,k1,k2,k3)*Pair(mu,k4) + 
     -       (0,1)*k1k2*EpsComp(nu,k2,k3,k5)*Pair(mu,k4) - 
     -       (0,1)*k2k4*EpsComp(nu,k1,k2,k3)*Pair(mu,k5) + 
     -       (0,1)*k2k3*EpsComp(nu,k1,k2,k4)*Pair(mu,k5) - 
     -       (0,1)*k2k4*EpsComp(mu,k2,k3,k5)*Pair(nu,k1) + 
     -       (0,1)*k2k3*EpsComp(mu,k2,k4,k5)*Pair(nu,k1) - 
     -       2*k2k5*k3k4*Pair(mu,k2)*Pair(nu,k1) + 
     -       2*k2k4*k3k5*Pair(mu,k2)*Pair(nu,k1) + 
     -       2*k2k3*k4k5*Pair(mu,k2)*Pair(nu,k1) + 
     -       4*k2k4*k2k5*Pair(mu,k3)*Pair(nu,k1) - 
     -       4*k2k4*k3k5*Pair(mu,k3)*Pair(nu,k1) - 
     -       4*k2k3*k2k4*Pair(mu,k5)*Pair(nu,k1) + 
     -       (0,1)*k2k5*EpsComp(mu,k1,k3,k4)*Pair(nu,k2) + 
     -       (0,1)*k1k2*EpsComp(mu,k3,k4,k5)*Pair(nu,k2) - 
     -       2*k2k5*k3k4*Pair(mu,k1)*Pair(nu,k2) + 
     -       2*k2k4*k3k5*Pair(mu,k1)*Pair(nu,k2) + 
     -       2*k2k3*k4k5*Pair(mu,k1)*Pair(nu,k2) + 
     -       4*k1k5*k3k4*Pair(mu,k2)*Pair(nu,k2) - 
     -       4*k1k4*k3k5*Pair(mu,k2)*Pair(nu,k2) - 
     -       4*k1k3*k4k5*Pair(mu,k2)*Pair(nu,k2) - 
     -       2*k1k5*k2k4*Pair(mu,k3)*Pair(nu,k2) - 
     -       2*k1k4*k2k5*Pair(mu,k3)*Pair(nu,k2) + 
     -       4*k1k4*k3k5*Pair(mu,k3)*Pair(nu,k2) + 
     -       2*k1k2*k4k5*Pair(mu,k3)*Pair(nu,k2) - 
     -       2*k1k5*k2k3*Pair(mu,k4)*Pair(nu,k2) + 
     -       2*k1k3*k2k5*Pair(mu,k4)*Pair(nu,k2) - 
     -       2*k1k2*k3k5*Pair(mu,k4)*Pair(nu,k2) + 
     -       4*k1k4*k3k5*Pair(mu,k4)*Pair(nu,k2) + 
     -       2*k1k4*k2k3*Pair(mu,k5)*Pair(nu,k2) + 
     -       2*k1k3*k2k4*Pair(mu,k5)*Pair(nu,k2) - 
     -       2*k1k2*k3k4*Pair(mu,k5)*Pair(nu,k2) - 
     -       (0,1)*Eps(k1,k2,k3,k4)*Pair(mu,k5)*Pair(nu,k2) - 
     -       (0,1)*Eps(k2,k3,k4,k5)*
     -        (2*k1k2*metric(mu,nu) - Pair(mu,k2)*Pair(nu,k1) - 
     -          Pair(mu,k1)*Pair(nu,k2)) - 
     -       (0,1)*k2k5*EpsComp(mu,k1,k2,k4)*Pair(nu,k3) - 
     -       (0,1)*k1k2*EpsComp(mu,k2,k4,k5)*Pair(nu,k3) + 
     -       4*k2k4*k2k5*Pair(mu,k1)*Pair(nu,k3) - 
     -       4*k2k4*k3k5*Pair(mu,k1)*Pair(nu,k3) - 
     -       2*k1k5*k2k4*Pair(mu,k2)*Pair(nu,k3) - 
     -       2*k1k4*k2k5*Pair(mu,k2)*Pair(nu,k3) + 
     -       4*k1k4*k3k5*Pair(mu,k2)*Pair(nu,k3) + 
     -       2*k1k2*k4k5*Pair(mu,k2)*Pair(nu,k3) - 
     -       4*k1k2*k2k5*Pair(mu,k4)*Pair(nu,k3) + 
     -       4*k1k4*k2k5*Pair(mu,k4)*Pair(nu,k3) + 
     -       4*k1k2*k3k5*Pair(mu,k4)*Pair(nu,k3) - 
     -       8*k1k4*k3k5*Pair(mu,k4)*Pair(nu,k3) + 
     -       (0,1)*k2k5*EpsComp(mu,k1,k2,k3)*Pair(nu,k4) + 
     -       (0,1)*k1k2*EpsComp(mu,k2,k3,k5)*Pair(nu,k4) - 
     -       2*k1k5*k2k3*Pair(mu,k2)*Pair(nu,k4) + 
     -       2*k1k3*k2k5*Pair(mu,k2)*Pair(nu,k4) - 
     -       2*k1k2*k3k5*Pair(mu,k2)*Pair(nu,k4) + 
     -       4*k1k4*k3k5*Pair(mu,k2)*Pair(nu,k4) - 
     -       4*k1k2*k2k5*Pair(mu,k3)*Pair(nu,k4) + 
     -       4*k1k4*k2k5*Pair(mu,k3)*Pair(nu,k4) + 
     -       4*k1k2*k3k5*Pair(mu,k3)*Pair(nu,k4) - 
     -       8*k1k4*k3k5*Pair(mu,k3)*Pair(nu,k4) + 
     -       4*k1k2*k2k3*Pair(mu,k5)*Pair(nu,k4) - 
     -       4*k1k4*k2k3*Pair(mu,k5)*Pair(nu,k4) + 
     -       ((0,-1)*k2k4*EpsComp(mu,k1,k2,k3) + 
     -          (0,1)*k2k3*EpsComp(mu,k1,k2,k4) - 
     -          4*k2k3*k2k4*Pair(mu,k1) + 
     -          (2*k1k4*k2k3 + 2*k1k3*k2k4 - 2*k1k2*k3k4 - 
     -             (0,1)*Eps(k1,k2,k3,k4))*Pair(mu,k2) + 
     -          4*(k1k2 - k1k4)*k2k3*Pair(mu,k4))*Pair(nu,k5)) - 
     -    512*k1k4*PropT2*(k2k3*k2k5*metric(mu,nu) - 
     -       ((k2k5 - k3k5)*Pair(mu,k2) + k2k3*Pair(mu,k5))*
     -        Pair(nu,k3) - Pair(mu,k3)*
     -        ((k2k5 - k3k5)*Pair(nu,k2) - 
     -          2*(k2k5 - k3k5)*Pair(nu,k3) + k2k3*Pair(nu,k5))))


      born_bmunu=dble(bmunutmp)
c      print*, bmunutmp

      end



      function Pair(index,lorentzvec)
      implicit none
      integer index
      real *8 lorentzvec(0:3)
      real *8 Pair
      
      Pair=lorentzvec(index)
      return
      end



      function Eps(k1,k2,k3,k4)
      implicit none
      double precision k1(0:3),k2(0:3),k3(0:3),k4(0:3)
      double precision p1(4),p2(4),p3(4),p4(4)
      double precision Eps,epf
      external epf
      integer mu
      do mu=1,3
         p1(mu)=k1(mu)
         p2(mu)=k2(mu)
         p3(mu)=k3(mu)
         p4(mu)=k4(mu)
      enddo
      p1(4)=k1(0)
      p2(4)=k2(0)
      p3(4)=k3(0)
      p4(4)=k4(0)

      Eps=Epf(p1,p2,p3,p4)
      return
      end



c http://wwwasd.web.cern.ch/wwwasd/cernlib/download/2002_source/src/mclibs/isajet/code/epf.F
      FUNCTION EPF(A,B,C,D)
C          CALCULATE TOTALLY ANTISYMMETRIC TENSOR EPSILON CONTRACTED
C          WITH FOUR 4-VECTORS.
      DIMENSION A(4),B(4),C(4),D(4)
      DOUBLE PRECISION EPF
      DOUBLE PRECISION A,B,C,D,CD,BCD
      CD(I,J)=C(I)*D(J)-C(J)*D(I)
      BCD(I,J,K)=B(I)*CD(J,K)-B(J)*CD(I,K)+B(K)*CD(I,J)
      EPF=A(1)*BCD(2,3,4)-A(2)*BCD(1,3,4)+A(3)*BCD(1,2,4)
     1-A(4)*BCD(1,2,3)
      RETURN
      END



      function EpsComp(mu,k1,k2,k3)
      integer mu,i
      double precision k1(0:3),k2(0:3),k3(0:3),kaux(0:3)
      double precision EpsComp,Eps
      external Eps
      do i=0,3
         kaux(i)=0.
      enddo
      kaux(mu)=1.
      EpsComp=Eps(kaux,k1,k2,k3)
      return
      end
