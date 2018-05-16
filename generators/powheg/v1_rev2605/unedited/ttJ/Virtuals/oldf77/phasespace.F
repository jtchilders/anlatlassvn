C-----------------------------------------------------------------------
C     $Modified: Mon Aug  9 16:08:00 2004 by puwer $
C     as intgration variables we use xg, xq the normalized energies
C     and the angles
C-----------------------------------------------------------------------
      subroutine evgggtt1(ka,kb,kq,kqb,kg,jacobian,ecms,x)
      implicit none
      double precision ka(4),kb(4),kq(4),kqh(4),kqb(4),kg(4)
      double precision x(5)
      double precision jacobian,ecms
      double precision Eg,Eq,xq,xg,zg,zq,phiq,phig,betaq
      double precision factor,xgp,xgm
      integer ipower,mu
#include "common.F"
#include "constants.F"
      xgp = (2.d0 - 4.d0 * mt/ecms)
      xgm = 0.d0
C
C     Using a ktmin we could improve the performance of the integration
C     by using:
C     xgm = 2.d0 * ktmin/ecms
C
      jacobian = 1.d0

      ipower = 2
      xg = xgm + (xgp-xgm) * x(1)**(ipower)
      jacobian = ipower * (xgp-xgm) * x(1)**(ipower-1)

      xq = (1.d0 - 2.d0 * mt/ecms) * x(2) + 2.d0 * mt/ecms
      jacobian = jacobian * (1.d0 - 2.d0 * mt/ecms)

      Eg = Ecms*xg/2.d0
      Eq = Ecms*xq/2.d0

      zg = 2.d0*x(3) - 1.d0
      jacobian = jacobian * 2.d0

      phiq = 2.d0 * pi * x(4)
      jacobian = jacobian * 2.d0 * pi

      phig = 2.d0 * pi * x(5)
      jacobian = jacobian * 2.d0 * pi

      betaq = dsqrt(1.d0 - mt**2/Eq**2)       

      zq = (2.d0 - 2.d0*xq - 2.d0 * xg + xg* xq )
     &     / ( xg * xq * betaq)

      if ( ( zq .gt. 1.d0) .or. (zq.lt. -1.d0) )  then
         jacobian = 0.d0
         return
      endif
 


      ka(1) = 0.d0
      ka(2) = ecms/2.d0
      ka(3) = 0.d0
      ka(4) = ecms/2.d0

      kb(1) = 0.d0
      kb(2) = -ecms/2.d0
      kb(3) = 0.d0
      kb(4) = ecms/2.d0

      kg(4) = Eg
      kg(1) = Eg * dcos(phig)* dsqrt(1.d0 - zg**2 )
c     ----- note that the 2 and 3 direction are exchanged
      kg(3) = Eg * dsin(phig)* dsqrt(1.d0 - zg**2 )
      kg(2) = Eg * zg

C     
C     The angles zq, phiq describe the orientation of kq
C     with respect to k3!
C     to reconstruct kq go first in the k3 co-ordinate 
      kq(4)  = Eq
      kqh(4) = Eq
      kqh(1) = kq(4) * betaq * dcos(phiq) * dsqrt(1.d0 - zq**2 )
      kqh(2) = kq(4) * betaq * dsin(phiq) * dsqrt(1.d0 - zq**2 )
      kqh(3) = kq(4) * betaq * zq
C     perform rotation and exchange 2<->3:
      kq(1) = kqh(1) * dcos(phig)*zg
     &     -  kqh(2) * dsin(phig)
     &     +  kqh(3) * dcos(phig)*dsqrt(1.d0 - zg**2)
      kq(3) = kqh(1) * dsin(phig)*zg
     &     +  kqh(2) * dcos(phig)
     &     +  kqh(3) * dsin(phig)*dsqrt(1.d0-zg**2)
      kq(2) = -kqh(1) * dsqrt(1.d0 - zg**2)
     &     +  kqh(3) * zg

      do mu=1,4
         kqb(mu) = ka(mu) + kb(mu) - kg(mu) - kq(mu)
      enddo
      
      
      end
