C
C     $Modified: Thu Aug  5 10:09:50 2004 by puwer $
C     Calculate the unpolarized matrix element squared for 
C     g(k1)g(k2) --> t(kq) tb(kqb) g(k3) g(k4)
C     All the momenta are ingoing
C
      function Mggggtt2(pset)
      implicit none
      double precision Mggggtt2
      double precision pset(4,-2:4)
      double precision k1(4),k2(4),k3(4),k4(4),kq(4),kqb(4)
      double precision q1(4),q2(4),r1(4),r2(4)
      double precision betaq,betaqb,sqrest(4),sqbrest(4),sq(4),sqb(4)
      double precision k3c(4),k4c(4),mq
      double precision dotp,ampsq
      integer mu,j,err

      do mu = 1,4
         k1(mu) = pset(mu,-1)
         k2(mu) = pset(mu,-2)
         kq(mu) = pset(mu,1)
         kqb(mu) = pset(mu,2)
         k3c(mu) = -pset(mu,3)
         k4c(mu) = -pset(mu,4)
      enddo
      mq = dsqrt(dotp(kq,kq))

C     Calculate the top spins for the helicity basis
      betaq = dsqrt(kq(1)**2 + kq(2)**2 + kq(3)**2)/kq(4)
      betaqb = dsqrt(kqb(1)**2 + kqb(2)**2 + kqb(3)**2)/kqb(4)
      sqrest(4) = 0.d0
      sqbrest(4) = 0.d0
      do j=1,3
         sqrest(j)  = kq(j)/betaq/kq(4)
         sqbrest(j) = kqb(j)/betaqb/kqb(4)
      enddo

C     ------ calculate the spins in the cms: ---------------
      call boost(kq,sqrest,sq,err)
      call boost(kqb,sqbrest,sqb,err) 

      
      call evalqr(kq,sq,kqb,sqb,q1,q2,r1,r2)

      Mggggtt2  = 2.d0 * ( ampsq(k1,k2,k3c,k4c,q1,q2,r1,r2)
     #     + ampsq(k1,k2,k3c,k4c,q2,q1,r1,r2) )

      end
