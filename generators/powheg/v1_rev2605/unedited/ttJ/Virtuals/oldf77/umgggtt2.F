C
C     $Modified: Wed Jul 14 10:33:35 2004 by puwer $
C     Calculate the unpolarized matrix element squared for 
C     g(k1)g(k2)->t(kq)tb(kqb)g(k3)
C     All the momenta are ingoing
C
      function UMgggtt2(pset)
      implicit none
      double precision UMgggtt2
      double precision pset(4,-2:3)
      double precision k1(4),k2(4),k3(4),kq(4),kqb(4)
      double precision q1(4),q2(4),r1(4),r2(4)
      double precision betaq,betaqb,sqrest(4),sqbrest(4),sq(4),sqb(4)
      double precision k3c(4),mq
      double precision dotp,mgggtt2
      integer mu,j,err
      do mu = 1,4
         k1(mu) = pset(mu,-1)
         k2(mu) = pset(mu,-2)
         kq(mu) = pset(mu,1)
         kqb(mu) = pset(mu,2)
         k3c(mu) = -pset(mu,3)
      enddo
      mq = dsqrt(dotp(kq,kq))
C     Calculate the top spins
      betaq = dsqrt(kq(1)**2 + kq(2)**2 + kq(3)**2)/kq(4)
      betaqb = dsqrt(kqb(1)**2 + kqb(2)**2 + kqb(3)**2)/kqb(4)
      sqrest(4) = 0
      sqbrest(4) = 0
      do j=1,3
         sqrest(j)  = kq(j)/betaq/kq(4)
         sqbrest(j) = kqb(j)/betaqb/kqb(4)
      enddo
C     ------ calculate the spins in the cms: ---------------
      call boost(kq,sqrest,sq,err)
      call boost(kqb,sqbrest,sqb,err) 
      
      call evalqr(kq,sq,kqb,sqb,q1,q2,r1,r2) 

      UMgggtt2  = 2.d0 * ( mgggtt2(k1,k2,k3c,q1,q2,r1,r2)
     #     + mgggtt2(k1,k2,k3c,q2,q1,r1,r2) )

      end
