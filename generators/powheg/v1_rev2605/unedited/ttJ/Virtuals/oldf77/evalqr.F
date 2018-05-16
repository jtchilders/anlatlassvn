C     $Modified: Thu Apr  1 15:37:26 2004 by puwer $
C----------------------------------------------------------------------
C
C     Calculate the qi,ri from the given sq,kq,...
C
C----------------------------------------------------------------------
      subroutine evalqr(kq,sq,kqb,sqb,q1,q2,r1,r2)
      implicit none
      double precision kq(4),kqb(4),sq(4),sqb(4)
      double precision q1(4),q2(4),r1(4),r2(4)
      double precision mq
      double precision dotp
      integer j
      mq = dsqrt(dotp(kq,kq))
      do j = 1,4
         q1(j) = ( kq(j) + mq * sq(j) ) / 2.d0
         q2(j) = ( kq(j) - mq * sq(j) ) / 2.d0
C     ------ note: the relation is reversed for the anti-quark:
         r1(j) = ( kqb(j) - mq * sqb(j) ) / 2.d0
         r2(j) = ( kqb(j) + mq * sqb(j) ) / 2.d0
      enddo
      end
