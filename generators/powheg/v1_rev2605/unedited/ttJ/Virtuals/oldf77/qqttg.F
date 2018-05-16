C$Modified: Wed Sep 0 13:20:47 2006: Sun May 31 17:43:24 2009 by uwer  by puwer $
      function ampqqttg(k1,k2,k3,q1,q2,r1,r2)
      implicit none
      double precision ampqqttg
      double precision k1(4),k2(4),k3(4),q1(4),q2(4),r1(4),r2(4)
      double precision kq(4),kqb(4)
      double precision mq2,dotp
      integer mu
      double precision res,rnc
      double complex spa,spb
      double complex a1mm,a2mm,a3mm,a4mm,a1mp,a2mp,a3mp,a4mp,
     #     a1pm,a2pm,a3pm,a4pm,a1pp,a2pp,a3pp,a4pp
      parameter(rnc=3.d0)
#include "ampl-qqttg-dec.F"
      do mu=1,4
         kq(mu) = q1(mu) + q2(mu)
         kqb(mu) = r1(mu) + r2(mu)
      enddo
      mq2 = dotp(kq,kq)
#include "ampl-qqttg.F"
      a4mm = - a1mm - a2mm - a3mm
      a4mp = - a1mp - a2mp - a3mp
      a4pm = - a1pm - a2pm - a3pm
      a4pp = - a1pp - a2pp - a3pp      
#ifdef DEBUG
      write(*,*) '----------------------------'
      write(*,*) 'C1 ', a2mm*conjg(a2mm)
      write(*,*) 'C2 ', a4mm*conjg(a4mm)
      write(*,*) 'C3 ', a1mm*conjg(a1mm)
      write(*,*) 'C4 ', a3mm*conjg(a3mm)
      write(*,*) 'C1 ', a2pm*conjg(a2pm)
      write(*,*) 'C2 ', a4pm*conjg(a4pm)
      write(*,*) 'C3 ', a1pm*conjg(a1pm)
      write(*,*) 'C4 ', a3pm*conjg(a3pm)
      write(*,*) 'C1 ', a2mp*conjg(a2mp)
      write(*,*) 'C2 ', a4mp*conjg(a4mp)
      write(*,*) 'C3 ', a1mp*conjg(a1mp)
      write(*,*) 'C4 ', a3mp*conjg(a3mp)
      write(*,*) 'C1 ', a2pp*conjg(a2pp)
      write(*,*) 'C2 ', a4pp*conjg(a4pp)
      write(*,*) 'C3 ', a1pp*conjg(a1pp)
      write(*,*) 'C4 ', a3pp*conjg(a3pp)
#endif      
      res = 
     #     + rnc * ( a3mm * conjg(a3mm) + a4mm * conjg(a4mm) ) 
     #     - 1.d0/rnc * ( a1mm * conjg( a1mm ) + a2mm * conjg( a2mm ) 
     #     + 2.d0 * a1mm * conjg(a2mm) + 2.d0 * a2mm * conjg(a1mm))
     #     + rnc * ( a3mp * conjg(a3mp) + a4mp * conjg(a4mp) ) 
     #     - 1.d0/rnc * ( a1mp * conjg( a1mp ) + a2mp * conjg( a2mp ) 
     #     + 2.d0 * a1mp * conjg(a2mp) + 2.d0 * a2mp * conjg(a1mp))
     #     + rnc * ( a3pm * conjg(a3pm) + a4pm * conjg(a4pm) ) 
     #     - 1.d0/rnc * ( a1pm * conjg( a1pm ) + a2pm * conjg( a2pm ) 
     #     + 2.d0 * a1pm * conjg(a2pm) + 2.d0 * a2pm * conjg(a1pm))
     #     + rnc * ( a3pp * conjg(a3pp) + a4pp * conjg(a4pp) ) 
     #     - 1.d0/rnc * ( a1pp * conjg( a1pp ) + a2pp * conjg( a2pp ) 
     #     + 2.d0 * a1pp * conjg(a2pp) + 2.d0 * a2pp * conjg(a1pp))
      ampqqttg = res * ( rnc**2 - 1.d0 )
      end
#include "spinor.F"
      subroutine evalqr(kq,kqb,q1,q2,r1,r2)
      implicit none
      double precision kq(4),kqb(4),sq(4),sqb(4),sqrest(4),sqbrest(4)
      double precision q1(4),q2(4),r1(4),r2(4)
      double precision mq,betaq,betaqb;
      double precision dotp
      integer j,i,err
      mq = dsqrt(dotp(kq,kq))

      betaq = dsqrt(kq(1)*kq(1) + kq(2)*kq(2) 
     . + kq(3)*kq(3))/kq(4)

      betaqb = dsqrt(kqb(1)*kqb(1) + kqb(2)*kqb(2) 
     .     + kqb(3)*kqb(3))/kqb(4)

      sqrest(4) = 0.d0
      sqbrest(4) = 0.d0
  
      do i=1,3
         sqrest(i)  = kq(i)/betaq/kq(4)
         sqbrest(i) = kqb(i)/betaqb/kqb(4)
      enddo

      call boost(kq,sqrest,sq,err)
      call boost(kqb,sqbrest,sqb,err) 

      do j = 1,4
         q1(j) = ( kq(j) + mq * sq(j) ) / 2.d0
         q2(j) = ( kq(j) - mq * sq(j) ) / 2.d0
C     ------ note: the relation is reversed for the anti-quark:
         r1(j) = ( kqb(j) - mq * sqb(j) ) / 2.d0
         r2(j) = ( kqb(j) + mq * sqb(j) ) / 2.d0
      enddo
      end
