C$Modified: Wed Dec 13 13:59:14 2006 by puwer $
      subroutine qqgtthel(k1,k2,k3,q1,q2,r1,r2,
     .     a1mm,a1mp,a1pm,a1pp,a2mm,a2mp,a2pm,a2pp,
     .     a3mm,a3mp,a3pm,a3pp,a4mm,a4mp,a4pm,a4pp)
      implicit none
      double complex a1mm,a2mm,a3mm,a4mm,a1mp,a2mp,a3mp,a4mp,
     #     a1pm,a2pm,a3pm,a4pm,a1pp,a2pp,a3pp,a4pp
      double precision k1(4), k2(4),k3(4),q1(4),q2(4),r1(4),r2(4)
      double precision kq(4),kqb(4)
      double precision mq2      !quark mass squared
      double precision dotp
      double complex spa,spb
      integer i
#include "auto/ampl-qqttg-dec.F"
C     
      do i=1,4
         kq(i)  = q1(i) + q2(i)
         kqb(i) = r1(i) + r2(i)
      enddo
      mq2 = dotp(kq,kq)
#include "auto/ampl-qqttg.F"  
      a4mm = - a1mm - a2mm - a3mm
      a4mp = - a1mp - a2mp - a3mp
      a4pm = - a1pm - a2pm - a3pm
      a4pp = - a1pp - a2pp - a3pp
      end
