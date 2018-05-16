C$Modified: Fri Jun  4 09:57:18 2004 by puwer $
      subroutine gggtthel(k1,k2,k3,q1,q2,r1,r2,
     .     ammm,ammp,ampm,ampp,apmm,apmp,appm,appp)
      implicit none
      double complex ammm,ammp,ampm,ampp,apmm,apmp,appm,appp
      double precision k1(4), k2(4),k3(4),q1(4),q2(4)
      double precision r1(4),r2(4),kq(4), kqb(4) 
      double precision dotp
      double complex spa,spb
      double precision m
#include "all_amp_auto-dec3g.F"
      call add(q1,q2,kq) 
      call add(r1,r2,kqb) 
      m = dsqrt(2.d0 * dotp(q1,q2))
#include "all_amp_auto3g.F"     
      end
