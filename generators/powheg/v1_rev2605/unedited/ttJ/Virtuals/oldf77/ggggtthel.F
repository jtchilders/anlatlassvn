C$Modified: Sun May 31 12:59:51 2009 by uwer $
C
C    This is automatically created code, produced by ./mkcode.csh
C    on Fri Jun  4 15:10:38 CEST 2004.
C
      subroutine ggggtthel(k1,k2,k3,k4,q1,q2,r1,r2,
     . ammmm,ammmp,ammpm,ammpp,ampmm,ampmp,amppm,amppp,
     . apmmm,apmmp,apmpm,apmpp,appmm,appmp,apppm
     . ,apppp)
      implicit none
      double precision k1(4), k2(4),k3(4),k4,q1(4),q2(4)
      double precision r1(4),r2(4),kq(4), kqb(4) 
      double precision dotp,den
      double complex ammmm,ammmp,ammpm,ammpp,ampmm,ampmp,
     . amppm,amppp,
     . apmmm,apmmp,apmpm,apmpp,appmm,appmp,apppm,apppp
      double complex spa,spb
      double precision m
      double complex s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
      double complex s11,s12,s13,s14,s15,s16,s17,s18
#include "all_amp_auto-dec4g.F"
      call add(q1,q2,kq) 
      call add(r1,r2,kqb) 
      m = dsqrt(dotp(kq,kq))
#include "all_amp_auto4g.F" 
      end
