C$Modified: Sun May 31 12:49:08 2009 by uwer $
C----------------------------------------------------------------------
C
C     Function to calculate squared matrixelements
C     for gg --> ttbar g 
C
C----------------------------------------------------------------------
      double precision function Mgggtt2(k1,k2,k3,q1,q2,r1,r2)
C     calculate |Tfi|^2 from the color ordered subamplitudes
      IMPLICIT DOUBLE complex (A)
      double precision k1(4), k2(4),k3(4),q1(4),q2(4),r1(4),r2(4)
      double precision rnc
      parameter (rnc=3.d0)
      double precision res
      double complex t1,t2,t3,s1
C      double complex 123,132,321,312,231,213
#include "permutations3g.F"
C     number of colors
C      rnc = 3.d0
C
      Mgggtt2 = rnc**2*(rnc**2-1.d0) * 
     #     ( t1 - 1.d0/rnc**2 * t2 + (rnc**2+1.d0)/rnc**4 * t3)
      return
      end
