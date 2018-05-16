C     $Modified: Sun May 31 12:52:53 2009 by uwer $
      function ampsq(k1,k2,k3,k4,q1,q2,r1,r2)
      implicit none
      double precision ampsq
      double precision N
      parameter (N=3.d0)
      double precision k1(4),k2(4),k3(4),k4(4),q1(4),q2(4),
     #     r1(4),r2(4)
      double complex AMMMM, AMMMP, AMMPM, AMMPP, AMPMM, AMPMP, AMPPM, 
     #     AMPPP, APMMM, APMMP, APMPM, APMPP, APPMM, APPMP, APPPM, 
     #     APPPP
      double complex sMMMM, sMMMP, sMMPM, sMMPP, sMPMM, sMPMP, sMPPM, 
     #     sMPPP, sPMMM, sPMMP, sPMPM, sPMPP, sPPMM, sPPMP, sPPPM, 
     #     sPPPP
#include "permutations-dec4g.F"
#include "permutations4g.F"
C      N = 3.d0
      ampsq = (N**2-1.d0)/N**3
     #     * (N**6 * C1 + H2 * N**4 + C3 * N**2/2.d0 - 3.d0*C4*N**2-C4)
      end
