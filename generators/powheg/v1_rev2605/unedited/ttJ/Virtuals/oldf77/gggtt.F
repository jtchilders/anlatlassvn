C
C    This is automatically created code, produced by ./mkcode.csh
C    on Mon Jan 19 17:29:16 CET 2004.
C
      function aMMM(k1,k2,k3,q1,q2,r1,r2)
      implicit none
      double complex aMMM
      double precision k1(4), k2(4),k3(4),q1(4),q2(4)
      double precision r1(4),r2(4),kq(4), kqb(4) 
      double precision dotp
      double complex spa,spb
      double precision m
#include "gggttbMMM-dec.F"
***#include "../common.F"
#ifdef DEBUG
      write(*,*) 'aMMM :'
#endif
      call add(q1,q2,kq) 
      call add(r1,r2,kqb) 
      m = dsqrt(2.d0 * dotp(q1,q2))
#ifdef DEBUG2 
      write(*,*) 'aMMM : m = ',m   
      write(*,*) 'aMMM : spa(k1,k2) = ', spa(k1,k2)  
      write(*,*) 'aMMM : spa(k1,k3) = ', spa(k1,k3)  
      write(*,*) 'aMMM : spa(k2,k3) = ', spa(k2,k3) 
      write(*,*) 'aMMM : spa(k1,q1) = ', spa(k1,q1)  
      write(*,*) 'aMMM : spa(k1,q2) = ', spa(k1,q2)  
      write(*,*) 'aMMM : spa(k2,q1) = ', spa(k2,q1)  
      write(*,*) 'aMMM : spa(k2,q2) = ', spa(k2,q2)  
      write(*,*) 'aMMM : spa(k3,q1) = ', spa(k3,q1)  
      write(*,*) 'aMMM : spa(k3,q2) = ', spa(k3,q2)  
      write(*,*) 'aMMM : spa(q1,q2) = ', spa(q1,q2)  
#endif
#include "gggttbMMM.F" 
#ifdef DEBUG 
      write(*,*) 'aMMM : return value', aMMM  
#endif 
      end
      function aMMP(k1,k2,k3,q1,q2,r1,r2)
      implicit none
      double complex aMMP
      double precision k1(4), k2(4),k3(4),q1(4),q2(4)
      double precision r1(4),r2(4),kq(4), kqb(4) 
      double precision dotp
      double complex spa,spb
      double precision m
#include "gggttbMMP-dec.F"
***#include "../common.F"
#ifdef DEBUG
      write(*,*) 'aMMP :'
#endif
      call add(q1,q2,kq) 
      call add(r1,r2,kqb) 
      m = dsqrt(2.d0 * dotp(q1,q2))
#ifdef DEBUG2 
      write(*,*) 'aMMP : m = ',m   
      write(*,*) 'aMMP : spa(k1,k2) = ', spa(k1,k2)  
      write(*,*) 'aMMP : spa(k1,k3) = ', spa(k1,k3)  
      write(*,*) 'aMMP : spa(k2,k3) = ', spa(k2,k3) 
      write(*,*) 'aMMP : spa(k1,q1) = ', spa(k1,q1)  
      write(*,*) 'aMMP : spa(k1,q2) = ', spa(k1,q2)  
      write(*,*) 'aMMP : spa(k2,q1) = ', spa(k2,q1)  
      write(*,*) 'aMMP : spa(k2,q2) = ', spa(k2,q2)  
      write(*,*) 'aMMP : spa(k3,q1) = ', spa(k3,q1)  
      write(*,*) 'aMMP : spa(k3,q2) = ', spa(k3,q2)  
      write(*,*) 'aMMP : spa(q1,q2) = ', spa(q1,q2)  
#endif
#include "gggttbMMP.F" 
#ifdef DEBUG 
      write(*,*) 'aMMP : return value', aMMP  
#endif 
      end
      function aMPP(k1,k2,k3,q1,q2,r1,r2)
      implicit none
      double complex aMPP
      double precision k1(4), k2(4),k3(4),q1(4),q2(4)
      double precision r1(4),r2(4),kq(4), kqb(4) 
      double precision dotp
      double complex spa,spb
      double precision m
#include "gggttbMPP-dec.F"
***#include "../common.F"
#ifdef DEBUG
      write(*,*) 'aMPP :'
#endif
      call add(q1,q2,kq) 
      call add(r1,r2,kqb) 
      m = dsqrt(2.d0 * dotp(q1,q2))
#ifdef DEBUG2 
      write(*,*) 'aMPP : m = ',m   
      write(*,*) 'aMPP : spa(k1,k2) = ', spa(k1,k2)  
      write(*,*) 'aMPP : spa(k1,k3) = ', spa(k1,k3)  
      write(*,*) 'aMPP : spa(k2,k3) = ', spa(k2,k3) 
      write(*,*) 'aMPP : spa(k1,q1) = ', spa(k1,q1)  
      write(*,*) 'aMPP : spa(k1,q2) = ', spa(k1,q2)  
      write(*,*) 'aMPP : spa(k2,q1) = ', spa(k2,q1)  
      write(*,*) 'aMPP : spa(k2,q2) = ', spa(k2,q2)  
      write(*,*) 'aMPP : spa(k3,q1) = ', spa(k3,q1)  
      write(*,*) 'aMPP : spa(k3,q2) = ', spa(k3,q2)  
      write(*,*) 'aMPP : spa(q1,q2) = ', spa(q1,q2)  
#endif
#include "gggttbMPP.F" 
#ifdef DEBUG 
      write(*,*) 'aMPP : return value', aMPP  
#endif 
      end
      function aMPM(k1,k2,k3,q1,q2,r1,r2)
      implicit none
      double complex aMPM
      double precision k1(4), k2(4),k3(4),q1(4),q2(4)
      double precision r1(4),r2(4),kq(4), kqb(4) 
      double precision dotp
      double complex spa,spb
      double precision m
#include "gggttbMPM-dec.F"
***#include "../common.F"
#ifdef DEBUG
      write(*,*) 'aMPM :'
#endif
      call add(q1,q2,kq) 
      call add(r1,r2,kqb) 
      m = dsqrt(2.d0 * dotp(q1,q2))
#ifdef DEBUG2 
      write(*,*) 'aMPM : m = ',m   
      write(*,*) 'aMPM : spa(k1,k2) = ', spa(k1,k2)  
      write(*,*) 'aMPM : spa(k1,k3) = ', spa(k1,k3)  
      write(*,*) 'aMPM : spa(k2,k3) = ', spa(k2,k3) 
      write(*,*) 'aMPM : spa(k1,q1) = ', spa(k1,q1)  
      write(*,*) 'aMPM : spa(k1,q2) = ', spa(k1,q2)  
      write(*,*) 'aMPM : spa(k2,q1) = ', spa(k2,q1)  
      write(*,*) 'aMPM : spa(k2,q2) = ', spa(k2,q2)  
      write(*,*) 'aMPM : spa(k3,q1) = ', spa(k3,q1)  
      write(*,*) 'aMPM : spa(k3,q2) = ', spa(k3,q2)  
      write(*,*) 'aMPM : spa(q1,q2) = ', spa(q1,q2)  
#endif
#include "gggttbMPM.F" 
#ifdef DEBUG 
      write(*,*) 'aMPM : return value', aMPM  
#endif 
      end
      function aPMM(k1,k2,k3,q1,q2,r1,r2)
      implicit none
      double complex aPMM
      double precision k1(4), k2(4),k3(4),q1(4),q2(4)
      double precision r1(4),r2(4),kq(4), kqb(4) 
      double precision dotp
      double complex spa,spb
      double precision m
#include "gggttbPMM-dec.F"
***#include "../common.F"
#ifdef DEBUG
      write(*,*) 'aPMM :'
#endif
      call add(q1,q2,kq) 
      call add(r1,r2,kqb) 
      m = dsqrt(2.d0 * dotp(q1,q2))
#ifdef DEBUG2 
      write(*,*) 'aPMM : m = ',m   
      write(*,*) 'aPMM : spa(k1,k2) = ', spa(k1,k2)  
      write(*,*) 'aPMM : spa(k1,k3) = ', spa(k1,k3)  
      write(*,*) 'aPMM : spa(k2,k3) = ', spa(k2,k3) 
      write(*,*) 'aPMM : spa(k1,q1) = ', spa(k1,q1)  
      write(*,*) 'aPMM : spa(k1,q2) = ', spa(k1,q2)  
      write(*,*) 'aPMM : spa(k2,q1) = ', spa(k2,q1)  
      write(*,*) 'aPMM : spa(k2,q2) = ', spa(k2,q2)  
      write(*,*) 'aPMM : spa(k3,q1) = ', spa(k3,q1)  
      write(*,*) 'aPMM : spa(k3,q2) = ', spa(k3,q2)  
      write(*,*) 'aPMM : spa(q1,q2) = ', spa(q1,q2)  
#endif
#include "gggttbPMM.F" 
#ifdef DEBUG 
      write(*,*) 'aPMM : return value', aPMM  
#endif 
      end
      function aPPM(k1,k2,k3,q1,q2,r1,r2)
      implicit none
      double complex aPPM
      double precision k1(4), k2(4),k3(4),q1(4),q2(4)
      double precision r1(4),r2(4),kq(4), kqb(4) 
      double precision dotp
      double complex spa,spb
      double precision m
#include "gggttbPPM-dec.F"
***#include "../common.F"
#ifdef DEBUG
      write(*,*) 'aPPM :'
#endif
      call add(q1,q2,kq) 
      call add(r1,r2,kqb) 
      m = dsqrt(2.d0 * dotp(q1,q2))
#ifdef DEBUG2 
      write(*,*) 'aPPM : m = ',m   
      write(*,*) 'aPPM : spa(k1,k2) = ', spa(k1,k2)  
      write(*,*) 'aPPM : spa(k1,k3) = ', spa(k1,k3)  
      write(*,*) 'aPPM : spa(k2,k3) = ', spa(k2,k3) 
      write(*,*) 'aPPM : spa(k1,q1) = ', spa(k1,q1)  
      write(*,*) 'aPPM : spa(k1,q2) = ', spa(k1,q2)  
      write(*,*) 'aPPM : spa(k2,q1) = ', spa(k2,q1)  
      write(*,*) 'aPPM : spa(k2,q2) = ', spa(k2,q2)  
      write(*,*) 'aPPM : spa(k3,q1) = ', spa(k3,q1)  
      write(*,*) 'aPPM : spa(k3,q2) = ', spa(k3,q2)  
      write(*,*) 'aPPM : spa(q1,q2) = ', spa(q1,q2)  
#endif
#include "gggttbPPM.F" 
#ifdef DEBUG 
      write(*,*) 'aPPM : return value', aPPM  
#endif 
      end
      function aPMP(k1,k2,k3,q1,q2,r1,r2)
      implicit none
      double complex aPMP
      double precision k1(4), k2(4),k3(4),q1(4),q2(4)
      double precision r1(4),r2(4),kq(4), kqb(4) 
      double precision dotp
      double complex spa,spb
      double precision m
#include "gggttbPMP-dec.F"
***#include "../common.F"
#ifdef DEBUG
      write(*,*) 'aPMP :'
#endif
      call add(q1,q2,kq) 
      call add(r1,r2,kqb) 
      m = dsqrt(2.d0 * dotp(q1,q2))
#ifdef DEBUG2 
      write(*,*) 'aPMP : m = ',m   
      write(*,*) 'aPMP : spa(k1,k2) = ', spa(k1,k2)  
      write(*,*) 'aPMP : spa(k1,k3) = ', spa(k1,k3)  
      write(*,*) 'aPMP : spa(k2,k3) = ', spa(k2,k3) 
      write(*,*) 'aPMP : spa(k1,q1) = ', spa(k1,q1)  
      write(*,*) 'aPMP : spa(k1,q2) = ', spa(k1,q2)  
      write(*,*) 'aPMP : spa(k2,q1) = ', spa(k2,q1)  
      write(*,*) 'aPMP : spa(k2,q2) = ', spa(k2,q2)  
      write(*,*) 'aPMP : spa(k3,q1) = ', spa(k3,q1)  
      write(*,*) 'aPMP : spa(k3,q2) = ', spa(k3,q2)  
      write(*,*) 'aPMP : spa(q1,q2) = ', spa(q1,q2)  
#endif
#include "gggttbPMP.F" 
#ifdef DEBUG 
      write(*,*) 'aPMP : return value', aPMP  
#endif 
      end
      function aPPP(k1,k2,k3,q1,q2,r1,r2)
      implicit none
      double complex aPPP
      double precision k1(4), k2(4),k3(4),q1(4),q2(4)
      double precision r1(4),r2(4),kq(4), kqb(4) 
      double precision dotp
      double complex spa,spb
      double precision m
#include "gggttbPPP-dec.F"
***#include "../common.F"
#ifdef DEBUG
      write(*,*) 'aPPP :'
#endif
      call add(q1,q2,kq) 
      call add(r1,r2,kqb) 
      m = dsqrt(2.d0 * dotp(q1,q2))
#ifdef DEBUG2 
      write(*,*) 'aPPP : m = ',m   
      write(*,*) 'aPPP : spa(k1,k2) = ', spa(k1,k2)  
      write(*,*) 'aPPP : spa(k1,k3) = ', spa(k1,k3)  
      write(*,*) 'aPPP : spa(k2,k3) = ', spa(k2,k3) 
      write(*,*) 'aPPP : spa(k1,q1) = ', spa(k1,q1)  
      write(*,*) 'aPPP : spa(k1,q2) = ', spa(k1,q2)  
      write(*,*) 'aPPP : spa(k2,q1) = ', spa(k2,q1)  
      write(*,*) 'aPPP : spa(k2,q2) = ', spa(k2,q2)  
      write(*,*) 'aPPP : spa(k3,q1) = ', spa(k3,q1)  
      write(*,*) 'aPPP : spa(k3,q2) = ', spa(k3,q2)  
      write(*,*) 'aPPP : spa(q1,q2) = ', spa(q1,q2)  
#endif
#include "gggttbPPP.F" 
#ifdef DEBUG 
      write(*,*) 'aPPP : return value', aPPP  
#endif 
      end
