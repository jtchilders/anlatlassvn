C
C    This is automatically created code, produced by ./mkcode.csh
C    on Thu Oct  9 11:27:14 MET DST 2003.
C
      function den(x)
      implicit none
      double precision den,x
      den=1.d0/x
      end
      function aMMMM(k1,k2,k3,k4,q1,q2,r1,r2)
      implicit none
      double complex aMMMM
      double precision k1(4), k2(4),k3(4),k4,q1(4),q2(4)
      double precision r1(4),r2(4),kq(4), kqb(4) 
      double precision dotp,den
      double complex spa,spb
      double precision m
      double complex s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
#include "ggggttMMMM-dec.F"
#ifdef DEBUG
      write(*,*) 'aMMMM :'
#endif
      call add(q1,q2,kq) 
      call add(r1,r2,kqb) 
      m = dsqrt(dotp(kq,kq))
#include "ggggttMMMM.F" 
#ifdef DEBUG 
      write(*,*) 'aMMMM : return value', aMMMM  
#endif 
      end
      function aMMMP(k1,k2,k3,k4,q1,q2,r1,r2)
      implicit none
      double complex aMMMP
      double precision k1(4), k2(4),k3(4),k4,q1(4),q2(4)
      double precision r1(4),r2(4),kq(4), kqb(4) 
      double precision dotp,den
      double complex spa,spb
      double precision m
      double complex s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
#include "ggggttMMMP-dec.F"
#ifdef DEBUG
      write(*,*) 'aMMMP :'
#endif
      call add(q1,q2,kq) 
      call add(r1,r2,kqb) 
      m = dsqrt(dotp(kq,kq))
#include "ggggttMMMP.F" 
#ifdef DEBUG 
      write(*,*) 'aMMMP : return value', aMMMP  
#endif 
      end
      function aMMPM(k1,k2,k3,k4,q1,q2,r1,r2)
      implicit none
      double complex aMMPM
      double precision k1(4), k2(4),k3(4),k4,q1(4),q2(4)
      double precision r1(4),r2(4),kq(4), kqb(4) 
      double precision dotp,den
      double complex spa,spb
      double precision m
      double complex s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
#include "ggggttMMPM-dec.F"
#ifdef DEBUG
      write(*,*) 'aMMPM :'
#endif
      call add(q1,q2,kq) 
      call add(r1,r2,kqb) 
      m = dsqrt(dotp(kq,kq))
#include "ggggttMMPM.F" 
#ifdef DEBUG 
      write(*,*) 'aMMPM : return value', aMMPM  
#endif 
      end
      function aMMPP(k1,k2,k3,k4,q1,q2,r1,r2)
      implicit none
      double complex aMMPP
      double precision k1(4), k2(4),k3(4),k4,q1(4),q2(4)
      double precision r1(4),r2(4),kq(4), kqb(4) 
      double precision dotp,den
      double complex spa,spb
      double precision m
      double complex s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
#include "ggggttMMPP-dec.F"
#ifdef DEBUG
      write(*,*) 'aMMPP :'
#endif
      call add(q1,q2,kq) 
      call add(r1,r2,kqb) 
      m = dsqrt(dotp(kq,kq))
#include "ggggttMMPP.F" 
#ifdef DEBUG 
      write(*,*) 'aMMPP : return value', aMMPP  
#endif 
      end
      function aMPMM(k1,k2,k3,k4,q1,q2,r1,r2)
      implicit none
      double complex aMPMM
      double precision k1(4), k2(4),k3(4),k4,q1(4),q2(4)
      double precision r1(4),r2(4),kq(4), kqb(4) 
      double precision dotp,den
      double complex spa,spb
      double precision m
      double complex s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
#include "ggggttMPMM-dec.F"
#ifdef DEBUG
      write(*,*) 'aMPMM :'
#endif
      call add(q1,q2,kq) 
      call add(r1,r2,kqb) 
      m = dsqrt(dotp(kq,kq))
#include "ggggttMPMM.F" 
#ifdef DEBUG 
      write(*,*) 'aMPMM : return value', aMPMM  
#endif 
      end
      function aMPMP(k1,k2,k3,k4,q1,q2,r1,r2)
      implicit none
      double complex aMPMP
      double precision k1(4), k2(4),k3(4),k4,q1(4),q2(4)
      double precision r1(4),r2(4),kq(4), kqb(4) 
      double precision dotp,den
      double complex spa,spb
      double precision m
      double complex s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
#include "ggggttMPMP-dec.F"
#ifdef DEBUG
      write(*,*) 'aMPMP :'
#endif
      call add(q1,q2,kq) 
      call add(r1,r2,kqb) 
      m = dsqrt(dotp(kq,kq))
#include "ggggttMPMP.F" 
#ifdef DEBUG 
      write(*,*) 'aMPMP : return value', aMPMP  
#endif 
      end
      function aMPPM(k1,k2,k3,k4,q1,q2,r1,r2)
      implicit none
      double complex aMPPM
      double precision k1(4), k2(4),k3(4),k4,q1(4),q2(4)
      double precision r1(4),r2(4),kq(4), kqb(4) 
      double precision dotp,den
      double complex spa,spb
      double precision m
      double complex s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
#include "ggggttMPPM-dec.F"
#ifdef DEBUG
      write(*,*) 'aMPPM :'
#endif
      call add(q1,q2,kq) 
      call add(r1,r2,kqb) 
      m = dsqrt(dotp(kq,kq))
#include "ggggttMPPM.F" 
#ifdef DEBUG 
      write(*,*) 'aMPPM : return value', aMPPM  
#endif 
      end
      function aMPPP(k1,k2,k3,k4,q1,q2,r1,r2)
      implicit none
      double complex aMPPP
      double precision k1(4), k2(4),k3(4),k4,q1(4),q2(4)
      double precision r1(4),r2(4),kq(4), kqb(4) 
      double precision dotp,den
      double complex spa,spb
      double precision m
      double complex s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
#include "ggggttMPPP-dec.F"
#ifdef DEBUG
      write(*,*) 'aMPPP :'
#endif
      call add(q1,q2,kq) 
      call add(r1,r2,kqb) 
      m = dsqrt(dotp(kq,kq))
#include "ggggttMPPP.F" 
#ifdef DEBUG 
      write(*,*) 'aMPPP : return value', aMPPP  
#endif 
      end
      function aPMMM(k1,k2,k3,k4,q1,q2,r1,r2)
      implicit none
      double complex aPMMM
      double precision k1(4), k2(4),k3(4),k4,q1(4),q2(4)
      double precision r1(4),r2(4),kq(4), kqb(4) 
      double precision dotp,den
      double complex spa,spb
      double precision m
      double complex s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
#include "ggggttPMMM-dec.F"
#ifdef DEBUG
      write(*,*) 'aPMMM :'
#endif
      call add(q1,q2,kq) 
      call add(r1,r2,kqb) 
      m = dsqrt(dotp(kq,kq))
#include "ggggttPMMM.F" 
#ifdef DEBUG 
      write(*,*) 'aPMMM : return value', aPMMM  
#endif 
      end
      function aPMMP(k1,k2,k3,k4,q1,q2,r1,r2)
      implicit none
      double complex aPMMP
      double precision k1(4), k2(4),k3(4),k4,q1(4),q2(4)
      double precision r1(4),r2(4),kq(4), kqb(4) 
      double precision dotp,den
      double complex spa,spb
      double precision m
      double complex s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
#include "ggggttPMMP-dec.F"
#ifdef DEBUG
      write(*,*) 'aPMMP :'
#endif
      call add(q1,q2,kq) 
      call add(r1,r2,kqb) 
      m = dsqrt(dotp(kq,kq))
#include "ggggttPMMP.F" 
#ifdef DEBUG 
      write(*,*) 'aPMMP : return value', aPMMP  
#endif 
      end
      function aPMPM(k1,k2,k3,k4,q1,q2,r1,r2)
      implicit none
      double complex aPMPM
      double precision k1(4), k2(4),k3(4),k4,q1(4),q2(4)
      double precision r1(4),r2(4),kq(4), kqb(4) 
      double precision dotp,den
      double complex spa,spb
      double precision m
      double complex s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
#include "ggggttPMPM-dec.F"
#ifdef DEBUG
      write(*,*) 'aPMPM :'
#endif
      call add(q1,q2,kq) 
      call add(r1,r2,kqb) 
      m = dsqrt(dotp(kq,kq))
#include "ggggttPMPM.F" 
#ifdef DEBUG 
      write(*,*) 'aPMPM : return value', aPMPM  
#endif 
      end
      function aPMPP(k1,k2,k3,k4,q1,q2,r1,r2)
      implicit none
      double complex aPMPP
      double precision k1(4), k2(4),k3(4),k4,q1(4),q2(4)
      double precision r1(4),r2(4),kq(4), kqb(4) 
      double precision dotp,den
      double complex spa,spb
      double precision m
      double complex s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
#include "ggggttPMPP-dec.F"
#ifdef DEBUG
      write(*,*) 'aPMPP :'
#endif
      call add(q1,q2,kq) 
      call add(r1,r2,kqb) 
      m = dsqrt(dotp(kq,kq))
#include "ggggttPMPP.F" 
#ifdef DEBUG 
      write(*,*) 'aPMPP : return value', aPMPP  
#endif 
      end
      function aPPMM(k1,k2,k3,k4,q1,q2,r1,r2)
      implicit none
      double complex aPPMM
      double precision k1(4), k2(4),k3(4),k4,q1(4),q2(4)
      double precision r1(4),r2(4),kq(4), kqb(4) 
      double precision dotp,den
      double complex spa,spb
      double precision m
      double complex s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
#include "ggggttPPMM-dec.F"
#ifdef DEBUG
      write(*,*) 'aPPMM :'
#endif
      call add(q1,q2,kq) 
      call add(r1,r2,kqb) 
      m = dsqrt(dotp(kq,kq))
#include "ggggttPPMM.F" 
#ifdef DEBUG 
      write(*,*) 'aPPMM : return value', aPPMM  
#endif 
      end
      function aPPMP(k1,k2,k3,k4,q1,q2,r1,r2)
      implicit none
      double complex aPPMP
      double precision k1(4), k2(4),k3(4),k4,q1(4),q2(4)
      double precision r1(4),r2(4),kq(4), kqb(4) 
      double precision dotp,den
      double complex spa,spb
      double precision m
      double complex s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
#include "ggggttPPMP-dec.F"
#ifdef DEBUG
      write(*,*) 'aPPMP :'
#endif
      call add(q1,q2,kq) 
      call add(r1,r2,kqb) 
      m = dsqrt(dotp(kq,kq))
#include "ggggttPPMP.F" 
#ifdef DEBUG 
      write(*,*) 'aPPMP : return value', aPPMP  
#endif 
      end
      function aPPPM(k1,k2,k3,k4,q1,q2,r1,r2)
      implicit none
      double complex aPPPM
      double precision k1(4), k2(4),k3(4),k4,q1(4),q2(4)
      double precision r1(4),r2(4),kq(4), kqb(4) 
      double precision dotp,den
      double complex spa,spb
      double precision m
      double complex s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
#include "ggggttPPPM-dec.F"
#ifdef DEBUG
      write(*,*) 'aPPPM :'
#endif
      call add(q1,q2,kq) 
      call add(r1,r2,kqb) 
      m = dsqrt(dotp(kq,kq))
#include "ggggttPPPM.F" 
#ifdef DEBUG 
      write(*,*) 'aPPPM : return value', aPPPM  
#endif 
      end
      function aPPPP(k1,k2,k3,k4,q1,q2,r1,r2)
      implicit none
      double complex aPPPP
      double precision k1(4), k2(4),k3(4),k4,q1(4),q2(4)
      double precision r1(4),r2(4),kq(4), kqb(4) 
      double precision dotp,den
      double complex spa,spb
      double precision m
      double complex s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
#include "ggggttPPPP-dec.F"
#ifdef DEBUG
      write(*,*) 'aPPPP :'
#endif
      call add(q1,q2,kq) 
      call add(r1,r2,kqb) 
      m = dsqrt(dotp(kq,kq))
#include "ggggttPPPP.F" 
#ifdef DEBUG 
      write(*,*) 'aPPPP : return value', aPPPP  
#endif 
      end
