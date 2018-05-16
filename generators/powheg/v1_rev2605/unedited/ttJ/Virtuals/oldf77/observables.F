C$Modified: Mon Aug  9 15:56:46 2004 by puwer $
      function obs4(pset)
      implicit none
      double precision obs4,pset(4,-2:4),pout(4,-2:3),p2set(4,-2:2)
      double precision beam(4),p1(4),p2(4),p3(4),p4(4),
     .     pgjet(4),beam_h(4),beam_t1(4),beam_t2(4)
      integer icount,mu,njet
      double precision phi3,phi4,kt3,kt4,kt34,eta3,eta4
      double precision evalkt,dot3p,norm3,rmax,rmin,X,datanh,help
      double precision DeltaPhi2
      integer durham4,durham2
#include "common.F"   
      datanh(x) = dlog((1.d0+x)/(1.d0-x))/2.d0
      if (iobs.eq.0) then
C     for testing:
C     no cut at all
         obs4 = 1.d0
         return
      endif
C
      if (iobs.eq.102) then
C
C     Use the definition of Ellis and Soper to define the `gluonic'
C     jets. The top's are assumed to be always tagged.
C
         do mu=1,4
            beam(mu) = pset(mu,-1)/pset(4,-1)
            p3(mu) = pset(mu,3)
            p4(mu) = pset(mu,4)
         enddo
C
C     We don't assume anything about the beam direction:
C
         beam_h(4) = 0.d0
         beam_h(1) = beam(3)
         beam_h(2) = beam(1)
         beam_h(3) = beam(2)
         


         call cross3p(beam,beam_h,beam_t1)
         help = norm3(beam_t1)
         beam_t1(1) = beam_t1(1) / help
         beam_t1(2) = beam_t1(2) / help
         beam_t1(3) = beam_t1(3) / help
         call cross3p(beam,beam_t1,beam_t2)
         
         kt3 = evalkt(beam,p3)
         kt4 = evalkt(beam,p4)
         eta3 = datanh(dot3p(p3,beam)/p3(4))
         eta4 = datanh(dot3p(p4,beam)/p4(4))
C     We can also use
C         eta3 = -dlog(dtan(0.5d0*dacos(dot3p(p3,beam)/p3(4)))) 
C         eta4 = -dlog(dtan(0.5d0*dacos(dot3p(p4,beam)/p4(4)))) 
C     to determine eta3 and eta4.
C
C     With the angle phi we have to be carefull to cover the
C     full allowed range from 0 to 2 Pi or -Pi to Pi.
C     Furthermore when calculating the difference we should make
C     sure that the difference is always .le. Pi!
         phi3 =  datan2(dot3p(beam_t2,p3),dot3p(beam_t1,p3))
         phi4 =  datan2(dot3p(beam_t2,p4),dot3p(beam_t1,p4))

         if (((eta3-eta4)**2+DeltaPhi2(phi3,phi4)).lt.1.d0) then
C     
C     Merge parton 3 and parton 4
C
            kt34 = (kt3+kt4) 
            if (kt34.gt.ktmin) then
               obs4=1.d0
               return
            else
               obs4=0.d0
               return
            endif
         else
C
C     Two gluons remain as individual jets
C
            if (rmax(kt3,kt4).gt.ktmin) then
               obs4 = 1.d0
               return
            else
               obs4 = 0.d0
               return
            endif
         endif
      endif
      write(*,*) 'iobs = ',iobs,' not definied, stop program'
      stop
      end
C********************************************************************
      function obs3(pset)
      implicit none
      double precision obs3,pset(4,-2:3),pout(4,-2:2)
      double precision evalkt
      double precision p1(4),p2(4),p3(4),beam(4)
      integer mu 
      integer durham3
#include "common.F"      
      if (iobs.eq.0) then 
C     for testing
         obs3 = 1.d0
         return
      endif
      if (iobs.eq.102) then
C     we demand a minimal kt>ktmin GeV
C     we assume that the tops are always tagged
         obs3 = 0.d0
         do mu=1,4
            beam(mu) = pset(mu,-1)/pset(4,-1)
**            p1(mu) = pset(mu,1)
**            p2(mu) = pset(mu,2)
            p3(mu) = pset(mu,3)
         enddo
**         if (evalkt(beam,p1).lt.ktmin) return
**         if (evalkt(beam,p2).lt.ktmin) return
         if (evalkt(beam,p3).lt.ktmin) return
         obs3 = 1.d0
         return
      endif

      write(*,*) 'iobs = ',iobs,'not definied, stop program'
      stop
      end
C-------------------------------------------------------------------
C     calculate transversal momentum of p2 with respect to p1
C     p1 is assumed to be normalized
C-------------------------------------------------------------------
      function evalkt(p1,p2)
      double precision evalkt
      double precision p1(4),p2(4)
      double precision dot3p
C     total momentum minus longitudinal gives the transversal:
      evalkt = dsqrt( dot3p(p2,p2) - dot3p(p1,p2)**2 )
      end
      
      function DeltaPhi2(phi1,phi2)
      implicit none
      double precision DeltaPhi2,phi1,phi2
#include "constants.F"
      DeltaPhi2 = dabs(phi1-phi2)
      if (DeltaPhi2.gt.pi) DeltaPhi2 = 2.d0 * pi - DeltaPhi2
      DeltaPhi2 = Deltaphi2**2
      end
