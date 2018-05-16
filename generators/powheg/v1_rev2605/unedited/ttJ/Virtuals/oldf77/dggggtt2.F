C$Modified: Wed Jan  3 14:38:29 2007 by puwer $
C*******************************************************************
      function dggggtt2(x,wgt)
      implicit none
      double precision dggggtt2,x(8),wgt
      double precision pset(4,-2:4),p(4,4)
      double precision pset2(4,-2:3)
      double precision jacobi        ! weight 
      double precision masses(4)
      DOUBLE PRECISION rt12,rt13,rt14,rt23,rt24,rt34,rm1,rm2,rm3,rm4
      double precision obs4, subtractions,mggggtt2,D,dotp
      double precision beam(4),ksum(4)
      double precision k1m(0:3), k2m(0:3), k3m(0:3),k4m(0:3),
     &     kqm(0:3),kqbm(0:3)
      double precision kq(4),kqb(4),k1(4),k2(4),ka(4),kb(4)
      double precision  SGG_GGTTB
      DOUBLE PRECISION EVALKT
      double precision dipoles
      double precision prefac
      integer cut4part,durham4
      integer mu,ncount,ihelp,i,j,itest
      integer i1,i2,i3,i4
      integer npoint
      double precision smin
      double precision min,max,msquare
      data npoint/1/,min/1d300/,max/-1d300/
      save npoint,min,max
#include "common.F"
#include "constants.F"
      common/event/pset
      parameter (prefac =   1.d0
C     Factors from phase space measure:
     &     /2.d0/(2.d0*pi)**8
C     Two identical gluons in the final state
     &     /2.d0
C     Average over color and spins:
     &     / 64.d0 /4.d0
C     Converting to pbarn
     &     * hcq * 1.d12 )

C
      ierr = 0
      dggggtt2 = 0.d0
      ncount = 4
      masses(1) = mt
      masses(2) = mt
      masses(3) = 0.d0
      masses(4) = 0.d0


      pset(1,-1) = 0.d0
      pset(2,-1) = sqrt_s/2.d0
      pset(3,-1) = 0.d0
      pset(4,-1) = sqrt_s/2.d0

      pset(1,-2) = 0.d0
      pset(2,-2) = -sqrt_s/2.d0
      pset(3,-2) = 0.d0
      pset(4,-2) = sqrt_s/2.d0


      call event4(s,masses,x,p,jacobi,ierr)
*      call eventn(4,s,masses,x,p,jacobi,ierr)



      if (ierr.lt.0) then
         write(logunit,*) 'dggggtt2: error in event4, drop event'
         call flush(logunit)
         return
      endif

      do j=1,4
         do mu=1,4
            pset(mu,j) = p(mu,j)
         enddo
      enddo
      
      smin = sqrt_s**2
      
      do i=-2,4
         do j=i+1,4
            if ((i*j.ne.0).and.(dotp(pset(1,i),pset(1,j)).lt.smin)) 
     &           then
               smin = dotp(pset(1,i),pset(1,j))
            endif
         enddo
      enddo

      if (smin.lt.1d-4) return

#ifdef DEBUG
      rt12 = 2.d0*dotp(p(1,1),p(1,2))
      rt13 = 2.d0*dotp(p(1,1),p(1,3))
      rt14 = 2.d0*dotp(p(1,1),p(1,4))
      rt23 = 2.d0*dotp(p(1,2),p(1,3))
      rt24 = 2.d0*dotp(p(1,2),p(1,4))
      rt34 = 2.d0*dotp(p(1,3),p(1,4))
      rm1 = mt
      rm2 = mt 
      rm3 = 0.d0
      rm4 = 0.d0
#endif

      dggggtt2 = obs4(pset)

      if (dggggtt2.ne.0.d0) then
#ifdef MADGRAPH
      do mu=1,4  
         k1m( mod(mu,4) ) = pset(mu,-1)
         k2m( mod(mu,4) ) = pset(mu,-2)
         kqm( mod(mu,4) ) = pset(mu,1)
         kqbm( mod(mu,4) ) = pset(mu,2) 
         k3m( mod(mu,4) ) = pset(mu,3)
         k4m( mod(mu,4) ) = pset(mu,4)
      enddo
C     Note the addional factor of 2, in Madgraph the symmetriefactor 1/2
C     is already included in the Matrixelement while we include it in the 
C     phase-space
      dggggtt2 = dggggtt2 * 4.d0 * 64.d0 * 2.d0
     &     *  SGG_GGTTB(k1m, k2m, k3m, k4m, kqm, kqbm)
#else
      dggggtt2 = dggggtt2 * mggggtt2(pset) * g2**4 
#endif
      endif
C
      msquare = dggggtt2
C
      dipoles = subtractions(pset)

C
      if (ierr.lt.0) then
         write(*,*) 'dggggtt2: error in subtractions, drop event'
         write(*,*) 'ierr =', ierr
         dggggtt2 = 0.d0
         call printevent(pset)
         write(*,*) 'sqrt_s = ', sqrt_s
         do i=1,8
            write(*,*) 'x(',i,')= ',x(i)
         enddo
         return
      endif

*----------------------------------------------------
*      write(*,*) 'M2'
*      write(*,*)  msquare
*      write(*,*) 'Dipoles'
*      write(*,*) dipoles
*----------------------------------------------------
C
      if (dggggtt2.eq.0.d0) then
         dggggtt2 = dipoles
      else
         if ( dggggtt2 .gt. 10**7*dabs(dggggtt2 + dipoles)) then
            dggggtt2 = 0.d0
            write(logunit,*) ' 7 digits loss in precision, drop event' 
            return
         else
            dggggtt2 = dggggtt2 + dipoles
         endif
      endif
      
C
C     That's it:
      dggggtt2 = dggggtt2 * prefac * jacobi / s 

      end
