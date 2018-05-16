C
C     $Modified: Sun Dec 31 13:48:30 2006 by puwer $
C
      double precision function ggttgg2c(p1,p2,kt,ktb,p3,p4)
      implicit none
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      double precision kt(0:3),ktb(0:3)
      double precision pset(4,-2:4)
      double precision mggggtt2
      integer i
      do i = 1,4
         pset(i,-1) = p1(mod(i,4))
         pset(i,-2) = p2(mod(i,4))
         pset(i,1) = kt(mod(i,4))
         pset(i,2) = ktb(mod(i,4))
         pset(i,3) = p3(mod(i,4))
         pset(i,4) = p4(mod(i,4))
      enddo
      ggttgg2C = mggggtt2(pset)
      end
C------------------------------------------------------
      subroutine setobs(ii)
      implicit none
      integer ii
#include "common.F"
#include "constants.F"
      iobs = ii
      write(*,'(A,$)') '# Observable set to iobs = '
      write(*,'(I3)') iobs
      end
C------------------------------------------------------
      subroutine setshat(x)
      implicit none 
      double precision x
#include "common.F"
#include "constants.F"
      s = x
      sqrt_s = dsqrt(s)
      end
C------------------------------------------------------
      subroutine setmt(x)
      implicit none 
      double precision x
#include "common.F"
#include "constants.F"
      mt = x
      end
C------------------------------------------------------
      subroutine setktmin(x)
      implicit none 
      double precision x
#include "common.F"
#include "constants.F"
      ktmin = x 
      end
C------------------------------------------------------
      subroutine printevent(pset)
      double precision pset(4,-2:4)
      integer i,mu
      write(*,*) '*** Event *********************************'
      do i=-2,4
         do mu=1,4
            if (i.ne.0) write(*,'(A,I2,A,I2,A,E25.15)') 
     &           'pset(',mu,',',i,')=', pset(mu,i)
         enddo
      enddo
      write(*,*) '*******************************************'
      end
C------------------------------------------------------
      subroutine setupalphas(value)
      implicit none
      double precision value
#include "common.F"      
#include "constants.F"
      alphas = value
      g2 = 4.d0*pi*alphas
      write(*,'(A,$)') '# alphas set to: '
      write(*,*) alphas
#ifdef MADGRAPH
c     IF MADGRAPH IS USED WE HAVE TO REINITIALZE IT...
      write(*,*) '* Init Madgraph '
      call initialize
#endif
      end
