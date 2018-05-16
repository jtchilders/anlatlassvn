C$Modified: Mon Jan  8 15:59:56 2007 by puwer $

C----------------------------------------------------------------------
      subroutine eventn(n,s,m,x,pset,jacobi,err)
C----------------------------------------------------------------------
C     Subroutine to generate a n-parton event with center-of-mass
C     energy squared given by s and the parton masses set through m(n)
C     The momentum configuration is returned in pset(4,n) (the energy
C     is stored in the fourth component :-( ), the jacobian for
C     the transformation is returned through jacobi.
C
      implicit none
      integer n,err
      double precision pset(4,n),bvec(4),m(n)
      double precision s,sqrts,x(3*n-4),jacobi
      double precision costheta(n),sintheta(n),
     &     phi(n),sinphi(n),cosphi(n),MM(n),MMQ(n),mu(n),P(n)
      integer i,j,irn,rho
      double precision lam
      double precision dotp
      double precision lambda,a,b,c
      logical init
      data init/.true./
      save init
#include "constants.F"      
      lambda(a,b,c) = a**2+b**2+c**2-2.d0*a*b-2.d0*b*c-2.d0*a*c

      err = 0

#ifdef VERBOSE
      if (init) then
         write(*,*) ' -------------------------------------'
         write(*,*) ' Simple event generator'
         write(*,*) ' written by P.Uwer '
         write(*,*) 
         write(*,*) ' [0,1]^(3n-4) --> (p1,...,pn)         '
         write(*,*) ' based on Byckling&Kajantie Chap 10   '
         write(*,*) ' -------------------------------------'
         init = .false.
      endif
#endif

      sqrts = dsqrt(s)

      do i=1,n
         mu(i) = 0.d0
         do j=1,i
            mu(i) = mu(i) + m(j) 
         enddo
      enddo


      irn = 1
C
C     Set the angles
C     
      do i=2,n
         phi(i) = 2.d0 * pi * x(irn)
         sinphi(i) = dsin(phi(i))
C     Note: We have to calculate the cosine explicitly!!!!
C     dsqrt(1-sin^2) would only give cos() > 0
         cosphi(i) = dcos(phi(i))
         costheta(i) = 2.d0 * x(irn+1) - 1.d0
         sintheta(i) = dsqrt(1.d0-costheta(i)**2)
         irn = irn + 2
      enddo
 
C
C     Set the invariant masses MMQ(i) = (p_1+...+p_i)^2
C     

      MM(n) = sqrts
      MMQ(n) = s
      MM(1) = m(1)
      MMQ(1) = MM(1)**2
      do i=n-1,2,-1
         MM(i) = mu(i) + x(irn) * ( MM(i+1) - mu(i+1) )
         MMQ(i) = MM(i)**2
         irn = irn + 1
      enddo
C
C     Caclulate the three momenta
C
      do i=2,n
         if (m(i).eq.0.d0) then
            P(i) = (MMQ(i)-MMQ(i-1)) / (2.d0 * MM(i))
         else
            P(i) = dsqrt(lambda(MMQ(i),MMQ(i-1),m(i)**2)) 
     .           / (2.d0 * MM(i))
         endif
      enddo

C
C     calculate the jacobian for the transformation:
C
      jacobi = 1.d0
      do i=2,n
         jacobi = jacobi * (2.d0*pi) * p(i)
      enddo
      do i=3,n
         jacobi = jacobi * ( MM(i) - mu(i) )
      enddo
      jacobi = jacobi/(2.d0*sqrts)

C
C     Finally reconstruct the event in the overall CMS:
C
C     Calculate first the 4-momenta pset(rho,i) in the restframe of
C     p_1+...+p_i:
C
      do i=2,n
         pset(1,i) = P(i) * cosphi(i)*sintheta(i)
         pset(2,i) = P(i) * sinphi(i)*sintheta(i)
         pset(3,i) = P(i) * costheta(i)
         pset(4,i) = dsqrt(m(i)**2+p(i)**2)
      enddo

      pset(4,1) = MM(2) - pset(4,2)
      pset(1,1) = -pset(1,2)
      pset(2,1) = -pset(2,2)
      pset(3,1) = -pset(3,2)

      do i=3,n
         bvec(4) = MM(i)-pset(4,i)
         do rho=1,3
            bvec(rho)=-pset(rho,i)
         enddo
         do j = 1,i-1
            call boost(bvec,pset(1,j),pset(1,j),err)
         enddo
      enddo

      end
C-----------------------------------------------------------------------      

























