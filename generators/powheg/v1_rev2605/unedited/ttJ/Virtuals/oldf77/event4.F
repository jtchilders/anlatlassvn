C$Modified: Mon Aug  9 11:32:16 2004 by puwer $
C----------------------------------------------------------------------
      subroutine event4(s,m,x,pset,jacobi,err)
C----------------------------------------------------------------------
      implicit none
      integer n,err
      parameter (n=4)
      double precision pset(4,n),bvec(4),m(n)
      double precision s,x(3*n-4),jacobi
      double precision costheta(n),sintheta(n),
     &     phi(n),sinphi(n),cosphi(n),MM(n),MMQ(n),mu(n),P(n)
      integer i,j,irn,rho
      double precision lam
      double precision dotp
      double precision lambda,a,b,c
#include "constants.F"
      lambda(a,b,c) = a**2+b**2+c**2-2.d0*a*b-2.d0*b*c-2.d0*a*c


      do i=1,n
         mu(i) = 0.d0
         do j=1,i
            mu(i) = mu(i) + m(j) 
         enddo
      enddo


      irn = 1
      
C
C     set the angles
C
      do i=2,n
         phi(i) = 2.d0 * pi * x(irn)
         sinphi(i) = dsin(phi(i))
C     We have to calculate the cosine explicitly!!!!
C     dsqrt(1-sin^2) would only give cos() > 0
         cosphi(i) = dcos(phi(i))
         costheta(i) = 2.d0 * x(irn+1) - 1.d0
         sintheta(i) = dsqrt(1.d0-costheta(i)**2)
         irn = irn + 2
      enddo
 
C
C     Set the invariant masses MMQ(i) = (p_1+...+p_i)^2
C     
      MM(n) = dsqrt(s)
      MMQ(n) = s
      MM(1) = m(1)
      MMQ(1) = MM(1)**2
      do i=n-1,2,-1
         MM(i) = mu(i) + x(irn) * ( MM(i+1) - mu(i+1) )
         MMQ(i) = MM(i)**2
         irn = irn + 1;
      enddo
C
C     Caclulate the three momenta
C
      do i=2,n
         if (m(i).eq.0.d0) then
            lam=(MMQ(i)-MMQ(i-1))**2
         else
            lam=lambda(MMQ(i),MMQ(i-1),m(i)**2)
         endif
         P(i) = dsqrt(lam)/2.d0/MM(i)
      enddo


C
C     calculate the jacobian for the transformation:
C
      jacobi = 1.d0
      do i=2,n
         jacobi = jacobi * 2.d0*pi*p(i)
      enddo
      do i=3,n
         jacobi = jacobi * ( MM(i) - mu(i) )
      enddo
      jacobi = 1.d0/2.d0/dsqrt(s)*jacobi

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
      do rho=1,3
          pset(rho,1) = -pset(rho,2)
      enddo
      
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
