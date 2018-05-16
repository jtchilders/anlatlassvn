C
c***************************************************************************
      function dotrr(p1,p2)
c***************************************************************************
c     dotrr(p1,p2) = p1.p2
c***************************************************************************

      implicit none
      double precision dotrr,p1(0:3),p2(0:3)

      dotrr = p1(0)*p2(0) - p1(1)*p2(1) - p1(2)*p2(2) - p1(3)*p2(3)

      end

c***************************************************************************
      double complex function dotcc(v1,v2)
c***************************************************************************

      implicit none
      double complex v1(0:3), v2(0:3)

      dotcc = v1(0)*v2(0)-v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3)

      end

c***************************************************************************
      double complex function dotcr(v2,v1)
c***************************************************************************

      implicit none
      double precision v1(0:3)
      double complex  v2(0:3)
   
      dotcr = v1(0)*v2(0)-v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3)
   
      end

c***************************************************************************
      double complex function dotrc(v1,v2)
c***************************************************************************

      implicit none
      double precision v1(0:3)
      double complex  v2(0:3)
   
      dotrc = v1(0)*v2(0)-v1(1)*v2(1)-v1(2)*v2(2)-v1(3)*v2(3)
   
      end

c***************************************************************************
      double complex function contract_Tjj(T,j1,j2)
c***************************************************************************

      implicit none
c contract complex rank 2 tensor T^{mu,nu} with two complex vectors j1 and j2
      complex*16 T(0:3,0:3), j1(0:3), j2(0:3), resv(0:3)
      integer mu
      do mu = 0,3
         resv(mu) = T(mu,0)*j2(0) - T(mu,1)*j2(1) 
     &            - T(mu,2)*j2(2) - T(mu,3)*j2(3)
      enddo
      contract_Tjj = resv(0)*j1(0) - resv(1)*j1(1) 
     &             - resv(2)*j1(2) - resv(3)*j1(3)
      return
      end
      
C-------------------------------------------------------------------

      FUNCTION SC3(CHII,A1,A2,A3,CHIF,ALPHA)
      IMPLICIT NONE 
      COMPLEX*16 ZI, SC3
      PARAMETER ( ZI=(0D0,1D0) )
      INTEGER  ALPHA, ALP, N, I
      COMPLEX*16  CHII(2), CHIF(2), CHIAUX(2), CHIDUM, ASLASH(2,2)
      COMPLEX*16  A1(0:3), A3(0:3), AUX(0:3,3)
      REAL*8  A2(0:3)
C
      N = 3
      DO I = 0,3
         AUX(I,3) = A3(I)
         AUX(I,2) = A2(I)
         AUX(I,1) = A1(I)
      ENDDO
      CHIAUX(1) = CHII(1)
      CHIAUX(2) = CHII(2)
      ALP = ALPHA
C
      DO 30 I = 1,N
         IF (ALP.GT.0) THEN
            ASLASH(1,1) = AUX(0,I) - AUX(3,I)
            ASLASH(1,2) = -AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,1) = -AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) + AUX(3,I)
         ELSE
            ASLASH(1,1) = AUX(0,I) + AUX(3,I)
            ASLASH(1,2) = AUX(1,I) - ZI*AUX(2,I)
            ASLASH(2,1) = AUX(1,I) + ZI*AUX(2,I)
            ASLASH(2,2) = AUX(0,I) - AUX(3,I)
         ENDIF
         CHIDUM = CHIAUX(1)*ASLASH(1,1) + CHIAUX(2)*ASLASH(2,1)
         CHIAUX(2) = CHIAUX(1)*ASLASH(1,2) + CHIAUX(2)*ASLASH(2,2)
         CHIAUX(1) = CHIDUM
         ALP = - ALP
 30   CONTINUE
C
      SC3  = CHIAUX(1) * CHIF(1) + CHIAUX(2) * CHIF(2)
      RETURN
      END

    
c***************************************************************************
      double complex function contract_Tjjj(T,j1,j2,j3)
c***************************************************************************

      implicit none
c contract complex rank 2 tensor T^{mu,nu,ka} with 
c	three complex vectors j1, j2, j3
      complex*16 T(0:3,0:3,0:3), j1(0:3), j2(0:3), j3(0:3)
      complex*16 resv(0:3),resvv(0:3,0:3)
      integer mu,nu
      do mu = 0,3
         do nu = 0,3
            resvv(mu,nu) = T(mu,nu,0)*j3(0) - T(mu,nu,1)*j3(1) 
     &                   - T(mu,nu,2)*j3(2) - T(mu,nu,3)*j3(3)
         enddo
      enddo 
      do mu = 0,3
         resv(mu) = resvv(mu,0)*j2(0) - resvv(mu,1)*j2(1) 
     &            - resvv(mu,2)*j2(2) - resvv(mu,3)*j2(3)
      enddo
      contract_Tjjj = resv(0)*j1(0) - resv(1)*j1(1) 
     &              - resv(2)*j1(2) - resv(3)*j1(3)
      return
      end
      
c***************************************************************************
      subroutine vcartx(p,vmass,vwidth,ncart,nsv , vc)
c***************************************************************************
c
c This subroutine computes an effective VECTOR wavefunction for an internal
c vector boson line. The propagator is inlcuded in the Feynman gauge.
c
c input:
c       real    p(0:3)         : four-momentum of vector boson
c       real    vmass          : mass          of vector boson
c       integer ncart = 0,1,2,3: cartesian polarization direction
c                                of vector boson
c       integer nsv  = -1 or 1 : +1 for final, -1 for initial
c
c output:
c       complex vc(6)          : vector wavefunction       epsilon^mu(v)
c     
      implicit none
      double complex vc(6), d
      double precision p(0:3),vmass,vwidth,q2
      integer ncart, mu, nsv

      q2 = p(0)**2-p(1)**2-p(2)**2-p(3)**2
      if (vmass.eq.0d0) then
         d = 1d0/q2
      else
         d = 1d0/dcmplx( q2-vmass**2, vmass*vwidth )
      endif
      do mu = 0,3
         vc(mu+1) = dcmplx(0d0,0d0)
      enddo
      if (ncart.eq.0) then
         vc(1) = d
      else
         vc(ncart+1) = -d
      endif
c
      vc(5) = dcmplx(p(0),p(3))*nsv
      vc(6) = dcmplx(p(1),p(2))*nsv
c
      return
      end
      
