C$Modified: Sun Apr 18 22:22:32 2004 by puwer $
C*************************************************************************
      subroutine printConf(p)
      implicit none
      double precision p(4,-2:3)
      integer mu,i
      do i=-2,-1
         write(*,*) '---- part.: ',i,'--------'
         do mu = 1,4
            write(*,*) p(mu,i)
         enddo
      enddo  
      do i=1,3
         write(*,*) '---- part.: ',i,'--------'
         do mu = 1,4
            write(*,*) p(mu,i)
         enddo
      enddo
      end
C*************************************************************************
      subroutine CheckKinematics(p)
      implicit none
      double precision p(4,-2:4)
      double precision kin1(4),kin2(4),kq(4),kqb(4),k1(4),k2(4)
      double precision ktot1(4),ktot2(4)
      double precision rcostheta,yij
      double precision dotp,rmin
      integer mu,nparton,i,j
      parameter(nparton=4)
#include "common.F"
      do mu =1,4
         kin1(mu) = p(mu,-1)
         kin2(mu) = p(mu,-2)
         kq(mu) = p(mu,1)
         kqb(mu)= p(mu,2)
         k1(mu) = p(mu,3)
         k2(mu) = p(mu,4)
         ktot1(mu) = p(mu,-1) + p(mu,-2)
         ktot2(mu) = p(mu,1) + p(mu,2)+p(mu,3) + p(mu,4)
      enddo

      write(*,*) '>>>>>>>> CheckKinematics >>>>>>>'
      write(*,*) ' massless momenta: '
      write(*,*) dotp(kin1,kin1)
      write(*,*) dotp(kin2,kin2)
      write(*,*) dotp(k1,k1)
      write(*,*) dotp(k2,k2)
      write(*,*) ' massiv momenta: '
      write(*,*) dsqrt(dotp(kq,kq))
      write(*,*) dsqrt(dotp(kqb,kqb))
      write(*,*) 'Momentum conservation'
      do mu=1,4
         write(*,*) kin1(mu) + kin2(mu) ,' = ',
     >        k1(mu)+k2(mu)+kq(mu)+kqb(mu)
      enddo
      write(*,*) 'Center of mass enery squared'
      write(*,*)  dotp(ktot1,ktot1),' = ',dotp(ktot2,ktot2) 
      write(*,*) 'dot products:'
      write(*,*) 'kin1.k3 = ', dotp(kin1,k1)
      write(*,*) 'kin2.k3 = ', dotp(kin2,k1)
      write(*,*) 'kin1.k4 = ', dotp(kin1,k2)
      write(*,*) 'kin2.k4 = ', dotp(kin2,k2)
      write(*,*) 'kq.k3 = ', dotp(kq,k1)
      write(*,*) 'kq.k4 = ', dotp(kq,k2)
      write(*,*) 'kqb.k3 = ', dotp(kqb,k1)
      write(*,*) 'kqb.k4 = ', dotp(kqb,k2)
      write(*,*) 'k3.k4 = ', dotp(k1,k2)
      DO i=1,nparton-1
         DO j=i+1,nparton
            write(*,*) i,j
            rcostheta = 1.d0 /
     &           dsqrt( p(4,i)**2  - dotp(p(1,i),p(1,i))  ) /
     &           dsqrt( p(4,j)**2  - dotp(p(1,j),p(1,j))  ) *
     &           (p(4,i) * p(4,j)  - dotp(p(1,i),p(1,j)) )
            yij = 2.d0 / shad * rmin( p(4,i)**2, p(4,j)**2) *
     &           ( 1.d0 - rcostheta )
            write(*,*) 'yij=',yij
         ENDDO
      ENDDO

      write(*,*) '<<<<<<<< CheckKinematics <<<<<<<'
      end
