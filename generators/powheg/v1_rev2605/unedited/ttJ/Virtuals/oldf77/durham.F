C     $Modified: Thu Jun  3 18:15:53 2004 by puwer $
      function durham4(rs,ycut,p4,p3)
      implicit none
      integer nparton
      parameter (nparton=4)
      integer durham4,durham3
      double precision ycut,rymin,rs
      double precision p4(4,-2:nparton),p3(4,-2:nparton-1),
     &     p2(4,-2:nparton-2)
      double precision r(nparton,nparton)
      double precision yij,rcostheta
      integer i,j,k,ihelp,jhelp,mu      
      double precision dotp,rmin
      rymin =1.d99

      do i=1,nparton
         do j=i,nparton
            r(i,j) = dotp(p4(1,i),p4(1,j))
            r(j,i) = r(i,j)
         enddo
      enddo

c     use the DURHAM cut, y_ij = 2 * min(E_i^2 ,E_j^2) * ( 1 - costheta_ij )
      DO i=1,nparton-1
         DO j=i+1,nparton
            rcostheta = 1.d0 /
     &           dsqrt( p4(4,i)**2  - r(i,i)  ) /
     &           dsqrt( p4(4,j)**2  - r(j,j)  ) *
     &           (p4(4,i) * p4(4,j)  - r(i,j) )
            yij = 2.d0 / rs * rmin( p4(4,i)**2, p4(4,j)**2) *
     &           ( 1.d0 - rcostheta )
            IF ( yij.lt.rymin ) THEN
               rymin = yij
               ihelp = i
               jhelp = j
            ENDIF
         ENDDO
      ENDDO

      if (rymin.gt.ycut) then
C     4-jet event
         durham4 = 4
         return
      endif
C
C     combine the momenta
C
      do mu=1,4
         p3(mu,-1) = p4(mu,-1)
         p3(mu,-2) = p4(mu,-2)
      enddo
      k=1
      do i = 1,4
         if ((i.ne.ihelp).and.(i.ne.jhelp)) then
            do mu=1,4
               p3(mu,k) = p4(mu,i)
            enddo
            k = k+1
         endif
      enddo
      do mu=1,4
         p3(mu,3) = p4(mu,ihelp) + p4(mu,jhelp)
      enddo
      durham4 = durham3(rs,ycut,p3,p2)
      end
C***********************************************************************
      function durham3(rs,ycut,p3,p2)
      implicit none
      integer durham3
      double precision ycut,rymin
      integer nparton
      parameter (nparton=3)
      double precision p3(4,-2:nparton),rs,p2(4,-2:nparton-1),
     .     r(nparton,nparton)
      double precision yij,rcostheta
      integer i,j,k,ihelp,jhelp,mu
      double precision dotp,rmin

      rymin = 1.d99
      do i=1,nparton
         do j=i,nparton
            r(i,j) = dotp(p3(1,i),p3(1,j))
            r(j,i) = r(i,j)
         enddo
      enddo
c     use the DURHAM cut, y_ij = 2 * min(E_i^2 ,E_j^2) * ( 1 - costheta_ij )
      DO i=1,nparton-1
         DO j=i+1,nparton
            rcostheta = 1.d0 /
     &           dsqrt( p3(4,i)**2  - r(i,i)  ) /
     &           dsqrt( p3(4,j)**2  - r(j,j)  ) *
     &           (p3(4,i) * p3(4,j)  - r(i,j) )
            yij= 2.d0 / rs * rmin( p3(4,i)**2, p3(4,j)**2) *
     &           ( 1.d0 - rcostheta )
            IF ( yij.lt.rymin )  THEN
               rymin = yij
               ihelp = i
               jhelp = j
            ENDIF
         ENDDO
      ENDDO
C      write(*,*) 'durham3: ', rymin 
      if (rymin.gt.ycut) then
C     4-jet event
         durham3 = 3
         return
      else
         durham3 = 2
         return
      endif
C     combine the momenta
      do mu=1,4
         p2(mu,-1) = p3(mu,-1)
         p2(mu,-2) = p3(mu,-2)
      enddo
      k=1
      do i = 1,nparton
         if ((i.ne.ihelp).and.(j.ne.ihelp)) then
            do mu=1,4
               p2(mu,k) = p3(mu,i)
            enddo
            k = k+1
         endif
      enddo
      do mu=1,4
         p2(mu,2) = p3(mu,ihelp)+p3(mu,jhelp)
      enddo
      end

C***********************************************************************
      function durham2(rs,ycut,p2,p1)
      implicit none
      integer durham2
      double precision ycut,rymin
      integer nparton
      parameter (nparton=2)
      double precision p2(4,-2:nparton),rs,p1(4),r(nparton,nparton)
      double precision yij,rcostheta
      integer i,j,k,ihelp,jhelp,mu
      double precision dotp,rmin
      rymin = 1.d99

      do i=1,nparton
         do j=i,nparton
            r(i,j) = dotp(p2(1,i),p2(1,j))
            r(j,i) = r(i,j)
         enddo
      enddo
c     use the DURHAM cut, y_ij = 2 * min(E_i^2 ,E_j^2) * ( 1 - costheta_ij )
      DO i=1,nparton-1
         DO j=i+1,nparton
            rcostheta = 1.d0 /
     &           dsqrt( p2(4,i)**2  - r(i,i)  ) /
     &           dsqrt( p2(4,j)**2  - r(j,j)  ) *
     &           (p2(4,i) * p2(4,j)  - r(i,j) )
            yij= 2.d0 / rs * rmin( p2(4,i)**2, p2(4,j)**2) *
     &           ( 1.d0 - rcostheta )
            IF ( yij.lt.rymin )  THEN
               rymin = yij
               ihelp = i
               jhelp = j
            ENDIF
         ENDDO
      ENDDO
C      write(*,*) 'durham2: ', rymin 
      if (rymin.gt.ycut) then
C     2-jet event
         durham2 = 2
         return
      else
         durham2 = 1
C     combine the momenta
         do mu=1,4
            p1(mu) = p2(mu,1)+p2(mu,2)
         enddo
      endif
      end
C***********************************************************************      
      function rmin(a,b)
      implicit none
      double precision rmin,a,b
      if (a.lt.b) rmin = a
      if (b.lt.a) rmin = b
      end
C***********************************************************************      
      function rmax(a,b)
      implicit none
      double precision rmax,a,b
      if (b.gt.a) rmax = b
      if (a.gt.b) rmax = a
      end
