C$Modified: Thu Apr  1 11:30:05 2004 by puwer $
      subroutine CollConf(i,j,ktq,z,pout)
      implicit none
      integer i,j,mu
      double precision ktq
      double precision pout(4,-2:4) 
      double precision delta,z
      double precision p(4,100) ! momenta
      double precision w        ! weight 
      double precision masses(100) ! masses
      double precision psum(4),pm(4),pi(4),pj(4)
      double precision k34(4), k3(4),k4(4),pvec(4),n(4),kt(4),kg(4)
      double precision kth(4),kgh(4),ksum(4)
      double precision p1(4),p2(4),kq(4),kqb(4),psum1(4),psum2(4)
      double precision beta,n1(4),n2(4),n3(4)
      double precision mt,sqrts, z1,z2,ndotp,help ,xdoof
      double precision theta, phi,xhelp
      double precision rpi,drandm, drnum,rn
      double precision dotp
      double complex spa, spb
      integer ihelp
      integer*4 iflag
      external rn
      rpi = dacos(-1.d0)
      sqrts = 1500.d0
      mt = 175.d0   
      
      if (j.lt.i) then
         ihelp = j
         j = i
         i = ihelp
      endif
      if (ktq.ge.0.d0) then
         write(*,*) 'ktq should be negative'
         stop
      endif


      if (i.lt.0) then 
C     create the outgoing momenta:
         masses(1) = mt
         masses(2) = mt
         masses(3) = 0.d0
         masses(4) = 0.d0
         call rambo(4, sqrts, masses,p, w)
         do mu=1,4
            pout(mu,1) = p(mu,1)
            pout(mu,2) = p(mu,2)
            pout(mu,3) = p(mu,3)
            pout(mu,4) = p(mu,4)
         enddo
C     
         if ((j.eq.1).or.(j.eq.2)) then
            beta = dsqrt(1.d0 - mt**2/pout(4,j))
         else
            beta =1.d0
         endif
         do mu=1,3
            n3(mu) =  pout(mu,j)/pout(4,j)/beta
         enddo
         n1(1) =  n3(2)/dsqrt(n3(1)**2+n3(2)**2)
         n1(2) = -n3(1)/dsqrt(n3(1)**2+n3(2)**2)
         n1(3) = 0.d0
         n2(1) = n3(2)*n1(3)-n3(3)*n1(2)
         n2(2) = n3(3)*n1(1)-n3(1)*n1(3)
         n2(3) = n3(1)*n1(2)-n3(2)*n1(1)
         theta = dacos(1.d0  + ktq/pout(4,j)/sqrts*2.d0)
         phi = 2.d0 * rpi * rn(4711.d0)
         pout(4,-1) = sqrts/2.d0
         pout(4,-2) = sqrts/2.d0
         if (i.eq.-1) then
            do mu = 1,3
               pout(mu,-1) = sqrts/2.d0 * (
     .              dsin(theta)*dcos(phi) * n1(mu)
     .              + dsin(theta)*dsin(phi) * n2(mu)
     .              + dcos(theta)*n3(mU)
     .              )
               pout(mu,-2) = - pout(mu,-1)
            enddo
         else
            do mu = 1,3
               pout(mu,-2) = sqrts/2.d0 * (
     .          dsin(theta)*dcos(phi) * n1(mu)
     .              + dsin(theta)*dsin(phi) * n2(mu)
     .              + dcos(theta)*n3(mU)
     .              )
               pout(mu,-1) = - pout(mu,-2)
            enddo
         endif
#ifdef DEBUG 
      do mu=1,4
         p1(mu) = pout(mu,-1)  
         p2(mu) = pout(mu,-2) 
         kq(mu) = pout(mu,1)  
         kqb(mu) = pout(mu,2)
         k3(mu) = pout(mu,3)
         k4(mu) = pout(mu,4)
         psum1(mu) = p1(mu) + p2(mu)
         psum2(mu) = kq(mu) + kqb(mu) + k3(mu) + k4(mu)
      enddo
      write(*,*) 's = ',dotp(psum1,psum1) 
      write(*,*) 's = ',dotp(psum2,psum2) 
      write(*,*) 'ktq = ', ktq
      write(*,*) 'dotp(p1,k3) =',dotp(p1,k3)
      write(*,*) 'dotp(p2,k3) =',dotp(p2,k3)
      write(*,*) 'dotp(p1,k4) =',dotp(p1,k4)
      write(*,*) 'dotp(p2,k4) =',dotp(p2,k4)
#endif      
      return
      endif

      if ((i.eq.3).and.(j.eq.4)) then
C
C     create first the two incoming momenta:
         masses(1) = 0.d0
         masses(2) = 0.d0
         call rambo(2, sqrts, masses,p, w)
         do mu=1,4
            pout(mu,-1) = p(mu,1)
            pout(mu,-2) = p(mu,2)
         enddo
         
         if (z.le.0) then
            z1 = rn(xdoof)
            xdoof = z1
            z2 = 1.d0 - z1
         else
            if (z.lt.1.d0) then
               z1 = z
               z2 = 1.d0 - z1
            else
               write(*,*) 'CollConf: z not in allowed range'
               stop
            endif
         endif
         masses(1) = mt
         masses(2) = mt
         masses(3) = dsqrt(-ktq/z1/z2)
         call rambo(3, sqrts, masses,p, w)
         do mu=1,4
            pout(mu,1) = p(mu,1)
            pout(mu,2) = p(mu,2)
            k34(mu) = p(mu,3)
         enddo
C     generate n:
         theta = 3.14159d0 * rn(xdoof)
         xdoof = theta
         phi   = 2.d0 * 3.14159d0 * rn(xdoof)
         xdoof = phi
         n(1) = dsin(theta)*dcos(phi)
         n(2) = dsin(theta)*dsin(phi)
         n(3) = dcos(theta)
         n(4) = 1.d0
         ndotp = dotp(n,k34)
         do mu = 1,4
            pvec(mu) = k34(mu) + ktq/2.d0/ndotp*n(mu)/z1/z2
         enddo
         kt(4) = 0.d0
         kt(1) = pvec(2)*n(3) - pvec(3)*n(2)
         kt(2) = pvec(3)*n(1) - pvec(1)*n(3)
         kt(3) = pvec(1)*n(2) - pvec(2)*n(1)
         xdoof = dsqrt(kt(1)**2 + kt(2)**2 + kt(3)**2)
         do mu=1,3
            kt(mu) = kt(mu)/xdoof*dsqrt(-ktq)
         enddo
         do mu = 1,4
            pout(mu,3) = z1*pvec(mu) + kt(mu) - ktq/z1*n(mu)/2.d0/ndotp
            pout(mu,4) = z2*pvec(mu) - kt(mu) - ktq/z2*n(mu)/2.d0/ndotp
            k3(mu) = pout(mu,3)
            k4(mu) = pout(mu,4)
         enddo
#ifdef DEBUG
         write(*,*) 'preset ktq: ',ktq
         write(*,*) '     kt.kt: ',dotp(kt,kt)
         write(*,*) '       p.p: ',dotp(p,p)
         write(*,*) '       n.n: ',dotp(n,n)
         write(*,*) '      kt.n: ',dotp(kt,n)
         write(*,*) '      kt.p: ',dotp(kt,pvec)
         write(*,*) '   2 k3.k4: ',2.d0 * dotp(k3,k4)
         write(*,*) '     k34^2: ',dotp(k34,k34)
         write(*,*) '-ktq/z1/z2: ', -ktq/z1/z2
         write(*,*) '      k3^2: ',dotp(k3,k3)
         write(*,*) '      k4^2: ',dotp(k4,k4)
         write(*,*) '       n.p: ',dotp(n,pvec)
         write(*,*) '     k34.n: ',dotp(n,k34)
         do mu=1,4
            write(*,*) k34(mu), ' = ', k3(mu) + k4(mu) 
         enddo
         do mu=1,4
            write(*,*) pout(mu,-1)+pout(mu,-2), ' = ', 
     >          pout(mu,1) + pout(mu,2) + pout(mu,3) + pout(mu,4) 
         enddo

#endif      
      return
      endif
      if ((i.lt.3).and.(j.gt.2)) then
C
C     create first the two incoming momenta:
         masses(1) = 0.d0
         masses(2) = 0.d0
         call rambo(2, sqrts, masses,p, w)
         do mu=1,4
            pout(mu,-1) = p(mu,1)
            pout(mu,-2) = p(mu,2)
         enddo
C
C     
         if (z.le.0) then
            z1 = rn(xdoof)
            xdoof = z1
            z2 = 1.d0 - z1
         else
            if (z.lt.1.d0) then
               z1 = z
               z2 = 1.d0 - z1
            else
               write(*,*) 'CollConf: z not in allowed range'
               stop
            endif
         endif
         masses(1) = mt
         masses(2) = dsqrt(-ktq/z1/z2+mt**2/z1)
         masses(3) = 0
         write(*,*) ktq,1.d0/z1/z2,masses(2)
         call rambo(3, sqrts, masses,p, w)
         if (i.eq.1) then
            do mu=1,4
               pout(mu,2) = p(mu,1)
            enddo
         else
            do mu=1,4
               pout(mu,1) = p(mu,1)
            enddo
         endif
         if (j.eq.3) then
            do mu=1,4
               pout(mu,4) = p(mu,3)
            enddo
         else
            do mu=1,4
               pout(mu,3) = p(mu,3)
            enddo
         endif
         do mu=1,4
            ksum(mu) = p(mu,2)
         enddo
C     generate n:
         theta = 3.14159d0 * rn(xdoof)
         xdoof = theta
         phi   = 2.d0 * 3.14159d0 * rn(xdoof)
         xdoof = phi
         n(1) = dsin(theta)*dcos(phi)
         n(2) = dsin(theta)*dsin(phi)
         n(3) = dcos(theta)
         n(4) = 1.d0
         ndotp = dotp(n,ksum)
         do mu = 1,4
            pvec(mu) = ksum(mu) 
     .           + n(mu)/2.d0/ndotp*(
     .           (ktq+z1**2*mt**2-mt**2)/z1
     .           + (ktq+z2**2*mt**2)/z2 )
         enddo
         kt(4) = 0.d0
         kt(1) = pvec(2)*n(3) - pvec(3)*n(2)
         kt(2) = pvec(3)*n(1) - pvec(1)*n(3)
         kt(3) = pvec(1)*n(2) - pvec(2)*n(1)
         xdoof = dsqrt(kt(1)**2 + kt(2)**2 + kt(3)**2)
         do mu=1,3
            kt(mu) = kt(mu)/xdoof*dsqrt(-ktq)
         enddo
         do mu = 1,4
            kq(mu) = z1*pvec(mu) + kt(mu) 
     .           - (ktq+z1**2*mt**2-mt**2)/z1*n(mu)/2.d0/ndotp
            kg(mu) = z2*pvec(mu) - kt(mu) 
     .           - (ktq+z2**2*mt**2)/z2 *n(mu)/2.d0/ndotp
         enddo
         if (i.eq.1) then
            do mu=1,4
               pout(mu,1) = kq(mu)
            enddo
         else
            do mu=1,4
               pout(mu,2) = kq(mu)
            enddo
         endif
         if (j.eq.3) then
            do mu=1,4
               pout(mu,3) = kg(mu)
            enddo
         else
            do mu=1,4
               pout(mu,4) = kg(mu)
            enddo
         endif
         write(*,*) 'preset ktq: ',ktq
         write(*,*) '     kt.kt: ',dotp(kt,kt)
         write(*,*) '       p.p: ',dotp(pvec,pvec)
         write(*,*) '       n.n: ',dotp(n,n)
         write(*,*) '      kt.n: ',dotp(kt,n)
         write(*,*) '      kt.p: ',dotp(kt,pvec)
         write(*,*) dotp(kq,kq) + 2.d0 * dotp(kq,kg), ' = ',
     .        -ktq/z1/z2+mt**2/z1
      endif

      end
