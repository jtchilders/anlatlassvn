      subroutine phi1_2m_nobw_ktmin(p1,m2,m3,ktmin,xth,xphi,
     $     p2,p3,wt)
c     massive particle p1 decaying into p2 mass m2 and p3, mass m3. 
c     Minimum transverse momentum of p2 (or p3) greater than ktmin.
c     Vectors returned p2 and p3 are in the same frame as p1 is supplied
c     Expression evaluate is 
c     ds2 d^4 p2 d^4 p3 (2 pi)^4 delta(p1-p2-p3)/(2 pi)^6
c     delta(p2^2-s2) delta(p3^2-s3)
      implicit none
      include 'pwhg_math.h'
      double precision p1(4),p2(4),p3(4),p3cm(4)
      double precision x3,xth,xphi,costh,sinth,phi,cphi,sphi,ktmin
      double precision wt,wt0,w3
      double precision s3max,s3min
      double precision m1,m2,m3,s1,s2,s3,lambda,lambda2
      integer j
      parameter(wt0=1d0/8d0/pi/(4*pi))
      real * 8 z,t,mod_p2,sinth_min,th_min,max_costh

      wt=0d0      
      w3=1d0
      s3=m3**2
      s2=m2**2
      s1=p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2  
      if (s1 .lt. 0d0) then
         write(*,*) '** error in phi1_2m_nobw_ktmin: s1 < 0 **'
         write(*,*) 'The POWHEG BOX continues but the problem '
     $        //'should be solved'
         return
      endif
      m1=dsqrt(s1)

      lambda2=((s1-s2-s3)**2-4d0*s2*s3)
      
      if (lambda2.le.0) then
         write(*,*) "** PROBLEMS in phi1_2m_nobw_ktmin: lambda2 < 0 **"
         write(*,*) 'The POWHEG BOX continues but the problem '
     $        //'should be solved'
         wt=0d0 
         return
      endif
c     mod_p2 is the modulus of the vector p2 or p3, one recoiling against the other
      lambda=sqrt(lambda2)
      mod_p2=lambda/(2*m1)      
      sinth_min = ktmin/mod_p2
      if (sinth_min.gt.1) then
         write(*,*) "** PROBLEMS in phi1_2m_nobw_ktmin: sinth_min>1 **"
         write(*,*) 'The POWHEG BOX continues but the problem '
     $        //'should be solved'
         return
      endif
      th_min = asin(sinth_min)
      max_costh = cos(th_min)

      t = 2d0*xth-1d0      
      w3=w3*2
      z=1.5d0*(t-t**3/3)
      w3 = w3*1.5d0*(1-t**2)
c      costh=2*max_costh*xth-max_costh      
      costh=max_costh*z
      w3=w3*max_costh

      phi=2d0*pi*xphi
      w3=w3*2d0*pi

      sinth=sqrt(1d0-costh**2)
      cphi=cos(phi)
      sphi=sin(phi)

      wt=wt0*w3*lambda/s1

      p3cm(4)=m1/2d0*(s1+s3-s2)/s1

c      p3cm(1)=m1/2d0*lambda/s1*sinth*sphi
c      p3cm(2)=m1/2d0*lambda/s1*sinth*cphi
c      p3cm(3)=m1/2d0*lambda/s1*costh

      p3cm(1)=mod_p2*sinth*sphi
      p3cm(2)=mod_p2*sinth*cphi
      p3cm(3)=mod_p2*costh

      call boost(m1,p1,p3cm,p3)
      do j=1,4
      p2(j)=p1(j)-p3(j)
      enddo
      if (  (p1(4) .lt. 0d0) 
     &     .or. (p2(4) .lt. 0d0) 
     &     .or. (p3(4) .lt. 0d0)) then  
         write(6,*) 'p1(4)',p1(4)
         write(6,*) 'p2(4)',p2(4)
         write(6,*) 'p3(4)',p3(4)
         write(6,*) 'p1sq',p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2,s1
         write(6,*) 'p2sq',p2(4)**2-p2(1)**2-p2(2)**2-p2(3)**2,s2
         write(6,*) 'p3sq',p3(4)**2-p3(1)**2-p3(2)**2-p3(3)**2,s3
         write(6,*) 'in phi1_2m.f'
      endif
      end



