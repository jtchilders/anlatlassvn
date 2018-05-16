      subroutine phi1_2(x3,x4,p1,p2,p3,
     .     m2,m3,wt)
c     massive particle p1 decaying into p2 mass m2 and p3 mass m3.
c     with invariant mass 
c     of particle two s2 and particle three s3 integrated over.
c     vectors returned p2 and p3 are in the same frame as p1 is supplied
c     Expression evaluate is 
c     d^4 p2/(2 pi)^4   d^4 p3 /(2 pi)^4  (2 pi)^4 delta(p1-p2-p3)
c     (2 pi) delta(p2^2-m2^2) (2 pi) delta(p3^2-m3^2)
      implicit none
      include 'nlegborn.h'
      include '../include/pwhg_flst.h'
      include '../include/pwhg_kn.h'
      include '../include/pwhg_math.h'
      include 'PhysPars.h'
      real * 8  p1(4),p2(4),p3(4),p3cm(4)
      real * 8  x3,x4,costh,sinth,phi,cphi,sphi
      real * 8  wt,wt0,mod_p3
      real * 8  one,two
      real * 8  m1,s1,s2,s3,lambda,m2,m3
      integer j
      parameter(wt0=1d0/8d0/pi)
      parameter(one=1d0)
      parameter(two=2d0)

      s1=p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2  
      if (s1 .lt. 0d0) stop 'phi1_2: s1 < 0 ' 
      m1=sqrt(s1)

c     The following should never happen
      if(m1.lt.m2+m3) then 
         write(*,*) 'p1', p1 
         write(*,*) 'm1,m2,m3', m1, m2, m3
         write(*,*) 'phi1_2: m1 < m2+m3 ?'  
         call exit(-1)
      endif

      costh=two*x3-one      
      phi=two*pi*x4
      sinth=sqrt(one-costh**2)
      cphi=cos(phi)
      sphi=sin(phi)
      s2=m2**2
      s3=m3**2
      lambda=((s1-s2-s3)**2-4d0*s2*s3)
      if (lambda .lt. 0d0) then 
         write(*,*) 'lambda = ',lambda
c         write(*,*) s1,s2,s3
         write(*,*) '**** phi1_2: lambda < 0 *****' 
c     most probably is gonna be minus a small number...  turn it into a positive one
         lambda=abs(lambda)
      endif
      lambda=sqrt(lambda)
      wt=wt0*lambda/s1
      
      mod_p3 = lambda/(2*m1)

      p3cm(4)=m1/two*(s1+s3-s2)/s1
      p3cm(1)=mod_p3*sinth*sphi
      p3cm(2)=mod_p3*sinth*cphi
      p3cm(3)=mod_p3*costh
      call boost(m1,p1,p3cm,p3)
      do j=1,4
         p2(j)=p1(j)-p3(j)
      enddo

!      s34=p2(4)**2-p2(1)**2-p2(2)**2-p2(3)**2  
!      s56=p3(4)**2-p3(1)**2-p3(2)**2-p3(3)**2  
!      if (s34 .lt. 0d0) then
!         write(*,*)'s34 is lt zero in phi12',s34,x1,x2,x3,x4
!      elseif (s56 .lt. 0d0) then
!         write(*,*)'s56 is lt zero,in phi12',s56,x1,x2,x3,x4
!      endif
      if (  (p1(4) .lt. 0d0) 
     & .or. (p2(4) .lt. 0d0) 
     & .or. (p3(4) .lt. 0d0)) then
         write(*,*) 'phi1_2: one of E1,E2,E3 < 0'
         call exit(-1)
      endif

      end



