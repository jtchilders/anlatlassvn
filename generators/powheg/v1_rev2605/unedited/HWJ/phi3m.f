      subroutine phi3m(xth,xphi,p0,p1,p2,m1,m2,wt)
c     massive particle p0 in rest frame 
c     decaying into p1 fixed mass m1 and p2 fixed mass m2.
c     vectors returned p1 and p2 are in the frame in which 
C     p0 is supplied
c result is 1/8/pi * 2|p|/sqrts  * domega/(4*pi)
c     factor of (2*pi)^4 included in definition of phase space
      implicit none
      include 'pwhg_math.h'
      real * 8  p0(4),p1(4),p2(4),p1cm(4)
      real * 8  xth,xphi,phi,s,roots,costh,sinth
      real * 8  wt0,wt
      real * 8  m1,m2,m1sq,m2sq,lambda,lambda2,smin,mod_p1
      integer j
c      common/smin/smin
      parameter(wt0=1d0/8d0/pi)
      wt=0d0

      s=p0(4)**2-p0(1)**2-p0(2)**2-p0(3)**2  

      smin=(m1+m2)**2
      if (s .lt. smin) then
         write(*,*) 'phi3m: s < smin' 
         write(*,*) 'POWHEG ABORTS'
         call pwhg_exit(-1)
      endif

      roots=sqrt(s)
      m1sq=m1**2
      m2sq=m2**2
      costh=2d0*xth-1d0    
      sinth=sqrt(1d0-costh**2)
      phi=2d0*pi*xphi

      lambda2=((s+m1sq-m2sq)**2-4d0*m1sq*s)

      if (lambda2 .lt. 0d0) then
         write(6,*) 'phi3m:lambda2=', lambda2
         write(*,*) 'POWHEG ABORTS'
         call pwhg_exit(-1)
      endif
      lambda=sqrt(lambda2)

      wt=wt0*lambda/s

      mod_p1 = lambda/(2*roots)

      p1cm(4)=roots/2d0*(s+m1sq-m2sq)/s
      p1cm(1)=mod_p1*sinth*sin(phi)
      p1cm(2)=mod_p1*sinth*cos(phi)
      p1cm(3)=mod_p1*costh

c      write(6,*) 'e',roots/2d0*(s+m1sq-m2sq)/s
c      write(6,*) 'p',roots/2d0*lambda/s

c      write(6,*) 'sinth',sinth
c      write(6,*) 'costh',costh
c      write(6,*) 'p1cm**2',p1cm(4)**2-p1cm(1)**2-p1cm(2)**2-p1cm(3)**2
c      pause

      call boost(roots,p0,p1cm,p1)
      do j=1,4
      p2(j)=p0(j)-p1(j)
      enddo

      
      if (  (p0(4) .lt. 0d0) 
     &     .or. (p1(4) .lt. 0d0) 
     &     .or. (p2(4) .lt. 0d0)) then  
         write(6,*) 'p0',p0(4),p0(4)**2-p0(1)**2-p0(2)**2-p0(3)**2,s
         write(6,*) 'p1',p1(4),p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2,m1sq
         write(6,*) 'p2',p2(4),p2(4)**2-p2(1)**2-p2(2)**2-p2(3)**2,m2sq
         write(6,*) 'in phi3m'
      endif
      end
      
