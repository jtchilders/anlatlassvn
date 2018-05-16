C
C     $Modified: Mon Apr  6 14:26:45 2009 by uwer $
C
C     some functions to calculate spinor products
C
#define WEINZIERL
#define STEFAN
      function f77spaversion()
      implicit none
      integer f77spaversion
#ifdef STEFAN
      f77spaversion = 0
#else
      f77spaversion = 1
#endif
      end
#ifdef DIXON
      function spa(p1,p2)
      implicit none
      double complex spa
      double complex spb
      double precision p1(4),p2(4)
      spa = conjg(spb(p2,p1))
      end

      function spb(p1,p2)
      implicit none
      double complex spb
      double precision p1(4),p2(4),dotp
      double precision sij,p1p,p2p,cphi,sphi
C     Formualae taken from Lance Dixon hep-ph/9601359
      sij = 2.d0 * dotp(p1,p2)
      p1p = p1(4) + p1(3)
      p2p = p2(4) + p2(3)
      cphi = ( p1(1)*p2p-p2(1)*p1p ) / dsqrt( p1p*p2p* dabs(sij) )
      sphi = ( p1(2)*p2p-p2(2)*p1p ) / dsqrt( p1p*p2p* dabs(sij) )
      spb = -dsqrt( dabs(sij) ) * (cphi *(1,0) - (0,1) * sphi)
      return
      end
#endif
#ifdef HLPSW
C taken from Harris, Laenen,Phaf,Sullivan,Weinzierl
      function spa(p1,p2)
      implicit none
      double complex spa
      double precision p1(4),p2(4)
      double precision pp,qp
      double complex pt,qt
      if ((p1(4).gt.0.d0) .and. (p2(4).gt.0.d0) ) then
         pp = p1(4) + p1(3)
         pt = p1(1) + (0.d0,1.d0) * p1(2)
         
         qp = p2(4) + p2(3)
         qt = p2(1) + (0.d0,1.d0) * p2(2)
         spa = 1.d0 / dsqrt(pp*qp) * (qp*pt - pp*qt)
      else
         stop
      endif
      end
C      
      function spb(p1,p2)
      implicit none
      double precision p1(4),p2(4)
      double complex spa
      double complex spb
      double precision dotp,sign
      if (dotp(p1,p2).ge.0) then
         sign = 1.d0
      else
         sign = -1.d0
      endif
      
      spb = conjg(spa(p2,p1)) * sign
      end
#endif
#ifdef WEINZIERL
      function spa(p1,p2)
      implicit none
      double complex spa
      integer n
      double precision p1(4),p2(4)
      double precision pplus,qplus,x,y,prefac
      double complex pt,qt,a,b
C     Formualae taken from  Stefans thesis

#ifdef STEFAN
      pplus = p1(4) + p1(3)
      qplus = p2(4) + p2(3)
#else
      pplus = p1(4) + p1(1)
      qplus = p2(4) + p2(1)
#endif
      if (pplus.eq.0.d0) write(*,*) 'p1= ',p1
      if (qplus.eq.0.d0) write(*,*) 'p2= ',p2

C     Original Version:
**      pt = p1(1) + (0.d0,1.d0) * p1(2)
**      qt = p2(1) + (0.d0,1.d0) * p2(2)
**      a = qplus*pt
**      b = pplus*qt
      prefac = 1.d0 / dsqrt(dabs(pplus*qplus))
#ifdef STEFAN
      x = ( qplus * p1(1) - pplus * p2(1) ) * prefac;
      y = ( qplus * p1(2) - pplus * p2(2) ) * prefac;
#else
      x = ( qplus * p1(2) - pplus * p2(2) ) * prefac;
      y = ( qplus * p1(3) - pplus * p2(3) ) * prefac;
#endif
#ifdef ORIGINALVERSION
      a = p1(1)*qplus + (0.d0,1.d0) * qplus * p1(2)
      b = p2(1)*pplus + (0.d0,1.d0) * pplus * p2(2)
#endif
C     
      n = 0
      if (p1(4).lt.0.d0) n = 1
      if (p2(4).lt.0.d0) n = n + 1

      if (n.eq.0) spa = x + y * (0,1.d0)
      if (n.eq.1) spa = y - x * (0,1.d0)
      if (n.eq.2) spa = -x - y * (0,1.d0)

#ifdef ORIGINALVERSION
C      spa = (0.d0,-1.d0)**n / dsqrt(dabs(pplus*qplus)) * (a - b)
      spa = 1.d0 / dsqrt(dabs(pplus*qplus)) * (a - b)
      if (n.eq.0) spa = spa 
      if (n.eq.1) spa = (0,-1.d0) * spa
      if (n.eq.2) spa = -1.d0 * spa
#endif
      end
C      
      function spb(p1,p2)
      implicit none
      double precision p1(4),p2(4)
      double complex spa
      double complex spb
      double precision sign
**      if (dotp(p1,p2).ge.0) then
      spb = conjg(spa(p2,p1))
      if ((p1(4)*p2(4)).lt.0.d0) spb = spb * (-1.d0)
      end
#endif

      function dotp(p1,p2)
      implicit none
      double precision dotp
      double precision p1(4),p2(4)
      integer i
c
      dotp = p1(4)*p2(4)
      do i=1,3
         dotp = dotp - p1(i)*p2(i)
      enddo
      return
      end

      
      function dot3p(p1,p2)
      implicit none
      double precision dot3p
      double precision p1(4),p2(4)
      integer i
c
      dot3p = 0.d0
      do i=1,3
         dot3p = dot3p + p1(i)*p2(i)
      enddo
      return
      end

      subroutine cross3p(p1,p2,res)
C     Calculate cross product
      implicit none
      double precision p1(4),p2(4),res(4)
      res(4) = 0.d0
      res(1) = p1(2)*p2(3)-p1(3)*p2(2)
      res(2) = p1(3)*p2(1)-p1(1)*p2(3)
      res(3) = p1(1)*p2(2)-p1(2)*p2(1)
      return
      end

      function spb_kleiss(p1,p2)
      implicit none
      double complex spb_kleiss
      double precision p1(4),p2(4) 
C     Formualae taken from Kleiss
      spb_kleiss = ( ( p1(4) - p1(1) ) * ( p2(4) + p2(1) ) * (1,0)
     #    - ( p1(2) * (1,0) + (0,1) * p1(3) )
     #      *(p2(2) * (1,0) - (0,1) * p2(3)) ) /
     #   dsqrt( (p1(4)-p1(1)) * (p2(4)+p2(1)) )
      end


      subroutine add(s1,s2,s)
      implicit none
      double precision s1(4),s2(4),s(4)
      integer j
      do j=1,4
         s(j) = s1(j) + s2(j)
      enddo
      end

      subroutine boost(p,x,new,ier)
C     This subroutine calculates the boost from the restframe
C     of p into a system where p has the form as given.
C     This boost is applied to x, the result is saved in new.
C     In particular, if x=(0,0,0,mq) the result is new=p!
      implicit none
      double precision p(4), x(4),new(4)
      double precision m,m2,e, help
      double precision dotp
      integer ier
      integer j
C      write(*,*) 'hallo'
      e = p(4)
      help = 0.d0
      do j=1,3
         help = help + p(j)*x(j)
      enddo   

      m2 = dotp(p,p)

      if (m2.lt.0) then
         write(*,*) 'boost: p.p<0, p=', p 
         ier = -3
         return
      endif

      m = dsqrt(m2)

      if (m/e.gt.0.99999999999d0) then
         do j=1,4
            new(j) = x(j)
         enddo
      else
         do j=1,3
         new(j) = x(j) + ( 1.d0/(1.d0+m/e)* help/e + x(4) )/m*p(j)
         enddo
         new(4) = ( e*x(4) + help ) / m
      endif
      ier = 0
      end

      subroutine boost2r(p,x,new)
C     This subroutine calculates the boost from the frame
C     in which p is given into the restframe of p, this
C     boost is applied to x, the result is stored in new
      implicit none
      double precision p(4), x(4),new(4)
      double precision m,mq,beta, e, help
      integer j
      e = p(4)
      mq = p(4)**2
      help = 0.d0
      do j=1,3
         mq = mq - p(j)**2
         help = help + p(j)*x(j)
      enddo
      m = dsqrt(mq)
      beta = dsqrt(1.d0 - mq/e**2)

      if (beta .le. 1.d-13) then
         do j=1,4
            new(j) = x(j)
         enddo
      else
      new(4) = e/m * ( x(4) - help/e )

      do j=1,3
         new(j) = x(j) - (1.d0 - e/m) * help/beta**2/e**2*p(j) 
     #        - x(4)/m * p(j)
      enddo
      endif
      end
      
      subroutine cspin(q1,q2,s)
      implicit none
      double precision q1(4),q2(4),s(4)
      double precision m, mq
      integer j
      
      mq=(q1(4)+q2(4))**2
      do j=1,3
         mq = mq - (q1(j)+q2(j))**2
      enddo
      m = dsqrt(mq)
      do j=1,4
         s(j) = (q1(j) - q2(j))/m
      enddo
      end

      subroutine parity(old,new)
C     perform a parity transformation
      implicit none
      double precision old(4), new(4)
      integer j
      new(4) = old(4)
      do j=1,3
         new(j) = -1.d0 *  old(j)
      enddo
      end

      subroutine cp(source,dest)
      implicit none
      double precision source(4), dest(4)
      integer j
      do j=1,4
         dest(j) = source(j)
      enddo
      end

      function norm3(x)
      implicit none
      double precision norm3
      double precision x(4) 
      norm3 = dsqrt(x(1)**2+x(2)**2+x(3)**2)
      end

      subroutine rotate(R,p)
      implicit none
      double precision R(3,3),p(4),help(4)
      integer mu,i
      do mu = 1,4
         help(mu) = p(mu)
      enddo
      do i=1,3
         p(i) = R(i,1) * help(1)
     .        + R(i,2) * help(2)
     .        + R(i,3) * help(3)
      enddo
      end












